# USE.MCMC - agent guide

## Project

R package implementing pseudo-absence sampling for species distribution models, with two sampling backends:

- **MCMC** - `R/mcmcSampling.R`, `R/paSamplingMcmc.R`, `R/acceptNextPoint.R`. Inner loop has a C++ implementation in `src/mcmc_loop.cpp` (Rcpp + RcppArmadillo) reached via `R/RcppExports.R`. Source installs therefore need a C++ toolchain.
- **Nearest-neighbor** - `R/paSamplingNn.R`, with helpers in `R/maxResNn.R` / `R/optimalDistanceThresholdNn.R`. This backend is **dimension-general** (arbitrary `dimensions`, default `c("PC1","PC2")`), not 2D-only: it draws uniform proposals over the PC bounding box (no `sf` grid), NN-remaps them, and excludes the presence region with a metric density criterion (no convex hull). Two load-bearing facts a "simplification" could break: (1) the NN *remap* is uniform-over-support in every dimension and needs **no** correction, but the support-distance *threshold* does: it is scaled by `.nnThresholdDimFactor` (default `√(d/2)`, ==1 at d=2; see `optimalDistanceThresholdNn`'s `dim.correction`). Keep the factor; dropping it over-rejects the low-density interior in d>2. (2) Box-uniform proposals are correct in all d but their acceptance rate falls ~`v_d/2^d`, so high-d is compute-bound: `n.candidates`/`n.tr` are the knobs, and an adaptive batch loop tops up toward `n.samples`. `../GaussNiche/virtualSpecies_nd_fn.R` exposes this as the `pa_nn` sampler.

Density estimation lives in `R/mclustDensityFunction.R` (mclust-backed GMM, with a vectorized `fast_gmm_density_batch` helper). User-facing documentation is in `vignettes/` (the MCMC insights vignette is the canonical narrative). Tests are in `tests/testthat/`.

### Engine dispatch (MCMC backend)

`mcmcSampling()` and `paSamplingMcmc()` take `engine = c("auto", "R", "cpp")`. `"auto"` (the default) picks the C++ path whenever both the density and proposal closures carry an `rcpp_spec` attribute. The attribute is attached by:

- `mclustDensityFunction()` (`R/mclustDensityFunction.R`) - `type = "mclust_density"` with the precomputed GMM parameters.
- `addHighDimGaussian()` (`R/addHighDimGaussian.R`) - `type = "gaussian_proposal"` with the proposal mean and covariance.

Custom user closures lack the attribute and stay on the R reference loop. `engine = "cpp"` errors when called with a closure that lacks `rcpp_spec`. Burn-in on both engines is single-pass Robbins-Monro adaptation (target acceptance 0.234, γ_t = 1/(t+1)^0.6); `max.burnin.cycles` is accepted but ignored with a soft deprecation warning. If you change the C++ entry point in `src/mcmc_loop.cpp`, regenerate `R/RcppExports.R` and `src/RcppExports.cpp` via `Rcpp::compileAttributes()` (and then `devtools::document()` for the man pages).

### Gotcha: roxygen export tag and `MIN_COV_CORRECTION`

`R/mcmcSampling.R` defines a top-level constant `MIN_COV_CORRECTION <- 1e-10`. It must stay **above** the `mcmcSampling()` roxygen block: if it lands between the `#' @export` tag and the function definition, roxygen2 binds the export to the constant and silently drops `mcmcSampling` from `NAMESPACE`. Keep the constant above the `#'` block (with its own one-line comment) when editing this file.

## Sister repositories

Two repos at `../` depend on or describe this package and are part of the same project. Treat them as in-scope when relevant: read freely, and edit them as part of a task when a change here has clear downstream impact. Mention what you touched.

### `../GaussNiche/` - downstream consumer

Applied scripts that build virtual species in PCA-defined environmental space and call USE.MCMC for pseudo-absence sampling. Loaded via `library(USE)`, treats this package as an installed dependency, **not** as sourced files.

Key files:
- `virtualSpecies_fn.R` - wrappers with a modular sampler interface; the main integration surface.
- `1_developing_framework.R` - pipeline development script.
- `2_testing_wrapper_function.R` - wrapper validation.

When you change an exported function signature, default argument, return shape, or sampler behavior in USE.MCMC, check call-sites in these three files and update them in lockstep. The recent `engine = c("auto", "R", "cpp")` argument on `paSamplingMcmc()` / `mcmcSampling()` defaults to `"auto"`, so existing GaussNiche callers keep working without changes; only mention `engine` to GaussNiche when intentionally forcing the R path for comparison.

### `../markov-chain-sampler-paper/` - methodology writeup

LaTeX paper (Wiley NJD class). Sections in `text/`:
- `introduction.tex`, `methods.tex`, `results.tex`, `discussion.tex`
- `appendix/appendix.tex`

Bibliography in `references.bib`; figures in `graphics/`; entry point is `main.tex`.

When a code change affects the paper's claims (methodology, defaults, results), make a surgical edit in the relevant `text/*.tex` file. Flag (don't replace) figures in `graphics/` that may now be stale. Don't recompile LaTeX, restructure sections, rewrite prose for style, or modify `references.bib` unless asked.

**Vignette → paper figure pipeline.** Ten figures in `../markov-chain-sampler-paper/graphics/` are direct knitr outputs from `insights-on-MCMC-pseudo-absence-sampling-vignette.Rmd`. Chunk labels:

- `plot-simple-denisty` → naive density figure
- `plot-density-threshold-sweep` → env threshold-quantile sweep
- `plot-pseudo-absence-density` → full pseudo-absence density
- `pseudo-absence-sweep-contours` → species-quantile sweep with contour rings
- `plot-methods-overview` → 2×2 compound figure assembling three sampling-function construction stages plus the sampled result (naive / env-threshold sweep / species-threshold sweep / sampled points in environmental space), badge-labelled a–d. Because panel d (`p.points.in.the.environment`) needs the full MCMC→remap→thin pipeline, this chunk lives *after* the `points-in-the-environment` chunk, not next to the density chunks (`text/methods.tex`, `\label{fig:methods overview}`)
- `plot-in-geo` → sampled points in geographical space
- `plot-pipeline-scheme` → six-panel workflow scheme (Fig. 1; replaces the legacy hand-drawn `scheme.png`)
- `generate-combined-plot` → posterior comparison across sampling methods (last chunk in the vignette; combines 2D and 3D method comparisons via cowplot)
- `gelman-plot` → Gelman-Rubin convergence
- `autocorreltation-plot` → autocorrelation analysis (chunk label is misspelled; preserved on purpose)

Paper figure name == knitr chunk-output name (`<chunk-label>-1.png`). The other vignette, `Insights-on-nearest-neighbor-search.Rmd`, does NOT feed any current paper figure: it documents an alternative method only mentioned in the appendix narrative. After re-knitting the MCMC vignette, run `make figures` from `../markov-chain-sampler-paper/` to copy the refreshed PNGs into the paper repo, then commit there. If you rename a chunk, also update both the matching entry in `../markov-chain-sampler-paper/Makefile` (`MCMC_FIGS`) and the `\includegraphics{graphics/<name>}` call in the relevant `text/*.tex`. `trace_plots.png` is hand-curated and stays outside this pipeline. `scheme.png` is the legacy hand-drawn workflow figure; it remains on disk for reference but is no longer referenced from the paper (replaced by `plot-pipeline-scheme-1.png`).

**Why both Rmds set `self_contained: false`** in their YAML: stock `html_vignette` defaults to `self_contained: TRUE`, which makes pandoc base64-embed images and then **delete** the `<vignette>_files/figure-html/` directory after rendering. That breaks the copy step: `make figures` would find nothing to copy. Keeping `self_contained: false` makes pandoc reference the PNGs externally, so they survive on disk for the Makefile to grab. If you ever revert this flag, the pipeline silently breaks (knit succeeds, `make figures` reports "Not found" for everything). Don't.

## Cross-repo workflow

- Read sister repos freely for context before answering questions about behavior or design.
- Keep edits scoped to what the current task requires: no opportunistic refactors in sister repos.
- After non-trivial R changes here, run `tests/testthat/`. After API-shape changes, also grep `../GaussNiche/` for the symbol.
- **Running tests / docs needs the rocker container + a GDAL bind.** The user's R package library (`~/R/rocker-rstudio/4.5`, R 4.5) only loads inside `/cluster/scratch/$USER/rocker_rstudio_4.5.sif`; the Euler module R (4.4.1) can't load it (GLIBC). `terra`/`sf` additionally need `libgdal.so.33`/`libproj.so.25`/`libgeos_c.so.1`/`libudunits2.so.0`, absent in the bare container, so inject a soname-matching set via `LD_LIBRARY_PATH` (GDAL+PROJ from the paraview manual install, GEOS+udunits from the module stack). `tools/check_nn_generalization.sh` is the ready-made SLURM job that does all of this and runs `devtools::test()` + `roxygenise()`; it last reported `PASS=271 FAIL=0`. (This works because the tests build rasters in-memory via `terra::rast(df, type="xyz")` and never touch GDAL drivers.) For terra/sf-free logic, `tools/nn_core_harness.R` validates the FNN-only NN sampling core (dimension-generality, the `√(d/2)` threshold, density flattening, exclusion direction) with no GDAL needed. See the `r-execution-environment` memory for the exact lib paths.
- Vignettes in `vignettes/` must stay consistent with the current API; update them alongside breaking changes.

## Keeping this guide up to date

This file is part of the working surface: treat it like code, not like prose that lives forever. As you work in the repo, keep it accurate:

- If you delete or rename a file path mentioned here (e.g. an `R/*.R` file, a vignette, a chunk label, a sister-repo file), update or remove the corresponding line in the same change.
- If you add a new sampling backend, a new engine knob, a new on-disk contract between repos, or a load-bearing convention (something a future agent could break with a "routine" edit), add a short note here: lead with the rule, then a one-line *why*.
- If you discover a gotcha by hitting it (e.g. the `MIN_COV_CORRECTION` / roxygen export issue), record it here so the next agent doesn't have to rediscover it.
- Conversely, prune entries that have become wrong or stale. Outdated guidance is worse than missing guidance.
- The bar for inclusion is *non-obvious from reading the code*. Don't restate what well-named functions or directory structure already convey.
- Don't add a "Changelog" or "Recent changes" section here: `git log` and `NEWS.md` already cover that. This guide describes the *current* state.

Cross-session preferences and one-off feedback (e.g. "the package is unpublished, skip migration prose") belong in the `memory/` system under `~/.claude/projects/.../memory/`, not here. Use this file for facts that travel with the codebase and apply to every agent that opens it.
