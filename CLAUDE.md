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

### Uniform-sampling contract (`species.model = NULL` / `species.cutoff.threshold = 1`)

`mclustDensityFunction(species.model = NULL)` builds the "uniform environment" target: density **1 inside** the environmental support (env GMM density ≥ `threshold`) and the floor outside, with no presence GMM subtracted. `paSamplingMcmc(species.cutoff.threshold = 1)` uses it to skip the presence-GMM fit and sample the environment uniformly. The load-bearing convention spans R and C++: the R closure returns `1` directly, and the `rcpp_spec` carries `species_cutoff = Inf` as the **sentinel** plus a *dummy* single-component species GMM. `combined_density()` in `src/mcmc_loop.cpp` detects `std::isinf(sp_cutoff)` and returns `1.0` directly — the dummy species arrays are valid-shaped marshalling ballast only and are never evaluated. If you add finite-cutoff validation in `mclustDensityFunction`, or refactor `combined_density`, preserve this `Inf`-sentinel path or `species.cutoff.threshold = 1` silently breaks. The `mcmc_loop_cpp` **signature is unchanged**, so this needed no `compileAttributes()` regen. Note the percentile direction: **higher `species.cutoff.threshold` = weaker exclusion** (it is fed to `quantile(species.densities, …)`; the target is `1 − sp_density/quantile`, so a high quantile excludes only the densest points), and at the limit `= 1` only the single max-density point would be excluded — uniform mode is the continuous endpoint of that, NOT an inversion.

### Gotcha: roxygen export tag and `MIN_COV_CORRECTION`

`R/mcmcSampling.R` defines a top-level constant `MIN_COV_CORRECTION <- 1e-10`. It must stay **above** the `mcmcSampling()` roxygen block: if it lands between the `#' @export` tag and the function definition, roxygen2 binds the export to the constant and silently drops `mcmcSampling` from `NAMESPACE`. Keep the constant above the `#'` block (with its own one-line comment) when editing this file.

## CRAN-readiness (keep it submittable)

**Not submitting yet.** We are keeping the package CRAN-*ready* but are **NOT** submitting to CRAN at this time. Do not run `devtools::release()` / `devtools::submit_cran()`, do not upload a tarball to the CRAN incoming queue, and do not bump the version "for release" — actual submission is a future, human-initiated decision. Your job is only to keep the invariants below green so that submission is a one-step action when the maintainer decides to.

This package targets CRAN; keep `R CMD check --as-cran` clean (no ERROR/WARNING; only the unavoidable "New submission" NOTE). Before committing anything that touches `DESCRIPTION`, `R/`, `man/`, `src/`, or `vignettes/`, preserve these invariants:

- **License is `GPL (>= 2)`** — the package vendors GPL code from the upstream USE package; never relicense (e.g. to Apache).
- **Tarball stays small (< 5 MB).** Never commit data/download caches (e.g. a `geodata`/WorldClim cache under `vignettes/`) or rendered `vignettes/*_files/` dirs — both are `.Rbuildignore`d / `.gitignore`d; don't un-ignore them. The animation GIFs under `vignettes/` are pkgdown-only and `.Rbuildignore`d.
- **No build-time internet in vignettes.** Use the bundled `Worldclim_tmp`; gate any Suggests / data-package use behind `requireNamespace(...)` (see the `rnaturalearthdata` guard in the nearest-neighbor vignette, and `geodata` chunks set `eval = FALSE`).
- **Docs/deps:** every exported function needs `@return` + documented arguments (incl. `...`); declare every package used anywhere — code, examples, **and vignettes** — in `Imports`/`Suggests`; keep R sources ASCII; `man/`+`NAMESPACE` are roxygen-generated, so edit the `#'` source and re-run `devtools::document()`, never the `.Rd`/`NAMESPACE` by hand.
- **Parallelism:** respect `_R_CHECK_LIMIT_CORES_` — never hardcode `parallel::detectCores()`; cap at 2 under check (see `plotDensityLines.R`).
- **Verify before committing:** `sbatch ../GaussNiche/sbatch/submit_cran_full_check.sh` (builds vignettes via the ETH proxy + runs `--as-cran`). The NOTEs/WARNING it leaves are environment-only — missing `qpdf`, proxy-blocked URL/ORCID checks (`CONNECT tunnel failed, 403`), "unable to verify current time", and the container's compile flags — all clear on CRAN.

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

**Vignette → paper figure pipeline.** This vignette (`insights-on-MCMC-pseudo-absence-sampling-vignette.Rmd`, `dev = "cairo_pdf"`) produces, as **vector PDFs**, the Methods and appendix-diagnostic figures that `../markov-chain-sampler-paper/` pulls in via `make figures` (paper figure name == `<chunk-label>-1.pdf`). Chunk labels it owns:

- `plot-simple-denisty` → naive density figure
- `plot-density-threshold-sweep` → env threshold-quantile sweep
- `plot-pseudo-absence-density` → full pseudo-absence density
- `pseudo-absence-sweep-contours` → species-quantile sweep with contour rings
- `plot-methods-overview` → 2×2 compound figure (naive / env-threshold sweep / species-threshold sweep / sampled points), badge-labelled a–d. Panel d (`p.points.in.the.environment`) needs the full MCMC→remap→thin pipeline, so this chunk lives *after* `points-in-the-environment` (`text/methods.tex`, `\label{fig:methods overview}`)
- `plot-in-geo` → sampled points in geographical space
- `plot-pipeline-scheme` → six-panel workflow scheme (Fig. 1; replaces the legacy hand-drawn `scheme.png`)
- `gelman-plot` → Gelman-Rubin convergence (appendix)
- `autocorreltation-plot` → autocorrelation analysis (appendix; chunk label misspelled, preserved on purpose)
- `trace_plots` → trace + posterior of the chains (appendix; now a real chunk producing `trace_plots-1.pdf`, no longer a static PNG)

After re-knitting, run `make figures` from `../markov-chain-sampler-paper/` (keep `self_contained: false`) to copy the refreshed PDFs, then commit there. If you rename a chunk, also update the matching `MCMC_FIGS` entry in the paper `Makefile` and the `\includegraphics{graphics/<name>}` call in `text/*.tex`. The other vignette, `Insights-on-nearest-neighbor-search.Rmd`, does NOT feed any paper figure (appendix narrative only).

**NOT produced by this vignette.** The paper's Results figures (`sampler-comparison-boxplots.pdf` 2-D / `sampler-comparison-boxplots-5d.pdf` 5-D), the appendix per-species PC matrices (`pcmatrix-5d-sp1..4.pdf`), and the geographic grid (`geo-grid-5d.pdf`) all come from **GaussNiche** (`../GaussNiche/run_2d_experiment.R` / `run_5d_experiment.R` via `experiment_plots.R`) and are copied into the paper's `graphics/` by hand — not via this vignette or `make figures`. The former `generate-combined-plot` posterior-comparison figure was REMOVED from the paper (replaced by those GaussNiche boxplots); `scheme.png` is the legacy hand-drawn figure, now unreferenced.

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
