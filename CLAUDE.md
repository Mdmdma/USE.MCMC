# USE.MCMC — agent guide

## Project

R package implementing pseudo-absence sampling for species distribution models, with two sampling backends:

- **MCMC** — `R/mcmcSampling.R`, `R/paSamplingMcmc.R`, `R/acceptNextPoint.R`
- **Nearest-neighbor** — `R/paSamplingNn.R`

Density estimation helpers live in `R/gaussianMixureDensityFunction.R` and `R/mclustDensityFunction.R`. User-facing documentation is in `vignettes/` (the MCMC insights vignette is the canonical narrative). Tests are in `tests/testthat/`.

## Sister repositories

Two repos at `../` depend on or describe this package and are part of the same project. Treat them as in-scope when relevant — read freely, and edit them as part of a task when a change here has clear downstream impact. Mention what you touched.

### `../GaussNiche/` — downstream consumer

Applied scripts that build virtual species in PCA-defined environmental space and call USE.MCMC for pseudo-absence sampling. Loaded via `library(USE)` — treats this package as an installed dependency, **not** as sourced files.

Key files:
- `virtualSpecies_fn.R` — wrappers with a modular sampler interface; the main integration surface.
- `1_developing_framework.R` — pipeline development script.
- `2_testing_wrapper_function.R` — wrapper validation.

When you change an exported function signature, default argument, return shape, or sampler behavior in USE.MCMC, check call-sites in these three files and update them in lockstep.

### `../markov-chain-sampler-paper/` — methodology writeup

LaTeX paper (Wiley NJD class). Sections in `text/`:
- `introduction.tex`, `methods.tex`, `results.tex`, `discussion.tex`
- `appendix/appendix.tex`

Bibliography in `references.bib`; figures in `graphics/`; entry point is `main.tex`.

When a code change affects the paper's claims (methodology, defaults, results), make a surgical edit in the relevant `text/*.tex` file. Flag — don't replace — figures in `graphics/` that may now be stale. Don't recompile LaTeX, restructure sections, rewrite prose for style, or modify `references.bib` unless asked.

**Vignette → paper figure pipeline.** Eight figures in `../markov-chain-sampler-paper/graphics/` are direct knitr outputs from `insights-on-MCMC-pseudo-absence-sampling-vignette.Rmd`. Chunk labels:

- `plot-simple-denisty` → naive density figure
- `plot-density-threshold-sweep` → env threshold-quantile sweep
- `plot-pseudo-absence-density` → full pseudo-absence density
- `pseudo-absence-sweep-contours` → species-quantile sweep with contour rings
- `plot-in-geo` → sampled points in geographical space
- `generate-combined-plot` → posterior comparison across sampling methods (last chunk in the vignette; combines 2D and 3D method comparisons via cowplot)
- `gelman-plot` → Gelman-Rubin convergence
- `autocorreltation-plot` → autocorrelation analysis (chunk label is misspelled; preserved on purpose)

Paper figure name == knitr chunk-output name (`<chunk-label>-1.png`). The other vignette, `Insights-on-nearest-neighbor-search.Rmd`, does NOT feed any current paper figure — it documents an alternative method only mentioned in the appendix narrative. After re-knitting the MCMC vignette, run `make figures` from `../markov-chain-sampler-paper/` to copy the refreshed PNGs into the paper repo, then commit there. If you rename a chunk, also update both the matching entry in `../markov-chain-sampler-paper/Makefile` (`MCMC_FIGS`) and the `\includegraphics{graphics/<name>}` call in the relevant `text/*.tex`. `scheme.png` and `trace_plots.png` are hand-curated and stay outside this pipeline.

**Why both Rmds set `self_contained: false`** in their YAML: stock `html_vignette` defaults to `self_contained: TRUE`, which makes pandoc base64-embed images and then **delete** the `<vignette>_files/figure-html/` directory after rendering. That breaks the copy step — `make figures` would find nothing to copy. Keeping `self_contained: false` makes pandoc reference the PNGs externally, so they survive on disk for the Makefile to grab. If you ever revert this flag, the pipeline silently breaks (knit succeeds, `make figures` reports "Not found" for everything). Don't.

## Cross-repo workflow

- Read sister repos freely for context before answering questions about behavior or design.
- Keep edits scoped to what the current task requires — no opportunistic refactors in sister repos.
- After non-trivial R changes here, run `tests/testthat/`. After API-shape changes, also grep `../GaussNiche/` for the symbol.
- Vignettes in `vignettes/` must stay consistent with the current API; update them alongside breaking changes.
