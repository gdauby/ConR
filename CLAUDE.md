# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`ConR` is an R package (GPL-2, roxygen2 7.2.3, no CRAN release since the rgdal/rgeos
retirement) that computes the population metrics behind IUCN Red List criteria A–D and
assigns threat categories, designed to run over hundreds/thousands of taxa at once.
Spatial work is done with `sf`/`terra` (the old sp/rgdal stack is gone).

## Development workflow

**Any development that is not truly minimal must happen on a dedicated branch, never
directly on `master`.** Branch from the current `master` (or `devel`), keep one branch per
concern, and open a PR to merge back — this is how the existing history works (`devel`,
`new_func`, `paral_test`, `proj_eoo_test` all merged via PRs). "Truly minimal" means a
typo, a one-line doc tweak, or a comment; anything that touches function behaviour,
signatures, `DESCRIPTION`/`NAMESPACE`, or adds files gets a branch.

```r
git switch -c fix-<short-topic> master
```

**Before creating a branch, check `git status`. If there is uncommitted work, stop and
tell the user what it is; resolve it together (commit, discard, etc.) before going further.**
Never stash it silently, and never work around it with `git worktree` or by creating extra
directories outside the repo — all work happens in this checkout.

Regenerate `NAMESPACE` and `man/` with `devtools::document()` in the same commit as the
roxygen change, and run `devtools::test()` before opening the PR.

## Commands

Run from the package root with `Rscript -e "..."` (devtools/roxygen2/testthat are not
in `Suggests`; install them if missing).

```r
devtools::load_all(".")                      # load without installing (main dev loop)
devtools::document()                         # regenerate NAMESPACE + man/*.Rd from roxygen
devtools::test()                             # run all tests
devtools::test(filter = "5")                 # run only tests/testthat/test-5.R
testthat::test_file("tests/testthat/test-3.R")
devtools::check()                            # full R CMD check
R CMD INSTALL .
pkgdown::build_site()                        # -> docs/ (deployed by .github/workflows/pkgdown.yaml on push to master)
```

There is no `tests/testthat.R` driver, so `R CMD check` does **not** run the test suite —
use `devtools::test()` explicitly.

Tests are numbered, not named after the function they cover: test-1 EOO.computing,
test-2 subpop.comp, test-3 criterion_B, test-4 locations.comp, test-5 pop.decline.fit,
test-6 pop.decline, test-7 criterion_A, test-8 criterion_C, test-9 pop.decline.test, test-10 AOH.estimation, test-11 criterion_D.
Several use `context()` (deprecated in testthat 3e); the package is not on testthat edition 3.

`git` may fail here with "detected dubious ownership" — fix once with
`git config --global --add safe.directory D:/MonDossierR/ConR_orig`.

## Architecture

Two layers, deliberately separated:

1. **Metric functions** compute biological/geographic indicators from raw data:
   `EOO.computing`, `AOO.computing`, `locations.comp`, `subpop.comp`, `severe_frag`,
   `AOH.estimation`, `AOO.decline`, `pop.decline` (+ `pop.decline.fit`, `.test`),
   `pop.fluctuation`, `EOO.sensitivity`.
2. **`criterion_*` functions** (`criterion_A/B/C/D`) orchestrate: they either accept
   pre-computed metrics or call the metric functions themselves, then delegate the
   threshold logic to the pure **`cat_criterion_*`** functions (`cat_criterion_a/b/c`),
   which take numeric vectors + threshold vectors and return category vectors.

Keep that split when adding features: thresholds and IUCN category logic belong in
`cat_criterion_*`; spatial/statistical computation belongs in the metric functions.

Multi-criteria consolidation happens *outside* the criterion functions, in the order used
by `vignettes/articles/ConR.Rmd`: merge the `criterion_*` data frames by `tax` →
`near.threatened()` per sub-criterion to split LC from NT → `cat_mult_criteria()` for the
consensus category and the main/auxiliary criteria. `cat_downlist()` downgrades categories
(e.g. regional assessments).

### Conventions that cut across the codebase

- **Occurrence input `XY`/`x`**: a `data.frame` whose *positions* matter and names do not —
  col 1 latitude (numeric, decimal degrees), col 2 longitude, col 3 taxon name. Functions
  rename them to `ddlat`, `ddlon`, `tax` internally.
- **`coord.check()`** is the mandatory front door for any occurrence-based function. It
  coerces/validates coordinates, sanitizes taxon names, optionally projects, splits by
  taxon, and returns `list(list_data, issue_close_to_anti, issue_long_span, issue_nrow,
  unique_occs)`. The `issue_*` slots feed the warning columns of the final results.
- **Projections** go through `proj_crs()`, never raw `st_crs()` strings. `proj_type = "cea"`
  (ESRI:54034, the default) plus shortcuts `"Antarctic"`, `"Africa_eac"`; anything else is an
  EPSG code or string passed through. EOO additionally has `mode = "spheroid"` (default,
  geodetic) vs `"planar"` (projects first).
- **Parallelism** is the same three-line idiom everywhere — copy it rather than inventing
  a new one:
  ```r
  cl <- activate_parallel(parallel, NbeCores)   # PSOCK + doSNOW
  `%d%` <- c_par(parallel)                      # %dopar% or %do%
  pro_res <- display_progress_bar(show_progress, max_pb = length(list_data))
  # foreach(..., .options.snow = pro_res$opts) %d% { ... }
  if (parallel) parallel::stopCluster(cl)
  ```
  Every user-facing loop takes `parallel`, `NbeCores`, `show_progress`.
- **Optional dependencies** are `Suggests`, guarded at the top of the function with
  `requireNamespace(..., quietly = TRUE)` + `stop("The 'x' package is required...")`.
  `lwgeom` gates `EOO.computing()`/`criterion_B()`; `alphahull` gates alpha-hull ranges.
  Tests are written to expect that error and `skip()` when the package is absent.
- Non-standard-evaluation symbols used inside `foreach`/`data.table` must be declared in
  `R/globals.R`.
- Removed functionality is kept as a `.Defunct()` stub (`IUCN.eval`, `map.res`) rather than
  deleted.

### Files worth knowing

- `R/cell.occupied.R` — grid-overlay engine shared by AOO and locations; repeated random
  grid translations keep the minimum number of occupied cells (`nbe_rep`).
- `R/locations.comp.R` — two methods: `"fixed_grid"` (`cell_size_locations`, km) and
  `"sliding_scale"` (`rel_cell_size` × max pairwise distance, Rivers et al. 2010); also
  handles threat polygons via `threat_list`/`threat_weight`.
- `R/pop.decline.fit.R` + `R/ICtab.mod.select.R` — nls.multstart/segmented model fitting
  (linear, quadratic, exponential, logistic, general_logistic, piecewise) and IC-based
  model selection feeding criteria A and C.
- `R/EOO.sensitivity.R` also defines `over.valid.poly()`; `R/over.valid.poly.R` is empty.

`man/*.Rd` and `NAMESPACE` are generated — edit the roxygen blocks and run
`devtools::document()`. `data/` holds the example datasets documented in `R/dataset.ex.R`
(`dataset.ex`, `example_criterionA/C`, `example_tutorial`, `land`, `Madagascar.protec`,
`Malagasy.amphibian`); `data-raw/` is gitignored. `dummy_dist()` generates throwaway
multi-species occurrence data for tests and examples.

The long-form tutorial is `vignettes/articles/ConR.Rmd` (built PDF alongside it) — it is
the reference for the intended end-to-end workflow.
