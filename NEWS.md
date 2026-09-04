
# ConR 2.2

* Removed the dependency on `FuzzyNumbers.Ext.2`, which had been archived from
  CRAN and prevented `ConR` from being installed from a clean CRAN library. The
  two functions used from it are now provided internally.
* Fixed the test of continuing decline for quadratic population trends in
  `pop.decline.test()`, which was evaluating the wrong curve: the squared model
  coefficient was used in place of the squared year.
* `pop.decline.test()` no longer fails with `object 'test' not found` when the
  fitted parabola reaches its minimum or maximum within the assessment period.
  Such a trend is increasing over part of the period and decreasing over the
  rest, so no direction can be reported as significant. It is now classified
  from the net change over the period as a non-significant decline or increase.
  This changes the estimated continuing decline (and therefore `criterion_C()`
  results) for taxa whose population trend is best described by such a model.


# ConR 2.1.0

*This entry was reconstructed retrospectively from the changes between the `v1.3.0` and
`v2.1.0` tags, which were released without release notes.*

Version 2 is a rewrite of the package. `ConR` 1.x provided a single assessment function
covering IUCN criterion B; version 2 covers criteria A, B, C and D, separates the
computation of population metrics from their conversion into IUCN categories, and replaces
the retired `sp`/`rgdal`/`rgeos` spatial stack with `sf`/`terra`. The package went from 6
to 32 exported functions, and from a single 3,700-line source file to one file per
function.

## Breaking changes

* `IUCN.eval()` is defunct. Use `criterion_B()`, which covers the same criterion and can
  either compute the underlying metrics itself or accept them as pre-computed arguments.
* `map.res()` is defunct; mapping of results in grid cells is not available in version 2.
* Writing results to an Excel file was removed, along with the `writexl` dependency and the
  `write_results` argument.
* Argument names were made consistently lower-case, and the arguments describing threat
  polygons were generalised beyond protected areas:

  | 1.x | 2.x |
  |---|---|
  | `Cell_size_AOO` | `cell_size_AOO` |
  | `Cell_size_locations` | `cell_size_locations` |
  | `Resol_sub_pop` | `resol_sub_pop` |
  | `Rel_cell_size` | `rel_cell_size` |
  | `protec.areas` | `threat_list` |
  | `method_protected_area` | `method_polygons` |
  | `ID_shape_PA` | `id_shape` (default now `"id_orig"`) |

* `EOO.computing()` lost the `Name_Sp`, `buff_width` and `write_results` arguments.
* Spatial outputs (`export_shp = TRUE`) are now `sf` objects instead of `sp` objects.

## New assessment functions

* `criterion_A()`, `criterion_B()`, `criterion_C()` and `criterion_D()` assess the four
  IUCN criteria.
* `cat_criterion_a()`, `cat_criterion_b()` and `cat_criterion_c()` convert population
  metrics into IUCN categories, and are usable on their own.
* `cat_mult_criteria()` derives the consensus category across criteria,
  `near.threatened()` separates Least Concern from Near Threatened, and `cat_downlist()`
  downgrades categories, for example for regional assessments.

## New population and distribution metrics

* `pop.decline()`, `pop.decline.fit()` and `pop.decline.test()` fit and compare models of
  population trend (linear, quadratic, exponential, logistic, generalised logistic and
  piecewise) and test for a continuing decline.
* `pop.fluctuation()` estimates extreme population fluctuations, and `AOO.decline()`
  estimates decline in the area of occupancy.
* `AOH.estimation()` estimates the area of habitat from a habitat map.
* `severe_frag()` assesses severe fragmentation, `subpop.radius()` and
  `subpop.estimation()` support subpopulation estimation, and `EOO.sensitivity()` tests how
  sensitive the extent of occurrence is to individual occurrences.

## Projections and spatial backend

* `sf` and `terra` replace `sp`, `rgdal`, `rgeos` and `geosphere`. This was necessary
  because `rgdal` and `rgeos` were retired, which is also why the package was removed from
  CRAN.
* A `proj_type` argument was added throughout, and `proj_crs()` resolves it. The default is
  the global cylindrical equal-area projection (`"cea"`, ESRI:54034); `"Antarctic"` and
  `"Africa_eac"` are provided as shortcuts, and any EPSG code is accepted.
* `EOO.computing()` gained a `mode` argument to compute the extent of occurrence either on
  the spheroid (`"spheroid"`, the default) or in a projected plane (`"planar"`).

## Other changes

* `coord.check()` centralises validation of occurrence coordinates, including detection of
  records close to the antimeridian and of taxa with too few unique coordinates.
* `dummy_dist()` generates example occurrence data for testing and examples.
* Added a package tutorial covering the full workflow from occurrences to consensus
  categories, available as an article on the package website.
* New example datasets: `example_criterionA`, `example_criterionC`,
  `example_criterionC_subpops`, `example_fluctuation` and `example_tutorial`.
* New dependencies: `terra`, `stars`, `units`, `data.table`, `dplyr`, `stringr`,
  `segmented`, `nls.multstart` and `lifecycle`. `knitr`, `lwgeom` and `rmapshaper` were
  added as suggested packages.

# ConR 1.3.0

* Adapt to changes in sp and rgdal
* Progress bar now visible when running in parallel
* Default map now retrieved from [rnaturalearth R package](https://CRAN.R-project.org/package=rnaturalearth)
* Various bug fixes
 




