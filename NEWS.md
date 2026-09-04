
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

# ConR 1.3.0

* Adapt to changes in sp and rgdal
* Progress bar now visible when running in parallel
* Default map now retrieved from [rnaturalearth R package](https://CRAN.R-project.org/package=rnaturalearth)
* Various bug fixes
 




