#' @title Internal function
#'
#' @description Test whether a function is monotone over an interval, by
#'   evaluating it on a regular grid and inspecting the sign of the successive
#'   differences. Monotonicity is non-strict: a constant function is both
#'   non-increasing and non-decreasing, and so returns `TRUE` either way.
#'
#' @param fun a function of one numeric argument
#' @param x.bound numeric of length two, the interval over which `fun` is
#'   evaluated
#' @param step numeric, the spacing of the evaluation grid. By default, it is 1
#' @param decreasing logical, whether to test for a non-increasing (`TRUE`, the
#'   default) or a non-decreasing (`FALSE`) function
#'
#' @details Replaces `FuzzyNumbers.Ext.2::is.decreasing()` and
#'   `FuzzyNumbers.Ext.2::is.increasing()`, which were the package's only use of
#'   that dependency. `FuzzyNumbers.Ext.2` was archived from CRAN in 2017, which
#'   made `ConR` uninstallable from a clean CRAN-only library.
#'
#' @return a logical of length one
#'
#' @keywords internal
#' @noRd
is_monotone <- function(fun, x.bound, step = 1, decreasing = TRUE) {

  y <- vapply(seq(x.bound[1], x.bound[2], by = step), fun, numeric(1))

  if (decreasing) all(diff(y) <= 0) else all(diff(y) >= 0)
}
