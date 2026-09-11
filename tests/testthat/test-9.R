library(ConR)

context("Test that pop.decline.test outputs are correct")

pop = c(10000, 9600, 9100, 8200, 7500, 7200, 7000)
pop1 = c(10000, 9900, 9800, 9900, 1000, 9700, 9800)
pop2 = c(7000, 7200, 7500, 8200, 9100, 9600, 1000)
yrs = c(1970, 1973, 1975, 1980, 1985, 1987, 1990)
modelos = c("linear", "quadratic", "exponential", "logistic", "general_logistic")

test_that("pop.decline.test", {

  result0 <- pop.decline(pop.size = pop, years = yrs, models = modelos, 
                         by.taxon = TRUE, show_progress = FALSE)
  testthat::expect_equal(pop.decline.test(result0), "signif.decline")

  result0 <- pop.decline(pop.size = pop, years = yrs, models = "logistic", 
                         by.taxon = TRUE, show_progress = FALSE)
  testthat::expect_equal(pop.decline.test(result0), "signif.decline")
  
  
  result0 <- suppressWarnings(pop.decline(pop.size = pop1, years = yrs, 
                                          models = modelos, 
                         by.taxon = TRUE, show_progress = FALSE))
  testthat::expect_equal(pop.decline.test(result0), "non.signif.decline")
  
  result <- pop.decline.test(result0)
  

})

test_that("pop.decline.test handles a non-monotone quadratic fit", {

  # A U-shaped trajectory: the population declines, bottoms out inside the
  # assessment window, then recovers. The fitted parabola therefore has its
  # vertex inside the window and is monotone in neither direction.
  # Regression test: this used to leave the object 'test' unassigned in the
  # quadratic branch of pop.decline.test(), so the function errored with
  # "object 'test' not found".
  pop.u <- c(10000, 8500, 7400, 7000, 7300, 8400, 9900)
  yrs.u <- c(1970, 1975, 1980, 1985, 1990, 1995, 2000)

  fit.u <- suppressWarnings(pop.decline(pop.size = pop.u, years = yrs.u,
                                        models = "quadratic",
                                        by.taxon = TRUE, show_progress = FALSE))

  res.u <- pop.decline.test(fit.u)

  testthat::expect_type(res.u, "character")
  testthat::expect_length(res.u, 1)
  # Not monotone over the window, so no direction can be claimed as significant
  testthat::expect_true(res.u %in% c("non.signif.decline", "non.signif.increase"))
})
