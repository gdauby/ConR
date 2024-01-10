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
