library(ConR)
library(testthat)

data(dataset.ex)

context("Test IUCN.eval")

test_that("IUCN.eval", {
  
  MyResults <- IUCN.eval(dataset.ex)
  
  testthat::expect_equal(dim(MyResults), c(6, 10))
  testthat::expect_output(str(MyResults), "data.frame", fixed = TRUE)
  
  
})



