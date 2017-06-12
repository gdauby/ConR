library(ConR)
library(testthat)

data(dataset.ex)
data(land)

context("Test that IUCN.eval outputs are correct length and objects")

test_that("IUCN.eval", {
  
  Results <- IUCN.eval(dataset.ex, country_map=land)
  
  expect_output(str(Results), "data.frame")
  expect_equal(dim(Results), c(6,9))
  
  expect_equal(404, Results[1,2])
  
})