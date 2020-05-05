library(ConR)
library(testthat)

data(dataset.ex)
data(land)

context("Test that locations.comp outputs are correct length and objects")

test_that("locations.comp", {
  
  locations <- locations.comp(dataset.ex)
  
  expect_output(str(locations), "List of 2")
  expect_output(str(locations[[2]]), "Named num [1:6]", fixed = TRUE)
  expect_output(str(locations[[1]]), "SpatialPolygons", fixed = TRUE)
  expect_equal(2, as.numeric(locations[[2]][5]))
  
})