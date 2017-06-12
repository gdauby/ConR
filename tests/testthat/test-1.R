
library(ConR)
library(testthat)

data(dataset.ex)
data(land)

context("Test that EOO.computing outputs are correct length and objects")

test_that("EOO.computing", {

  EOO <- EOO.computing(dataset.ex)
  
  expect_output(str(EOO), "data.frame")
  expect_equal(2635042, EOO[1,1])
  expect_equal(dim(EOO), c(6,1))
  
  EOO <- EOO.computing(dataset.ex, export_shp = T)
  
  expect_output(str(EOO), "List of 2")
  expect_is(EOO[[1]][[2]], "SpatialPolygons")
  
  EOO <- EOO.computing(dataset.ex, method.range = "alpha.hull")
  
  expect_equal(154205, EOO[1,1])

  
})