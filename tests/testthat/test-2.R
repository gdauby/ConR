library(ConR)
library(testthat)

data(dataset.ex)
data(land)

context("Test that subpop.comp outputs are correct length and objects")

test_that("EOO.computing", {
  
  SUB <- subpop.comp(dataset.ex, Resol_sub_pop=25)
  
  expect_output(str(SUB), "List of 6")
  expect_output(str(SUB[[1]]), "List of 2")
  expect_is(SUB[[1]][[2]], "SpatialPolygons")
  expect_equal(1, length(SUB[[1]][[1]]))
  
})