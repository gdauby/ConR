library(ConR)
library(testthat)

data(Madagascar.protec)
data(Malagasy.amphibian)

context("Test IUCN.eval with protected areas")

test_that("IUCN.eval", {
  
  MyResults <- IUCN.eval(Malagasy.amphibian, protec.areas = Madagascar.protec)
  
  expect_equal(dim(MyResults), c(201, 12))
  expect_output(str(MyResults), "data.frame", fixed = TRUE)
  
  
})



