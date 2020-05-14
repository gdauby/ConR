library(ConR)
library(testthat)

data(dataset.ex)
data(land)

dataset.ex <- dataset.ex[dataset.ex$tax == "Berlinia bruneelii",]
context("Test that locations.comp outputs are correct length and objects")

test_that("locations.comp", {
  
  locations <- locations.comp(dataset.ex)
  
  testthat::expect_output(str(locations), "List of 2")
  testthat::expect_output(str(locations[[2]]), "Named num", fixed = TRUE)

})