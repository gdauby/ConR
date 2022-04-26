
library(ConR)
data(dataset.ex)
dataset.ex <- dataset.ex[dataset.ex$tax == "species_1",]

context("Test that locations.comp outputs are correct length and objects")

test_that("locations.comp", {
  
  locations <- locations.comp(XY = dataset.ex)
  
  testthat::expect_equal(class(locations), "list")
  testthat::expect_output(str(locations$locations_poly), "sf", fixed = TRUE)

})