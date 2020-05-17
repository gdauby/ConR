
library(ConR)
data(dataset.ex)
dataset.ex <- dataset.ex[dataset.ex$tax == "species_1",]

context("Test that locations.comp outputs are correct length and objects")

test_that("locations.comp", {
  
  locations <- locations.comp(dataset.ex, show_progress = FALSE)
  
  testthat::expect_equal(class(locations), "list")
  testthat::expect_output(str(locations[[1]][[1]]), "SpatialPolygonsDataFrame", fixed = TRUE)

})