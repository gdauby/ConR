
library(ConR)



# data(land)

dummy_ex <- 
  data.frame(ddlat = rnorm(10)*10, 
             ddlon = rnorm(10)*10, taxa = rep("taxa", 10))

context("Test that EOO.computing outputs are correct length and objects")

test_that("EOO.computing", {

  EOO <- EOO.computing(dummy_ex, show_progress = FALSE, write_shp = FALSE, write_results = FALSE)
  
  expect_output(str(EOO), "data.frame")
  expect_equal(dim(EOO), c(1,1))
  
  
})