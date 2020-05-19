library(ConR)
data(dataset.ex)

dummy_ex <- 
  data.frame(ddlat = rnorm(10)*10, 
             ddlon = rnorm(10)*10, taxa = rep("taxa", 10))

context("Test that IUCN.eval outputs are correct length and objects")

test_that("IUCN.eval", {
  
  Results <- IUCN.eval(dummy_ex)
  
  expect_equal(class(Results), "data.frame")
  expect_equal(dim(Results), c(1,10))
  

})