library(ConR)

dummy_ex <- 
  data.frame(ddlat = rnorm(10)*10, 
             ddlon = rnorm(10)*10, taxa = rep("taxa", 10))

context("Test that subpop.comp outputs are correct length and objects")

test_that("subpop.comp", {
  
  SUB <- subpop.comp(dummy_ex, Resol_sub_pop=25)
  
  expect_equal(length(SUB), 2)
  expect_equal(length(SUB[[1]]), 1)

})