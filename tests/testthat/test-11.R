library(ConR)

context("Test that criterion_D outputs are correct")

AOO <- c(1, 2, 5, 10, 15, 20, 25)
locs <- c(1, 1, 2, 3, 4, 5, 6)

test_that("criterion_D uses Loc.threshold for the number of locations", {

  # AOO below 20 km2 or 5 locations or fewer qualify under D2
  res <- criterion_D(AOO = AOO, n.Locs = locs, subcriteria = "D2",
                     AOO.threshold = 20, Loc.threshold = 5)
  testthat::expect_equal(res$D.AOO, c(rep("VU", 5), "LC or NT", "LC or NT"))
  testthat::expect_equal(res$D2.Loc, c(rep("VU", 6), "LC or NT"))
  testthat::expect_equal(res$category_D, c(rep("VU", 6), "LC or NT"))

  # the number of locations can be used without AOO
  res <- criterion_D(n.Locs = locs, subcriteria = "D2", Loc.threshold = 2)
  testthat::expect_equal(res$category_D, c(rep("VU", 3), rep("LC or NT", 4)))
})
