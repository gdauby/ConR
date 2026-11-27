library(ConR)

context("Test that criterion_C outputs are correct")

data(example_criterionC_subpops)
data(example_criterionC)
assess.yr <- 2000
GL <- 10
modelos <- c("linear", "exponential", "logistic", "general_logistic")


test_that("criterion_C", {
  
    # Example with subpopulations
    result <- suppressWarnings(criterion_C(x = example_criterionC_subpops,
                                           years = NULL,
                                           assess.year = assess.yr,
                                           project.years = NULL,
                                           generation.time = GL,
                                           subpop.size = NULL,
                                           models = modelos,
                                           subcriteria = c("C1", "C2")))

    testthat::expect_equal(class(result), "data.frame")
    testthat::expect_equal(dim(result)[2], 19)
    testthat::expect_equal(result$tax, 
                           paste0("species ", 1:9), fixed = TRUE)
    testthat::expect_equal(result$assessment.year, rep(assess.yr, 9), fixed = TRUE)
    testthat::expect_equal(result$predictive.model, 
                           c(rep("exponential", 5),"logistic","exponential",
                             "linear", "linear"), fixed = TRUE)
    testthat::expect_equal(round(result$assess.pop.size, 0), 
                           c(227,555,49,8000,11000,1245000,229,13900,97))
    testthat::expect_equal(result$any.decline, 
                           c(rep("Decreasing", 6),"Increasing","Decreasing",
                             "Increasing"), fixed = TRUE)
    testthat::expect_equal(result$cont.decline, 
                           c(rep("Decreasing", 6),"Increasing","Decreasing",
                             "Stable"), fixed = TRUE)
    testthat::expect_equal(round(result$reduction_3gen,1), 
                           c(8.5,35.1,63.7,13.0,51.1,94.4,-9.0,39.6,35.3))
    testthat::expect_equal(result$max.subpop.size, 
                           c(100,500,45,6000,9000,1200000,100,7300,82))
    testthat::expect_equal(round(result$prop.subpop.size,1), 
                           c(44.1,90.1,91.8,75.0,81.8,96.4,43.7,52.5,84.5))
    testthat::expect_equal(result$C1, 
                           c("LC or NT","EN","CR","VU",rep("LC or NT", 5)))
    testthat::expect_equal(result$C2, 
                           c("EN","VU","CR","LC or NT",rep("LC or NT", 5)))
    testthat::expect_equal(result$category_C, 
                           c("EN","EN","CR","VU",rep("LC or NT", 5)))
    testthat::expect_equal(result$category_C_code, 
                           c("C2ai","C1","C1+C2ai+C2aii","C1",rep("", 5)))


    # Same example, but using the argument `prop.mature`
    result <- suppressWarnings(criterion_C(x = example_criterionC_subpops,
                                               years = NULL,
                                               assess.year = assess.yr,
                                               project.years = NULL,
                                               generation.time = GL,
                                               prop.mature = 0.85,
                                               subpop.size = NULL,
                                               models = modelos,
                                               subcriteria = c("C1", "C2")))

    testthat::expect_equal(class(result), "data.frame")
    testthat::expect_equal(dim(result)[2], 19)
    testthat::expect_equal(result$tax, 
                           paste0("species ", 1:9), fixed = TRUE)
    testthat::expect_equal(result$assessment.year, rep(assess.yr, 9), fixed = TRUE)
    testthat::expect_equal(result$predictive.model, 
                           c(rep("exponential", 5),"logistic","exponential",
                             "linear", "linear"), fixed = TRUE)
    testthat::expect_equal(round(result$assess.pop.size, 0), 
                           c(193,472,42,6800,9350,1058250,195,11815,82))
    testthat::expect_equal(result$any.decline, 
                           c(rep("Decreasing", 6),"Increasing","Decreasing",
                             "Increasing"), fixed = TRUE)
    testthat::expect_equal(result$cont.decline, 
                           c(rep("Decreasing", 6),"Increasing","Decreasing",
                             "Stable"), fixed = TRUE)
    testthat::expect_equal(round(result$reduction_3gen,1), 
                           c(8.5,35.1,63.7,13.0,51.1,94.4,-9.0,39.6,35.3))
    testthat::expect_equal(result$max.subpop.size, 
                           c(100,500,45,6000,9000,1200000,100,7300,82))
    testthat::expect_equal(round(result$prop.subpop.size,1), 
                           c(44.1,90.1,91.8,75.0,81.8,96.4,43.7,52.5,84.5))
    testthat::expect_equal(result$C1, 
                           c("LC or NT","EN","CR","VU","VU",rep("LC or NT", 4)))
    testthat::expect_equal(result$C2, 
                           c("EN","VU","CR","LC or NT",rep("LC or NT", 5)))
    testthat::expect_equal(result$category_C, 
                           c("EN","EN","CR","VU","VU",rep("LC or NT", 4)))
    testthat::expect_equal(result$category_C_code, 
                           c("C2ai","C1","C1+C2ai+C2aii","C1","C1",rep("", 4)))

    # Example without subpopulations (cannot assess subcriteria C2)
    result <- suppressWarnings(criterion_C(x = example_criterionC,
                                           years = NULL,
                                           assess.year = assess.yr,
                                           project.years = NULL,
                                           generation.time = GL,
                                           subpop.size = NULL,
                                           models = c("quadratic", modelos),
                                           subcriteria = c("C1")))
    
    testthat::expect_equal(class(result), "data.frame")
    testthat::expect_equal(dim(result)[2], 14)
    testthat::expect_equal(result$tax, 
                           paste0("species", 1:9), fixed = TRUE)
    testthat::expect_equal(result$assessment.year, rep(assess.yr, 9), fixed = TRUE)
    testthat::expect_equal(result$predictive.model, 
                           c("quadratic",rep("exponential", 4),"logistic",
                             "quadratic", "linear","exponential"), fixed = TRUE)
    testthat::expect_equal(round(result$assess.pop.size, 0), 
                           c(227,555,49,8000,11000,1245000,229,13900,97))
    testthat::expect_equal(result$any.decline, 
                           c(rep("Decreasing", 6),"Increasing","Decreasing",
                             "Increasing"), fixed = TRUE)
    testthat::expect_equal(result$cont.decline, 
                           # species 1 and 7 are the quadratic fits. Their
                           # fitted parabola has its vertex inside the
                           # assessment window, so it is monotone in neither
                           # direction and no trend is significant.
                           c("Stable", rep("Decreasing", 5), "Stable",
                             "Decreasing", "Stable"), fixed = TRUE)
    testthat::expect_equal(round(result$reduction_3gen,1), 
                           c(8.5,35.1,63.7,13.0,51.1,94.4,-9.0,39.6,35.3))
    testthat::expect_equal(result$C1, 
                           c("LC or NT","EN","CR","VU",rep("LC or NT", 5)))
    testthat::expect_equal(result$category_C, 
                           c("LC or NT","EN","CR","VU",rep("LC or NT", 5)))
    testthat::expect_equal(result$category_C_code, 
                           c("","C1","C1","C1",rep("", 5)))
  
})

test_that("cat_criterion_c assigns VU under C2a(ii)", {

  C2_df <- data.frame(tax = c("sp1", "sp2", "sp3", "sp4"),
                      assess.pop.size = c(5000, 5000, 1000, 5000),
                      any.decline = c("Decreasing", "Decreasing", "Decreasing", "Increasing"),
                      max.subpop.size = c(5000, 3000, 1000, 5000),
                      prop.subpop.size = c(100, 60, 100, 100),
                      mean.fluctuation = 0)

  res <- suppressWarnings(cat_criterion_c(C2_df = C2_df))

  testthat::expect_equal(res$ranks_C, c("VU", "LC or NT", "EN", "LC or NT"))
  testthat::expect_equal(res$cats_code, c("C2aii", "", "C2aii", ""))
})

test_that("criterion_C works without generation.time for several taxa", {

  args <- list(x = example_criterionC, years = NULL, assess.year = assess.yr,
               project.years = NULL, subpop.size = NULL,
               models = c("linear", "exponential"), subcriteria = "C1",
               show_progress = FALSE)

  testthat::expect_error(
    res <- suppressWarnings(do.call(criterion_C, args)), NA)
  testthat::expect_equal(nrow(res), 9)

  # a generation length below 3 years falls back on the same 3, 5 and 10
  # years used when no generation length is given
  res1 <- suppressWarnings(do.call(criterion_C, c(args, generation.time = 1)))
  testthat::expect_equal(res, res1)
})
