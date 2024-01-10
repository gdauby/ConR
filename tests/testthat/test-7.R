library(ConR)

context("Test that criterion_A outputs are correct")

data(example_criterionA)
assess.yr <- 2000
GL <- 10
explo <- 5

test_that("criterion_A", {
  
  # one species, two observations in time, one subcriterion
  pop <- c("1970" = 10000, "2000" = 6000)
  result <- suppressWarnings(criterion_A(x = pop, 
                                         years = as.double(names(pop)), 
                                         assess.year = assess.yr, 
                                         subcriteria = c("A2"), 
                                         generation.time = GL))
  
  testthat::expect_equal(class(result), "data.frame")
  testthat::expect_equal(dim(result)[2], 7)
  testthat::expect_equal(result$tax, "species 1", fixed = TRUE)
  testthat::expect_equal(result$assessment.year, assess.yr, fixed = TRUE)
  testthat::expect_equal(result$reduction_A12, 40, fixed = TRUE)
  testthat::expect_equal(result$A2, "VU", fixed = TRUE)
  testthat::expect_equal(result$category_A_code, "A2", fixed = TRUE)
  

  # one species, more observations and subcriteria
  pop <- c("1970"= 10000, "1980"= 8900, "1990"= 7000, "2000"= 6000, "2030"= 4000)
  result <- criterion_A(x = pop, 
                        years = as.double(names(pop)), 
                        assess.year = assess.yr, 
                        project.years = c(2010, 2020, 2030),
                        generation.time = GL)
  
  testthat::expect_equal(class(result), "data.frame")
  testthat::expect_equal(dim(result)[2], 14)
  testthat::expect_equal(result$tax, "species 1", fixed = TRUE)
  testthat::expect_equal(result$predictive.model, "exponential", fixed = TRUE)
  testthat::expect_equal(result$assessment.year, assess.yr, fixed = TRUE)
  testthat::expect_equal(result$reduction_A12, 40, fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A3, 1), 33.3, fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A4, 1), 40.9, fixed = TRUE)
  testthat::expect_equal(result$category_A, "VU", fixed = TRUE)
  testthat::expect_equal(result$category_A_code, "A2+A3+A4", fixed = TRUE)
  testthat::expect_equal(result$A1, "LC or NT", fixed = TRUE)
  testthat::expect_equal(result$A2, "VU", fixed = TRUE)
  testthat::expect_equal(result$A3, "VU", fixed = TRUE)
  testthat::expect_equal(result$A4, "VU", fixed = TRUE)
  

  # Another example: subcriterion A2 and exploitation (A2d)
   pop <- c("1980"= 9000, "1985"= 7500, "1990"= 6000)
   result <- suppressWarnings(criterion_A(x = pop,
                                          years = as.double(names(pop)),
                                          assess.year = assess.yr,
                                          subcriteria = c("A2"),
                                          generation.time = GL,
                                          exploitation = explo))
   
   testthat::expect_equal(class(result), "data.frame")
   testthat::expect_equal(dim(result)[2], 9)
   testthat::expect_equal(result$tax, "species 1", fixed = TRUE)
   testthat::expect_equal(result$predictive.model, "linear", fixed = TRUE)
   testthat::expect_equal(result$assessment.year, 1990, fixed = TRUE)
   testthat::expect_equal(result$reduction_A12, 65, fixed = TRUE)
   testthat::expect_equal(result$A2, "EN", fixed = TRUE)
   testthat::expect_equal(result$category_A_code, "A2", fixed = TRUE)
   testthat::expect_equal("basis_d" %in% names(result), TRUE)
   testthat::expect_equal(result$basis_d, 
                          "Extra reduction: 5%; no change in ranking")
   
  # IUCN (2019) example for criterion A (as in the guideliness)
  result <- suppressWarnings(criterion_A(example_criterionA,
                                         years = seq(1970, 2000, by = 2),
                                         assess.year = assess.yr,
                                         project.years = seq(2002, 2030, by = 2),
                                         generation.time = GL))
  
  testthat::expect_equal(class(result), "data.frame")
  testthat::expect_equal(dim(result)[2], 13)
  testthat::expect_equal(result$tax, 
                         paste0("species ", 1:6), fixed = TRUE)
  testthat::expect_equal(result$assessment.year, rep(assess.yr, 6), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A12, 0), 
                         c(33,30,30,34,26,29), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A3, 0), 
                         c(40,50,44,34,57,10), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A4, 0), 
                         c(58,44,42,33,53,31), fixed = TRUE)
  testthat::expect_equal(result$category_A, 
                         c("EN","EN","VU","VU","EN","VU"), fixed = TRUE)
  testthat::expect_equal(result$category_A_code, 
                         c("A4","A3","A3+A4","A2+A3+A4","A3+A4","A4"))
  testthat::expect_equal(result$A1, rep("LC or NT", 6))
  testthat::expect_equal(result$A2, 
                         c("VU","VU","LC or NT","VU","LC or NT","LC or NT"))
  testthat::expect_equal(result$A3, 
                         c("VU","EN","VU","VU","EN","LC or NT"))
  testthat::expect_equal(result$A4, 
                         c("EN","VU","VU","VU","EN","VU"))
  

  # Same IUCN data and options but assuming different generation lengths
  result <- suppressWarnings(criterion_A(example_criterionA,
                                         years = seq(1970, 2000, by = 2),
                                         assess.year = assess.yr,
                                         project.years = seq(2002, 2030, by = 2),
                                         generation.time = c(2,5,10,15,30,50)))
  testthat::expect_equal(class(result), "data.frame")
  testthat::expect_equal(dim(result)[2], 14)
  testthat::expect_equal(result$tax, 
                         paste0("species ", 1:6), fixed = TRUE)
  testthat::expect_equal(result$assessment.year, rep(assess.yr, 6), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A12, 0), 
                         c(26,9,30,45,26,29), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A3, 0), 
                         c(31,-13,44,45,10,1), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A4, 0), 
                         c(33,8,42,46,34,30), fixed = TRUE)
  testthat::expect_equal(result$category_A, 
                         c("VU","LC or NT","VU","VU","VU","LC or NT"), fixed = TRUE)
  testthat::expect_equal(result$category_A_code, 
                         c("A3+A4","A1+A2+A3+A4","A3+A4","A2+A3+A4","A4","A1+A2+A3+A4"))
  testthat::expect_equal(result$A1, 
                         rep("LC or NT",6))
  testthat::expect_equal(result$A2, 
                         c(rep("LC or NT",3),"VU",rep("LC or NT",2)))
  testthat::expect_equal(result$A3, 
                         c("VU","LC or NT","VU","VU","LC or NT","LC or NT"))
  testthat::expect_equal(result$A4, 
                         c("VU","LC or NT","VU","VU","VU","LC or NT"))

  # Same data but with different options
  result <- suppressWarnings(criterion_A(example_criterionA,
                                         years = NULL,
                                         assess.year = assess.yr + 10,
                                         generation.time = GL))
  
  testthat::expect_equal(class(result), "data.frame")
  testthat::expect_equal(dim(result)[2], 14)
  testthat::expect_equal(result$tax, 
                         paste0("species ", 1:6), fixed = TRUE)
  testthat::expect_equal(result$assessment.year, rep(assess.yr + 10, 6), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A12, 0), 
                         c(54,21,32,33,42,30), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A3, 0), 
                         c(14,79,51,33,62,9), fixed = TRUE)
  testthat::expect_equal(round(result$reduction_A4, 0), 
                         c(58,73,48,34,60,31), fixed = TRUE)
  testthat::expect_equal(result$category_A, 
                         c("EN","EN","EN","VU","EN","VU"), fixed = TRUE)
  testthat::expect_equal(result$category_A_code, 
                         c("A2+A4","A3+A4","A3","A2+A3+A4","A3+A4","A2+A4"))
  testthat::expect_equal(result$A1, 
                         c("VU",rep("LC or NT",5)))
  testthat::expect_equal(result$A2, 
                         c("EN","LC or NT",rep("VU",4)))
  testthat::expect_equal(result$A3, 
                         c("LC or NT","EN","EN","VU","EN","LC or NT"))
  testthat::expect_equal(result$A4, 
                         c("EN","EN","VU","VU","EN","VU"))

})
