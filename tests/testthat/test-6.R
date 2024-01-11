library(ConR)

context("Test that pop.decline outputs are correct")

pop = c(10000, 9050, 8250, 7500, 7200, 6950)
pop1 = c(10000, NA, 8200, NA, NA, 6000)
yrs = c(1970, 1975, 1980, 1985, 1990, 2000)
tax = c("species A", "species B")
dataset.ex = matrix(c(pop, pop1), nrow = length(tax),
  dimnames = list(tax, yrs), byrow = TRUE)
modelos <- c("linear","exponential","quadratic")
yrs.prj <- c(1960, 2050)



test_that("pop.decline", {
  
  # only one species, all models (default)
  result <- suppressWarnings(pop.decline(pop, yrs))
    
  testthat::expect_equal(class(result), "list")
  testthat::expect_output(str(result$predictions$species1), "data.frame", fixed = TRUE)
  testthat::expect_equal(result$predictions$species1$years, yrs)
  testthat::expect_equal(result$predictions$species1$pop.size, pop)
  testthat::expect_equal(round(result$predictions$species1$predicted, 0),
                         c(9670, 9042, 8456, 7907, 7394, 6465))

  
  # two species or more, exponential models with projections
  result <- pop.decline(dataset.ex, models = modelos, project.years = yrs.prj)

  testthat::expect_equal(class(result), "list")
  testthat::expect_equivalent(length(result), 4)
  testthat::expect_output(str(result$predictions$`species A`), "data.frame", fixed = TRUE)
  testthat::expect_output(str(result$predictions$`species B`), "data.frame", fixed = TRUE)
  testthat::expect_equal(result$predictions$`species A`$years,
                         sort(c(yrs, yrs.prj)), fixed = TRUE)
  testthat::expect_equal(result$predictions$`species B`$pop.size,
                         c(NA, pop1, NA), fixed = TRUE)
  testthat::expect_equal(round(result$predictions$`species B`$predicted, 0),
                         c(11768,9914,9100,8353,7666,7037,5928,2516), fixed = TRUE)
  
  # two species or more, different outputs
  result <- pop.decline(dataset.ex, models = modelos, project.years = yrs.prj, 
                        output = "model.fit")
  testthat::expect_equal(class(result), "list")
  testthat::expect_equal(length(result), 1)
  testthat::expect_equivalent(lengths(result), 2)
  testthat::expect_output(str(result$model.fit$`species A`), "nls", fixed = TRUE)
  testthat::expect_output(str(result$model.fit$`species B`), "nls", fixed = TRUE)
  testthat::expect_equal(attributes(result$model.fit$`species A`)$best.model.name, 
                          "exponential", fixed = TRUE)
  testthat::expect_equal(attributes(result$model.fit$`species B`)$best.model.name, 
                         "exponential", fixed = TRUE)
  
  result <- pop.decline(dataset.ex, models = modelos, project.years = yrs.prj, 
                        output = "model.selection")
  testthat::expect_equal(class(result), "list")
  testthat::expect_equal(length(result), 1)
  testthat::expect_equivalent(lengths(result), 2)
  testthat::expect_output(str(result$model.selection$`species A`), 
                          "data.frame", fixed = TRUE)
  testthat::expect_output(str(result$model.selection$`species B`), 
                          "data.frame", fixed = TRUE)

  # Few observations
  testthat::expect_warning(pop.decline(c(10000, 8200, 6000), 
                                       c(1970, 1985, 2000),
                                       models = "all", project.years = 2030))
  testthat::expect_error(pop.decline(c(10000, 6000), c(1970, 2000)))
  
})
