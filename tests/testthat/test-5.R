library(ConR)

context("Test that pop.decline.fit outputs are correct")

pop = c(10000, 9100, 8200, 7500, 7000)
yrs = c(1970, 1975, 1980, 1985, 1990)
dataset.ex <- cbind.data.frame(pop.size = pop, years = yrs)
modelos <- c("linear","exponential","quadratic")
yrs.prj <- c(1960, 2005)


test_that("pop.decline.fit", {
  
  result <- pop.decline.fit(dataset.ex, 
                            models = modelos,
                            plot.fit = FALSE,
                            project.years = yrs.prj)
  
  testthat::expect_equal(class(result), "list")
  testthat::expect_output(str(result$best.model), "nlsModel", fixed = TRUE)

  testthat::expect_equal(attributes(result$best.model)$best.model.name, 
                         "exponential", fixed = TRUE)
  
  testthat::expect_output(str(result$model.selection.result), 
                          "data.frame", fixed = TRUE)
  testthat::expect_equal(row.names(result$model.selection.result),
                         modelos, fixed = TRUE)

  testthat::expect_equal(result$predictions$years,
                         sort(c(yrs, yrs.prj)), fixed = TRUE)
  
  testthat::expect_equal(round(result$predictions$predicted, 0),
                         c(11967,9960,9086,8290,7562,6899,5238), fixed = TRUE)
  
  testthat::expect_error(pop.decline.fit(dataset.ex[1:2, ]))

  testthat::expect_error(pop.decline.fit(dataset.ex[ ,1]))
  
    
})
