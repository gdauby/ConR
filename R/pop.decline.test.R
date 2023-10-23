#' @title Test Estimated Continuing Decline
#'
#' @description Based on statistical models fitted to the population size data,
#'   this function assess if the model parameter estimates suggest a continuing
#'   decline. 
#'
#' @param x the object containing the model fitted to population data. Tipically 
#'   the result from `ConR` function `pop.decline.fit` 
#' @param best.name name of the the model fitted to data. Not used if `x` is the
#'   result from function `pop.decline.fit`
#' @param assess.year the year for which the assessment should be performed
#'   
#' @details The function extracts the confidence interval of the parameters of
#'   the model selected to describe the trend in population size. If the
#'   confidence interval of the parameters suggests a significant decline, the
#'   trend is classified as 'significantly decreasing'. If there is evidence of
#'   a decline, but the evidence is not significant, the trend is classified as
#'   'non-significantly decreasing'. For instance, if the inclination parameter
#'   of the linear model is negative and its confidence interval does not
#'   include zero, the trend is 'significantly decreasing'; if the inclination
#'   is negative and the confidence interval include zero, the trend is
#'   'non-significantly decreasing'.
#'   
#'   For the particular case of the piecewise-regression model, the function
#'   returns the classification of the population trend for each time interval.
#'   
#'   If the number of observations is too small (n. obs. <7) to obtain
#'   confidence intervals for models with more than two parameters (e.g.
#'   genealized logistic model), the assessment is carried out empirically, by
#'   assessing if the model predictions provides values that are sucessively
#'   declining. In this case, the trend is just classified as 'increasing' or
#'   'decreasing', and thus test of 'estimated continuing decline' becomes the
#'   same as the test of 'continuing decline at any rate' (sub-criterion B2).
#'   
#'   All significance tests of population trends assume a confidence level of
#'   0.95 (the default of `stats` function `confint()`). In the particular case
#'   of singular gradients of model fit, the function progressively decreases
#'   the confidence level from 0.95 until 0.75 until it gets estimates lower and
#'   upper confidence intervals.
#' 
#' @author Renato A. Ferreira de Lima
#'
#' @references 
#' IUCN 2019. Guidelines for Using the IUCN Red List Categories and Criteria.
#' Version 14. Standards and Petitions Committee. Downloadable from:
#' http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @importFrom FuzzyNumbers.Ext.2 is.decreasing
#' @importFrom FuzzyNumbers.Ext.2 is.increasing
#' @importFrom segmented slope
#' 
#' @export pop.decline.test
#'
#' @examples
#' ## Creating vectors with the population data and time intervals 
#' #(adapted from the IUCN 2019 workbook for Criterion A, available 
#' #at: https://www.iucnredlist.org/resources/criterion-a)
#' pop = c(10000, 9600, 9100, 8200, 7500, 7200, 7000)
#' yrs = c(1970, 1973, 1975, 1980, 1985, 1987, 1990)
#' 
#' ## Fitting data with different models and setting
#' best.model = pop.decline.fit(pop.size = pop, years = yrs, 
#' models = c("linear","exponential","logistic"))
#' pop.decline.test(x = best.model, assess.year = 1990)
#' 
#' best.model = pop.decline.fit(pop.size = pop, years = yrs, models = c("general_logistic"))
#' pop.decline.test(x = best.model, assess.year = 1990)
#' 
#' best.model = pop.decline.fit(pop.size = pop, years = yrs, models = c("quadratic"))
#' pop.decline.test(x = best.model, assess.year = 1990)
#' 
#' best.model = pop.decline.fit(pop.size = pop, years = yrs, models = c("piecewise"))
#' pop.decline.test(x = best.model, assess.year = 1990)
#' 
pop.decline.test <- function(x, 
                            best.name = NULL,
                            assess.year = NULL) {
  
  if(is.null(best.name)) {
    
    best.name <- attributes(x$best.model)$best.model.name[1]
    
    if(is.null(best.name)) {
      
      stop("Please provide one of the name of the model selected to fit population data")
      
    }
  }
  
  params <- stats::coef(x$best.model)
  ys <- x$data$Year[1:which.min(assess.year - x$data$Year)] - 
          min(x$data$Year[1:which.min(assess.year - x$data$Year)], na.rm = TRUE) 
  CI <- suppressMessages(
          try(stats::confint(x$best.model), TRUE))
  
  if(class(CI)[1] == "try-error") {
    CI <- suppressMessages(
      # try(stats::confint(x$best.model), TRUE))
      try(stats::confint(x$best.model, "b"), TRUE))
    
    if(class(CI)[1] != "try-error")
      CI <- t(as.matrix(CI))
      rownames(CI) <- "b"
  }    
  
  if(class(CI)[1] == "try-error") {
    
    seq.ys <- seq(min(ys),  max(ys), by = 1)
    preds <- stats::predict(x$best.model, newdata = data.frame(ys = seq.ys))
    mod <- stats::lm(I(diff(preds) / head(preds, 1)) ~ 1)
    ci.diff <- suppressMessages(stats::confint(mod))
    
    if(stats::coef(mod) < 0)
      test <- if(ci.diff[1]<0 & ci.diff[2]<0) "decrease" else "not.decreasing"
    
    if(stats::coef(mod) >= 0)
      test <- if(ci.diff[1]>0 & ci.diff[2]>0) "increase" else "not.increasing"
    
    
  } else {
    
    if(best.name %in% c("linear", "exponential", "logistic", "general_logistic")) {
  
      if(any(is.na(CI["b",]))) {
        
        p.values <- rev(seq(0.75, 0.95, 0.025))
        i = 1
        CI1 <- CI
        while(any(is.na(CI1["b",]))) {
          # CI1 <- stats::confint(x$best.model, level = p.values[i])
          CI1 <- 
            try(stats::confint(x$best.model, "b", level = p.values[i]), TRUE)
          
          if(class(CI1)[1] != "try-error") {
            CI1 <- t(as.matrix(CI1))
            rownames(CI1) <- "b"
          }
          i = i + 1
        }
        
        if(is.na(CI["b",][1])) 
          CI["b",][1] <- CI1["b",][1]
        
        if(is.na(CI["b",][2])) 
          CI["b",][2] <- CI1["b",][2]
        
      }
      
      if(params["b"] < 0)
        test <- if(CI["b",][1]<0 & CI["b",][2]<0) "signif.decline" else "non.signif.decline"
      
      if(params["b"] >= 0)
        test <- if(CI["b",][1]>0 & CI["b",][2]>0) "signif.increase" else "non.signif.increase"
      
    }
    
    if(best.name == "quadratic") { 
      
      f <- function(x) params["a"] + params["b"]*x + x*I(params["c"]^2)
      decrease <- FuzzyNumbers.Ext.2::is.decreasing(fun = f, x.bound = range(ys), step = 1)
      increase <- FuzzyNumbers.Ext.2::is.increasing(fun = f, x.bound = range(ys), step = 1)
      
      if(params["b"] < 0 & decrease)
        test <- if(CI["b",][1]<0 & CI["b",][2]<0) "signif.decline" else "non.signif.decline"
      
      if(params["b"] > 0 & increase)
        test <- if(CI["b",][1]>0 & CI["b",][2]>0) "signif.increase" else "non.signif.increase"
      
      # vertex <- (-params["b"]/(2*params["c"]))
      # root1 <- (-params["b"] - sqrt(params["b"]^2 + 4*params["c"]*params["a"]))/(2*params["c"]) 
      # root2 <- (-params["b"] + sqrt(params["b"]^2 + 4*params["c"]*params["a"]))/(2*params["c"]) 
      
    }
    
    if(best.name == "piecewise") { 
      
      ys.groups <- findInterval(ys, CI[,1])
      params.CI <- segmented::slope(x$best.model)[[1]]
      tests <- vector("list", length(unique(ys.groups)))
      periods <- vector("list", length(unique(ys.groups)))
      
      if (any(is.na(params.CI[,4])) | any(is.nan(params.CI[,4]))) {
        
        for(i in 1:length(unique(ys.groups))) {
          grp <- unique(ys.groups)[i]
          periods[[i]] <- paste0(" (",paste0(range(x$data$Year[ys.groups %in% grp]), collapse="-"), ")")
          
          if(params.CI[,1][i] < 0)
            tests[[i]] <- "decrease"
          
          if(params.CI[,1][i] > 0)
            tests[[i]] <- "increase"
        }
      } else {
        
        for(i in 1:length(unique(ys.groups))) {
          grp <- unique(ys.groups)[i]
          periods[[i]] <- paste0(" (",paste0(range(x$data$Year[ys.groups %in% grp]), collapse="-"), ")")
          
          if(params.CI[,1][i] < 0)
            tests[[i]] <- 
              if(params.CI[,4][i]<0 & params.CI[,5][i]<0) "signif.decline" else "non.signif.decline"
          
          if(params.CI[,1][i] > 0)
            tests[[i]] <- 
              if(params.CI[,4][i]>0 & params.CI[,5][i]>0) "signif.increase" else "non.signif.increase"
        }
      }
      
      test <- paste0(
                paste0(do.call(c, tests), do.call(c, periods)), 
                     collapse = "|")
    }
  }
  
  return(test)
} 
