#' @title Compute AIC/AICc Table and Select Best Model
#'
#' @description Compute the values of AIC and AICc, compare them and return the
#'   model comparison, the model with best support from the data and its name.
#'
#' @return a list with the model selection table and the best model.
#'
#' @param x a named list with the models to be compared
#' @param correct numerical. The value of number of observations divided by the
#'   number of model parameters to apply the correction to the AIC
#' @param cutoff numerical. The threshold value of delta-AIC to select equally
#'   plausible models
#' @param parsimony logical. Should the most parsimonious model be select in
#'   case of two or more equally plausible models? Default to TRUE.
#'
#' @details This function was adapted from the functions `AICtab` and `AICctab`
#'   from package `bbmle` (Bolker 2017), and should provide the same outputs.
#'   The main difference of `ICtab.mod.select` is that it some alternatives for
#'   allowing the comparison of models with few or more degrees of freedom then
#'   the number of observations, which return infinite values of logLik of the
#'   correction term in the AICc comparison. These alternatives are
#'   statistically incorrect, but they allow model ranking for selection within
#'   the needs of package `ConR`. These alternatives are provided because the
#'   assessment of species conservation status often relies on very few
#'   observations (e.g. population size estimates).
#'
#'   For the process of selecting the best model, we followed some basic steps.
#'   By default, if two or more models had delta-AIC smaller than the cutoff
#'   provided, the more parsimonious model (i.e. the model with less parameters)
#'   is selected. However, this decision can be changed (i.e. slect the model
#'   with best fit) by setting the argument 'parsimony' to FALSE. Next, if more
#'   than one model is selected (i.e. both have the same number of parameters),
#'   the selection process give preference to models that are not linear or
#'   quadratic, which tend provide projections that reach zero much faster and
#'   that can generate non-realistic projections depending on the population
#'   data or on the years chosen for the projection period.
#'  
#' @author Renato A. Ferreira de Lima
#'
#' @references 
#'  D. Anderson and K. Burnham (2004). Model selection and multi-model
#'  inference. Second edition. Springer-Verlag, New York.
#'  
#'  Ben Bolker and R Development Core Team (2017). bbmle: Tools for General Maximum Likelihood
#'   Estimation. R package version 1.0.20. https://CRAN.R-project.org/package=bbmle
#' 
ICtab.mod.select <- function(x, correct = 40, cutoff = log(8), parsimony = TRUE) {
  
  x1 <- x[!is.na(x)]
  if(is.null(names(x1))) 
     names(x1) = paste("model ", 1:length(x1), sep= "")
    
  nobs <- sapply(x1, stats::nobs)
  dfs <- sapply(x1, function(x) attr(stats::logLik(x), "df"))
  nks <- nobs/dfs 
  logLiks <- sapply(x1, function(x) c(stats::logLik(x)))
  
  if(any(is.infinite(logLiks)))
    logLiks[is.infinite(logLiks)] = min(logLiks[!is.infinite(logLiks)]) - 2
  
  aics <- (-2*logLiks) + 2*dfs
  corr <- 2*dfs*(dfs + 1)/(nobs - dfs - 1)
  
  if(any(is.infinite(corr))) 
    corr[is.infinite(corr)] <- 2 * dfs[is.infinite(corr)] *(dfs[is.infinite(corr)]) /
        (nobs[is.infinite(corr)] - dfs[is.infinite(corr)])
  
  aiccs <- aics + abs(corr)

  if(any(nks >= correct)) {
    
    dIC <- round(aics - min(aics, na.rm = TRUE), 2)
    tab <- data.frame(df = dfs)
    tab <- cbind(dAIC = dIC, tab)
    ICs <- tab
    
  } else {
    
    dIC <- round(aiccs - min(aiccs, na.rm = TRUE), 2)
    tab <- data.frame(df = dfs)
    tab <- cbind(dAICc = dIC, tab)
    ICs <- tab
    
  }
  
  if(any(ICs[, 1] <= cutoff)) {
    
    if(all(ICs[, 1] <= cutoff)) {
      
      if(parsimony) {
        
        best <- x1[which(ICs$df == min(ICs$df))][[1]]
        best.name <- names(x1[which(ICs$df == min(ICs$df))])
        
      } else {
        
        best <- x1[which(ICs[, 1] == min(ICs[, 1]))][[1]]
        best.name <- names(x1[which(ICs[, 1] == min(ICs[, 1]))])

      }
      
    } else {
      
      id <- which(ICs[, 1] <= cutoff)
      
      if(parsimony) {
        
        id <- id[which(ICs$df[id] == min(ICs$df[id]))]
        
      } else {
        
        id <- id[which(ICs[, 1][id] == min(ICs[, 1][id]))]
        
      }
      
      best.name <- names(x1[id])
      
      if(length(best.name) > 1) {
        
        best.name = best.name[!grepl("linear|quadratic", best.name)]
        warning("A linear or quadratic model was also selected among the candidate models but was not considered as the best model")
      }
      
      best <- x1[best.name][[1]]
      
    }  
  } else {
    
    best <- x1[which(ICs[, 1] == min(ICs[, 1]))][[1]]
    best.name <- names(x1[which(ICs[, 1] == min(ICs[, 1]))])
    
  }
  
  return(list(ICtab = ICs,
              best.model = best,
              best.model.name = best.name))
}
