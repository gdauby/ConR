#' @title Calculating Population Reduction
#'
#' @description Based on the fit of statistical models to population data, the 
#' function estimates the decline on the number of mature individuals across 
#' time, expressed in percentage.
#'
#' @param pop.size a vector, data frame or matrix containing the (estimated)
#'   number of mature individuals of species (i.e. population size). If a data
#'   frame or matrix, rows are the species and columns are the population sizes.
#' @param years a vector containing the years for which the population sizes is
#'   available
#' @param taxa a vector containing the name of the species in ```pop.size```
#' @param models a vector containing the names of the statistical models to be
#'   fitted to species population data
#' @param project.years a vector containing the years for which the number of
#'   mature individuals should be predicted using the best candidate statistical
#'   model
#' @param output a character or vector containing the desired output from the
#'   function. The options are: "predictions", "model.fit", "model.selection" 
#'   and "best.model". By default, the function returns only the precitions.
#' @param by.taxon logical. Should the output list be organized by the selected
#'   output options (i.e. predictions, model.fit, model.selection and
#'   best.model) for all taxa or should it cointain one taxon per taxon with
#'   all selected outputs? Defaults to FALSE (list organized by outputs and not
#'   taxa).
#' @inheritParams activate_parallel
#' @param show_progress logical. Whether progress informations should be
#'   displayed. TRUE by default
#' @param ... other parameters to be passed as arguments for functions 
#' ```pop.decline.fit``` and ```ICtab.mod.select```
#'
#' @return a named list
#'
#' @details
#' By default, the function compares the fit of six statistical models to the
#' population trends, namely: linear, quadratic, exponential, logistic,
#' generalized logistic and piece-wise. But users can use different combinations
#' of those models using the argument ```models```, according to the specificity
#' of their study region or groups of organism. See ```pop.decline.fit``` for
#' more details and assumptions on how those models are fitted to population
#' data and how the candidate models are selected.
#' 
#' @seealso \link[ConR]{pop.decline.fit}
#' 
#' @author Renato A. Ferreira de Lima
#'
#' @references 
#' 
#' IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @examples
#' ## Creating vectors with the population data and time intervals 
#' #(adapted from the IUCN 2019 workbook for Criterion A, available 
#' #at: https://www.iucnredlist.org/resources/criterion-a)
#' 
#' pop = c(10000, 9050, 8250, 7500, 7200, 6950)
#' pop1 = c(10000, NA, 8200, NA, NA, 6000)
#' yrs = c(1970, 1975, 1980, 1985, 1990, 2000)
#' tax = c("species A", "species B")
#' pops = matrix(c(pop, pop1), nrow = length(tax), 
#'   dimnames = list(tax, yrs), byrow = TRUE)
#' 
#' ## Fitting data with different models and settings
#' # only one species, all models (default)
#' pop.decline(pop, yrs)
#' 
#' # two species or more
#' pop.decline(pops)
#' 
#' # two species or more, less models 
#' pop.decline(pops, models = c("linear", "quadratic"))
#' pop.decline(pops, models = "exponential")
#' 
#' # two species or more, exponential models with projections
#' pop.decline(pops, models = "exponential", project.years = c(1960, 2050))
#' pop.decline(pops, models = "exponential", project.years = c(1973, 2005))
#' 
#' # two species or more, different outputs 
#' pop.decline(pops, models = "exponential", output = "model.fit")
#' pop.decline(pops, models = c("linear", "quadratic", "exponential"), 
#'  output = "model.selection")
#' 
#' ## Another examples 
#' # Few observations (warning or no model fit below 3 observations) 
#' pop.decline(pop.size = c(10000, 8200, 6000), years = c(1970, 1985, 2000), 
#'  models = "all", project.years = 2030)
#'
#' \dontrun{
#' # Not enough observations (error)
#' pop.decline(pop.size = c(10000, 6000), years = c(1970, 2000))
#' }
#' 
#' 
#' @importFrom foreach foreach
#' 
#' @export pop.decline
#' 
pop.decline <- function(pop.size = NULL,
                        years = NULL,
                        taxa = NULL,
                        models = "all", 
                        project.years = NULL,
                        output = "all",
                        by.taxon = FALSE,
                        parallel = FALSE,
                        NbeCores = 2,
                        show_progress = TRUE,
                        ...) {

  if(is.null(pop.size))
    stop("Please provide a numeric vector, data frame or matrix of population sizes")

  if(!inherits(pop.size, c("matrix", "numeric", "data.frame")))
    stop("Please provide a numeric vector, data frame or matrix of population sizes")
  
  if(inherits(pop.size, "matrix")) 
    pop.size <- as.data.frame(pop.size)

  if(inherits(pop.size, "numeric")) {
    nomes <- names(pop.size)
    pop.size <- as.data.frame(matrix(pop.size, nrow = 1))

    if(is.null(taxa)) {
      row.names(pop.size) <- NULL
    } else {
      row.names(pop.size) <- taxa[1]
    }
    
    if(is.null(nomes)) {
      colnames(pop.size) <- years
    } else {
      colnames(pop.size) <- nomes
    }
  }
  
  if(is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(pop.size)[grepl("[0-9]", names(pop.size))]))
    
    if(is.null(anos))
      stop("Please provide at least two years with estimates of population sizes")

  } else {
    
    anos <- years
    
  }
  
  if(dim(pop.size)[2] > length(as.numeric(anos)))
    pop.size <- pop.size[,grepl(paste0(years, collapse = "|"), names(pop.size)),
                         drop = FALSE]
  
  years <- anos
  
  if(dim(pop.size)[2] < 3 | length(years) < 3)
    stop("Too few time intervals (years) to fit trends of population reduction")
  
  if(is.null(taxa)) {
    nomes <- row.names(pop.size)
    
    if(length(nomes) == 1 & nomes[1] == "1") 
      nomes <- "species1"
    
    taxa <- nomes
    
  } else {
    taxa <- paste0("species", 1:dim(pop.size)[1])
  }
  
  if(!is.null(output)) {
    if (any(!output %in% c("predictions", "model.fit", 
                       "model.selection", "best.model", "all")))
      stop("Please provide one or more of the following output options: 'predictions', 'model.fit', 'model.selection", "best.model' or 'all'")
  }
  
  proj.yrs <- NULL
  if(!is.null(project.years))
    proj.yrs <- 
      sort(unique(project.years[!project.years %in% colnames(pop.size)]))
  
  pop_data <- split(pop.size, f = taxa)
  
  cl <- activate_parallel(parallel = parallel, NbeCores = NbeCores)
  `%d%` <- c_par(parallel = parallel)
  
  pro_res <- display_progress_bar(show_progress = show_progress,
                                  max_pb = length(pop_data))
  opts <- pro_res$opts
  pb <- pro_res$pb
  
  models.fit <- foreach::foreach(
    x = 1:length(pop_data),
    .options.snow = opts
  ) %d% {
    
    if (!parallel & show_progress) setTxtProgressBar(pb, x)
    
    df.x <- data.frame(pop.size = as.numeric(pop_data[[x]]), 
                       years = as.numeric(years))
    project.x <- years[is.na(df.x$pop.size)]
    if(length(project.x) == 0)
      project.x <- NULL
    
    if(!is.null(proj.yrs))
      project.x <- sort(unique(c(project.x, proj.yrs)))
    
    df.x <- df.x[!is.na(df.x$pop.size), , drop = FALSE]

    pred.sp <- pop.decline.fit(x = df.x, 
                               models = models,
                               project.years = project.x,
                               plot.fit = FALSE
                               ,...)
    pred.sp
  }
  
  if(parallel) parallel::stopCluster(cl)
  if(show_progress) close(pb)
  
  names(models.fit) <- taxa
  
  if (by.taxon) {
    
    if (is.null(output)) {

      res.tax <- list(predictions = NULL)
      result <- vector("list", length(taxa))
      result <- lapply(result, function(x) x <- res.tax)
      
      res <- lapply(models.fit, function(x) x$predictions)
      res <- lapply(res, function(x) x[,c("years", "pop.size", "predicted")])
      for(i in seq_along(result)) 
        result[[i]]$predictions <- res[[i]]

      names(result) <- taxa
      
    } else {

      res.tax <- list(predictions = NULL, model.fit = NULL, model.selection = NULL,
                      best.model = NULL)
      #res.tax <- res.tax[output] 
      result <- vector("list", length(taxa))
      result <- lapply(result, function(x) x <- res.tax)
      
      if ("predictions" %in% output | "all" %in% output) {
        res <- lapply(models.fit, function(x) x$predictions)
        res <- lapply(res, function(x) x[,c("years", "pop.size", "predicted")])
        for(i in seq_along(result)) 
          result[[i]]$predictions <- res[[i]]
      }
   
      if ("model.fit" %in% output | "all" %in% output) {
        fits <- lapply(models.fit, function(x) x$best.model)
        for(i in seq_along(result)) 
          result[[i]]$model.fit <- fits[[i]]
      }
      
      if ("model.selection" %in% output | "all" %in% output) {
        if ("model.selection.result" %in% names(models.fit[[1]])) {
          selection <- lapply(models.fit, function(x) x$model.selection.result)
          for(i in seq_along(result)) 
            result[[i]]$model.selection <- selection[[i]]
        }
      }
      
      if ("best.model" %in% output | "all" %in% output) {
        best <- lapply(models.fit, 
                       function(x) attributes(x$best.model)$best.model.name[1])
        for(i in seq_along(result)) 
          result[[i]]$best.model <- best[[i]]
      }
      
      names(result) <- taxa
    }  
    
  } else {
    
    if (is.null(output)) {
      
      result <- list()
      
      res <- lapply(models.fit, function(x) x$predictions)
      res <- lapply(res, function(x) x[,c("years", "pop.size", "predicted")])
      result$predictions <- res
      
    } else {
      
      result <- list()
      
      if ("predictions" %in% output | "all" %in% output) {
        res <- lapply(models.fit, function(x) x$predictions)
        res <- lapply(res, function(x) x[,c("years", "pop.size", "predicted")])
        result$predictions <- res
      }
      
      if ("model.fit" %in% output | "all" %in% output) {
        fits <- lapply(models.fit, function(x) x$best.model)
        result$model.fit <- fits
      }
      
      if ("model.selection" %in% output | "all" %in% output) {
        if ("model.selection.result" %in% names(models.fit[[1]])) {
          selection <- lapply(models.fit, function(x) x$model.selection.result)
          result$model.selection <- selection
        } else {
          res <- vector("list", length(taxa))
          names(res) <- taxa
          result$model.selection <- res
        }
      }
      
      if ("best.model" %in% output | "all" %in% output) {
        best <- lapply(models.fit, 
                       function(x) attributes(x$best.model)$best.model.name[1])
        result$best.model <- best
      }
    }
  }
  
  return(result)
} 
