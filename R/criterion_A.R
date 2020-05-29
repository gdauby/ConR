#' @title Assess IUCN Criterion A
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion A, which is based on population size reductions (Criteria A1, 
#'  A2, A3 and A4)
#'
#' @param x a named vector or a data frame containing the population sizes (e.g.
#'   number of mature individuals of the species) from the oldest to the most 
#'   recent estimate.
#' @param years a vector containing the years for which the population sizes
#'   is available (i.e. time series). Can be NULL if x contains the years as names.
#' @param assess.year the year for which the assessment should be performed.
#' @param project.years a vector containing the years for which population sizes
#'   were or should be projected.
#' @param generation.time a value or vector of generation lengths, i.e. the
#'   average age of parents of the current cohort (IUCN 2019).
#' @param models a vector containing the names of the models to be fitted to
#'   species population data to perform projections.
#' @param subcriteria a vector containing the sub-criteria that should be
#'   included in the assessment (i.e. A1, A2, A3 and/or A4).
#' @param data.type a character corresponding to the type of data (IUCN 2019):
#'   "observation", "index" or "AOO_EOO" (only these types are currently
#'   implemented)
#' @param nature.evidence a character corresponding to nature of evidence (IUCN
#'   2019): "observed", "estimated", "projected", "inferred" or "suspected"
#' @param A1.threshold numeric vector with the A1 thresholds to convert decline
#'   estimates into categories. Default is the thresholds recommended by IUCN.
#' @param A234.threshold numeric vector with the A2, A3 and A4 thresholds to
#'   convert decline estimate into categories. Default is the thresholds
#'   recommended by IUCN.
#' @param all.cats logical. Should the categories from all criteria be returned
#'   and not just the consensus categories?
#'
#' @return a data frame containing, for each of the species provided, the year of 
#'   assessment, the time interval of the assessment (include past and future estimates,
#'   if any), the population sizes in the intervalof assessment, the reduction of the 
#'   population size using the chosen sub-criteria (A1, A2, A3 and A4), the model used
#'   to obtain the projections of population size (if used), the IUCN categories 
#'   associated with these sub-criteria and the consensus category for criterion A.
#'
#' @details TO BE COMPLETED
#'
#'   Some important notes. The function can return the predictions of population
#'   estimates for years not in the observed data.
#'   
#'   If `years` is a subset of all the years contained in `x`, them `x` is filtered
#'   based on `years`. So, make sure you have selected the right years.
#'   
#'   If the year of assessment is not given, the most recent year is taken instead.
#'
#' @author Lima, R.A.F. & Dauby, G.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @export criterion_A
#'
#' @examples
#' 
#' ## Simplest example: one species, two observations in time, one subcriteria
#'  pop = c("1970" = 10000, "2000" = 6000)
#'  criterion_A(x = pop,
#'   years = c(1970, 2000), 
#'   assess.year = 2000,
#'   project.years = NULL,
#'   subcriteria = c("A2"),
#'   generation.time = 10)
#'   
#' ## Another example: one species, more observations and subcriteria
#' pop = c("1970" = 10000, "1980" = 8900, "1990" = 7000, "2000" = 6000, "2030" = 4000)
#' criterion_A(x = pop,
#'   years = c(1970, 1980, 1990, 2000, 2030), 
#'   assess.year = 2000,
#'   project.years = 2030,
#'   subcriteria = c("A1", "A2", "A3", "A4"),
#'   generation.time = 10)
#'   
#' ## Another example...
#'  pop = c("1980" = 9000, "1985" = 7500, "1990" = 6000)
#'  criterion_A(x = pop,
#'   years = c(1980, 1985, 1990), 
#'   assess.year = 2000,
#'   project.years = NULL,
#'   subcriteria = c("A2"),
#'   generation.time = 10)
#'
#' ## The data and criterion A assessment as described in IUCN (2019) - https://www.iucnredlist.org/resources/criterion-a
#' data(example_criterionA)
#' criterion_A(example_criterionA,
#'   years = seq(1970, 2000, by = 2), 
#'   assess.year = 2000,
#'   project.years = seq(2002, 2030, by = 2),
#'   subcriteria = c("A1", "A2", "A3", "A4"),
#'   generation.time = 10)
#'
#' ## Same data but with different options (need to use predictions from statistical models)
#' criterion_A(example_criterionA,
#'   years = NULL, 
#'   assess.year = 2010,
#'   project.years = NULL,
#'   subcriteria = c("A1", "A2", "A3", "A4"),
#'   generation.time = 10)
#'   
criterion_A = function(x, 
                       years = NULL, 
                       assess.year = NULL, 
                       project.years = NULL,
                       generation.time = NULL,
                       models = c("linear", "quadratic", "exponential", "logistic", "general_logistic","piecewise"),
                       subcriteria = c("A1", "A2", "A3", "A4"),
                       data.type = NULL,
                       nature.evidence = NULL,   
                       A1.threshold = c(50, 70, 90),
                       A234.threshold = c(30, 50, 80),
                       all.cats = TRUE) {
  
  if(is.null(x))
    stop("Please provide at least two estimates of population sizes")
  
  if(!any(subcriteria %in% c("A1", "A2", "A3", "A4")))
    stop("Please provide at least one sub-criteria for the assessment: A1, A2, A3 and/or A4")
  
  if(is.vector(x)) {
    
    x = as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, names(x))),
                      stringsAsFactors = FALSE)

  }
  
  if(is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    
    if(is.null(anos)) 
      stop("Please provide at least two years with estimates of population sizes") 
    
    years <- anos
    warning("The years of the population sizes were not given and were taken from the input population data", call. = FALSE)
  }
  
  if(length(years) < 2)
    stop("At least two years are needed to perform the assessment")
  
  if(!is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    all.yrs <- years
    if(!is.null(project.years)) 
      all.yrs = unique(c(all.yrs, project.years))
    
    if(!is.null(anos) & any(!anos %in% all.yrs)) {
      
      if(class(x[,1]) %in% c("factor", "character")) {
        
        x <- cbind.data.frame(x[1], x[,-1][ , anos %in% all.yrs])
        
      } else {
        
        x <- x[ , anos %in% all.yrs]
        
      }
    }
  } 
  
  if(is.null(assess.year)) {
    
    assess.year <- years[which.min(abs(years - as.numeric(format(Sys.Date(), "%Y"))))]
    warning("Year of assessment not given: assuming to be the most recent year of the assessment period")
    
  }
  
  if(!assess.year %in% years) {
    
    assess.year <- years[which.min(abs(years - assess.year))]
    warning(paste0("Year of assessment not in the provided time series: assuming the closest year: ",
            assess.year))
    
  }
  
  if(is.null(generation.time)) {
    
    prev.year <- assess.year - 10
    proj.year <- assess.year + 10
    warning("Generation length not given: assuming 10 years. Please, check if this is accurate for your species")
    
  } else {
    
    if((3 * generation.time) >= 10) {
      
      prev.year <- assess.year - 3 * generation.time
      proj.year <- assess.year + 3 * generation.time
      
      if((proj.year - assess.year)>100) {
        
        proj.year <- assess.year + 100
        warning("Maximum year to project population sizes is more than 100 years into the future: assuming 100 years after the year of assessment")
        
      }
      
    } else {
      
      prev.year = assess.year - 10
      proj.year = assess.year + 10
      warning("Three times the generation length is smaller than 10 years: assuming 10 years")
      
    }
  }
  
  if(any(subcriteria %in% c("A1", "A2")) & !any(subcriteria %in% c("A3", "A4"))) proj.year = assess.year
  
  if(!any(subcriteria %in% c("A1", "A2")) & any(subcriteria %in% c("A3", "A4"))) prev.year = assess.year
  
  if(is.null(project.years)) {
    
    all.yrs <- prev.year:proj.year
    yrs <- all.yrs[all.yrs %in% years]
    int = median(diff(years), na.rm=TRUE)
    if (!prev.year %in% yrs)
      yrs <- unique(c(seq(prev.year, min(yrs), by= int), yrs))
    if (!proj.year %in% yrs)
      yrs <- unique(c(yrs, seq(max(yrs), proj.year, by= int)))

  } else {
    
    yrs <- unique(c(years, project.years))
    int = median(diff(years), na.rm=TRUE)
    if (prev.year < min(yrs))
      yrs <- unique(c(seq(prev.year, min(yrs), by= int), yrs))
    if (proj.year > max(yrs))
      yrs <- unique(c(yrs, seq(max(yrs), proj.year, by= int)))
    
    if(length(project.years) == 1 & (max(project.years) - max(years[!years %in% project.years]) >= int))
      yrs = c(yrs[1:which.max(years[!years %in% project.years])],
              seq(max(years[!years %in% project.years]) + int, project.years - int, by = int),
              project.years)
  }  
  
  ## if the first column indicate the taxon name
  if(class(x[,1]) %in% c("factor", "character")) {
    
    names(x)[-1] <- gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    pop_data <- split(x[ ,grepl("[0-9]", names(x))], f = x[,1])
    
  } else {
    
    names(x) <- gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    nomes <- paste("species", 1:dim(x)[1])
    pop_data <- split(x, f = nomes)
    
  }
  
  best.models = NULL
  
  if(any(!yrs %in% names(pop_data[[1]]))) {   # Predictions based on population trends
    
    if(length(x) < 3) 
      stop("Too few year intervals to fit a model to population trends")
    
    best.models = pop_data
    
    ## Renato: Gilles, il faut peut-etre mettre ici la boucle en dplyr et/ou paralell
    for(i in 1:length(pop_data)) {
      
      pred.sp <- pop.decline.fit(pop.size = pop_data[[i]], 
                                 years = years, 
                                 models = models,
                                 project.years = yrs,
                                 plot.fit = FALSE,
                                 max.count = 50)
      pred.pop.size <- pred.sp$data$Observed
      pred.pop.size[is.na(pred.sp$data$Observed)] <- 
        pred.sp$data$Predicted[is.na(pred.sp$data$Observed)]
      names(pred.pop.size) <- pred.sp$data$Year
      pop_data[[i]] <- pred.pop.size
      best.models[[i]] <-  attributes(pred.sp$best.model)$best.model.name
      
    }
  }

  assess.period = paste(unique(sort(c(range(yrs), assess.year))), collapse="-")
  pop_data1 = lapply(pop_data, function(y) y[names(y) %in% yrs])
  ps.interval = sapply(pop_data1, function(y) paste(
                      unique(c(as.character(y[which(yrs %in% min(yrs))]),
                        as.character(y[which(yrs %in% assess.year)]),
                        as.character(y[which(yrs %in% max(yrs))]))), collapse = "-"))
 
  ## Population reduction using IUCN criteria
  Results = data.frame(
    species = names(pop_data),
    assessment.year = assess.year,
    assessment.period = assess.period,
    assessment.pop.sizes = ps.interval,
    #reduction_A12 = 100 * as.numeric(unlist(reduction_A12)),
    #reduction_A3 = 100 * as.numeric(unlist(reduction_A3)),
    #reduction_A4 = 100 * as.numeric(unlist(reduction_A4)),
    # category_cA = NA,
    # category_cA_code = NA,
    stringsAsFactors = FALSE
  )
  row.names(Results) = NULL
  
  if(!is.null(best.models))
    Results$predictive.model = as.character(unlist(best.models))
  
  # criteria A1/A2
  if("A1" %in% subcriteria | "A2" %in% subcriteria) {
    
    if(length(pop_data) == 1 ) {
      
      Results$reduction_A12 <-
        100 *(1 - (tail(as.numeric(pop_data[[1]]), 1) /        
               head(as.numeric(pop_data[[1]]), 1)))

    } else {
    
      Results$reduction_A12 <-
        100 * sapply(pop_data, function(y) 
          1 - (as.numeric(y[which(names(y) %in% assess.year)]) /
                 as.numeric(y[which(names(y) %in% min(years))])))
    
    }
  }
  
  # criteria A3
  if("A3" %in% subcriteria) {
    
    Results$reduction_A3 <-
      100 * sapply(pop_data, function(y)
        1 - (as.numeric(y[which(names(y) %in% max(yrs))]) /
               as.numeric(y[which(names(y) == assess.year)])))
    
  }
  
  # criteria A4
  if("A4" %in% subcriteria) {
    
    ### Gilles: cannot be computed if no generation.time provided : any assumed value?
    ### Renato: I think now it should be right
    #anos1 <- years[(1 + which(years == min(yrs))):(which(years == assess.year) - 1)]
    anos1 <- yrs[(1 + which(yrs == min(yrs))):(which(yrs == assess.year) - 1)]
    
    if(is.null(generation.time)) {
      
      ids <- which(yrs %in% (anos1 + 10))
      #anos2 <- yrs[which(yrs %in% (anos1 + 10))]
      
    } else {
      
      ids <- which(yrs %in% (anos1 + 3 * generation.time))
      #anos2 <- yrs[which(yrs %in% (anos1 + 3 * generation.time))]
      
    }
    
    if(length(ids) > 0) {
      
      anos2 <- yrs[ids]
      Results$reduction_A4 <-
        100 * sapply(pop_data, function(y)
          max(1 - (as.numeric(y[names(y) %in% as.character(anos2)]) /
                     as.numeric(y[names(y) %in% as.character(anos1)])), na.rm = TRUE))
      
    }
  }  
  
  ## specific function to categorize taxa based on reductions values.
  all_ranks <- cat_criterion_a(
    A1_val = if("A1" %in% subcriteria) Results$reduction_A12 else NULL,
    A2_val = if("A2" %in% subcriteria) Results$reduction_A12 else NULL,
    A3_val = if("A3" %in% subcriteria) Results$reduction_A3 else NULL,
    A4_val = if("A4" %in% subcriteria) Results$reduction_A4 else NULL, 
    A1.threshold = A1.threshold, 
    A234.threshold = A234.threshold,
    all.cats = all.cats
  )
  
  if(all.cats & !is.null(all_ranks$all_cats))
    Results = cbind.data.frame(Results, all_ranks$all_cats,
                               deparse.level = 0,
                               stringsAsFactors = FALSE)
  Results$category_cA <- all_ranks$ranks_A
  Results$category_cA_code <- all_ranks$cats_code
  
  return(Results)
}    
