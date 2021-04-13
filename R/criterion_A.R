#' @title Assess IUCN Criterion A
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion A, which is based on population size reductions (Criteria A1, 
#'  A2, A3, and A4)
#'
#' @param x a vector (one species) or a data frame (multiple species/
#'   subpopulations) containing the population size per year, from the oldest 
#'   to the most recent population estimate.
#' @param years a vector containing the years for which the population sizes
#'   are available (i.e. time series). It can be NULL if x contains the years as names.
#' @param assess.year numeric. The year for which the assessment should be performed.
#' @param project.years a vector containing the years for which population sizes
#'   were or should be projected.
#' @param generation.time a value or vector of generation lengths, i.e. the
#'   average age of parents of the current cohort (IUCN 2019).
#' @param models a vector containing the names of the models to be fitted to
#'   species population size to perform projections.
#' @param subcriteria a vector containing the sub-criteria that should be
#'   included in the assessment (i.e. A1, A2, A3 and/or A4).
#' @param data.type a character corresponding to the type of data (IUCN 2019):
#'   "observation", "index" or "AOO_EOO" (only these types are currently
#'   implemented)
#' @param nature.evidence a character corresponding to the nature of evidence (IUCN
#'   2019): "observed", "estimated", "projected", "inferred" or "suspected"
#' @param A1.threshold numeric vector with the A1 thresholds to convert decline
#'   estimates into categories. Default values are the thresholds recommended by the IUCN.
#' @param A234.threshold numeric vector with the A2, A3, and A4 thresholds to
#'   convert decline estimate into categories. Default values are the thresholds
#'   recommended by the IUCN.
#' @param all.cats logical. Should the categories from all criteria be returned
#'   and not just the consensus categories?
#' @param parallel logical. Should calculations be parallelized? Default to 
#'   FALSE.
#' @param NbeCores  integer. Number of cores for parallel computing. Default 
#'   to 2.
#' @param show_progress logical. Should the progress bar be displayed? Default
#'  to TRUE.
#' @param ... other parameters to be passed as arguments for function `pop.decline.fit`    
#'
#' @return A data frame containing, for each taxon, the year of assessment, the
#'   time interval of the assessment (include past and future estimates, if
#'   any), the population size in the interval of assessment, the reduction of
#'   the population size using the chosen sub-criteria (A1, A2, A3, and A4), the
#'   model used to obtain the projections of population size (if used), the IUCN
#'   categories associated with these sub-criteria and the consensus category
#'   for criterion A.
#'
#' @details As described in IUCN (2019), the choice between criteria A1 or A2 depends on 
#'   three conditions: the reduction must be reversible, the causes of the reduction 
#'   must be understood, and the threats must have ceased. "If any of the three conditions 
#'   (reversible and understood and ceased) are not met (...), then A2 should be used 
#'   instead of A1" (IUCN, 2019).
#'
#'   Some important notes. The function can return the predictions of population
#'   estimates for years not in the observed data, based on the fit of a set of different 
#'   statistical models. As stated in IUCN (2019), the model used to make the 
#'   predictions can result in very different estimates. So, it is preferable that
#'   the user choose one or two of the models based on the best available information 
#'   on types of threat (i.e. patterns of exploitation or habitat loss), life history 
#'   and ecology of the taxon being evaluated or any other processes that may contribute 
#'   to population decline. See IUCN (2019) for more details on the assumptions of each model.
#'   The selection of models based solely on their fit to population size should only be used 
#'   for larger time series (Number of observations > 10). 
#'   
#'   Some more technical notes. If `years` is a subset of all the years contained in `x`, 
#'   then `x` is filtered based on `years`. So, make sure you have selected the right years.
#'   If the year of assessment is not given, the most recent year is taken instead. The function 
#'   accepts a single generation length for all species or species-specific generation lengths. 
#'   In the latter case, it is necessary to provide exactly one value for each species analyzed.
#'   Currently, only one assessment year can be assigned for all taxa. Similarly, only
#'   one vector of years with population size available. Thus, it is advised not to mix
#'   taxa with great differences in generation length. 
#'   
#' @author Lima, R.A.F. & Dauby, G.
#' 
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#'
#' @examples
#' ## Simplest example: one species, two observations in time, one subcriterion
#'  pop = c("1970" = 10000, "2000" = 6000)
#'  criterion_A(x = pop,
#'   years = c(1970, 2000), 
#'   assess.year = 2000,
#'   project.years = NULL,
#'   subcriteria = c("A2"),
#'   generation.time = 10)
#'   
#' ## Another example: one species, more observations and subcriteria
#' pop = c("1970" = 10000, "1980" = 8900, "1990" = 7000, "2000" = 6000,
#'         "2030" = 4000)
#' criterion_A(x = pop,
#'   years = c(1970, 1980, 1990, 2000, 2030), 
#'   assess.year = 2000,
#'   project.years = c(2010, 2020, 2030),
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
#' ## The data and criterion A assessment as described in IUCN (2019)
#' #available at: https://www.iucnredlist.org/resources/criterion-a
#' data(example_criterionA)
#' criterion_A(example_criterionA,
#'   years = seq(1970, 2000, by = 2), 
#'   assess.year = 2000,
#'   project.years = seq(2002, 2030, by = 2),
#'   subcriteria = c("A1", "A2", "A3", "A4"),
#'   generation.time = 10)
#'
#' ## Same data and options but assuming different generation length for each taxon
#' criterion_A(example_criterionA,
#'   years = seq(1970, 2000, by = 2), 
#'   assess.year = 2000,
#'   project.years = seq(2002, 2030, by = 2),
#'   subcriteria = c("A1", "A2", "A3", "A4"),
#'   generation.time = c(2,5,10,15,30,50))
#'   
#' ## Same data but with different options (need to use predictions from statistical models)
#' criterion_A(example_criterionA,
#'   years = NULL, 
#'   assess.year = 2010,
#'   project.years = NULL,
#'   subcriteria = c("A1", "A2", "A3", "A4"),
#'   generation.time = 10)
#'   
#' @importFrom utils txtProgressBar setTxtProgressBar head
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export criterion_A
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
                       all.cats = TRUE,
                       parallel = FALSE,
                       NbeCores = 2,
                       show_progress = TRUE,
                       ...) {
  
  if (is.null(x))
    stop("Please provide at least two estimates of population sizes")
  
  if (!any(subcriteria %in% c("A1", "A2", "A3", "A4")))
    stop("Please provide at least one sub-criterion for the assessment: A1, A2, A3 and/or A4")
  
  if (is.matrix(x))
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  if(is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    
    if(is.null(anos)) 
      stop("Please provide at least two years with estimates of population sizes") 
    
    years <- anos
    warning("The years of the population sizes were not given and were taken from the input population data", call. = FALSE)
  }
  
  if(is.vector(x)) {
    
    if(is.null(names(x))) {
      
      x <- as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, years)),
                        stringsAsFactors = FALSE)
      
    } else {
      
      x <- as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, names(x))),
                        stringsAsFactors = FALSE)
      
    }
    
    x <- cbind.data.frame(data.frame(species = "species 1"), x)
    
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

    if(dim(x)[1] != length(generation.time)) {
      
      if(length(unique(generation.time)) > 1)
        stop("Number of generation lengths is different from the number of taxa in the assessment. Please provide one value for all taxa or one value for each taxa")
      
      if(length(unique(generation.time)) == 1) {
        
        generation.time <- rep(generation.time, dim(x)[1])
        warning("Only one generation length provided for two or more taxa: assuming the same generation length for all taxa")
        
      }
    }
    
    prev.year <- assess.year - 3 * generation.time
    proj.year <- assess.year + 3 * generation.time
    
    if(any((3 * generation.time) < 10)) {
      
      prev.year[(3 * generation.time) < 10] <- assess.year - 10
      proj.year[(3 * generation.time) < 10] <- assess.year + 10
      warning("Three times the generation length was smaller than 10 years for one or more species: assuming 10 years")
      
    }   
    
    if(any(proj.year - assess.year > 100)) {
      
      proj.year[proj.year - assess.year > 100] <- assess.year + 100
      warning("Maximum projection of population sizes is more than 100 years into the future: assuming 100 years after the year of assessment")
      
    }
  }

  if(any(subcriteria %in% c("A1", "A2")) & !any(subcriteria %in% c("A3", "A4"))) proj.year = rep(assess.year, length(proj.year))
  
  if(!any(subcriteria %in% c("A1", "A2")) & any(subcriteria %in% c("A3", "A4"))) prev.year = rep(assess.year, length(prev.year))
  
  if(is.null(project.years)) {
    
    all.yrs <- lapply(1:length(prev.year),
                      function(i)
                        prev.year[i]:proj.year[i])
    
    yrs <- lapply(1:length(all.yrs),
                  function(i)
                    all.yrs[[i]][all.yrs[[i]] %in% years])
    
    int <- stats::median(diff(years), na.rm = TRUE)
    
    miss.prev <- sapply(1:length(prev.year),
                        function(i)
                          ! prev.year[i] %in% yrs[[i]])
    
    if(any(miss.prev))
      yrs[miss.prev] <- 
      lapply(1:length(yrs[miss.prev]), 
             function(i) unique(c(seq(prev.year[miss.prev][i], min(yrs[miss.prev][[i]]), by= int), yrs[miss.prev][[i]])))
    
    miss.proj <- sapply(1:length(proj.year), 
                        function(i) !proj.year[i] %in% yrs[[i]])
    if(any(miss.proj))
      yrs[miss.proj] <- 
      lapply(1:length(yrs[miss.proj]), 
             function(i) unique(c(yrs[miss.proj][[i]], seq(max(yrs[miss.proj][[i]]), proj.year[miss.proj][i], by= int))))
    
    #all.yrs <- prev.year:proj.year
    #yrs <- all.yrs[all.yrs %in% years]
    # if (!prev.year %in% yrs)
    #   yrs <- unique(c(seq(prev.year, min(yrs), by= int), yrs))
    # if (!proj.year %in% yrs)
    #   yrs <- unique(c(yrs, seq(max(yrs), proj.year, by= int)))
    
  } else {
    
    yrs <-
      rep(list(unique(c(
        years, project.years
      ))), length(prev.year))
    
    int <- median(diff(years), na.rm = TRUE)
    
    miss.prev <- sapply(1:length(prev.year),
                        function(i)
                          ! prev.year[i] %in% yrs[[i]])
    
    min.prev <- sapply(1:length(prev.year),
                       function(i)
                         prev.year[i] < min(yrs[[i]], na.rm = TRUE))
    
    if(any(miss.prev)) {
      
      ids1 = which(miss.prev + min.prev == 1)
      yrs[ids1] <- lapply(1:length(yrs[ids1]), 
                          function(i) sort(unique(c(prev.year[ids1[i]], yrs[[ids1[i]]]))))
      ids2 = which(miss.prev + min.prev == 2)
      yrs[ids2] <- lapply(1:length(yrs[ids2]), 
                          function(i) unique(c(prev.year[ids2[i]], 
                                               seq(prev.year[ids2[i]], min(yrs[[ids2[i]]]), by= int), 
                                               yrs[[ids2[i]]]))) 
      
    }
    
    miss.proj <- sapply(1:length(proj.year),
                        function(i)
                          ! proj.year[i] %in% yrs[[i]])
    
    max.proj <- sapply(1:length(proj.year),
                       function(i)
                         max(yrs[[i]], na.rm = TRUE) < proj.year[i])
    
    if(any(miss.proj)) {
      
      ids1 = which(miss.proj + max.proj == 1)
      yrs[ids1] <- lapply(1:length(yrs[ids1]), 
                          function(i) sort(unique(c(yrs[[ids1[i]]], proj.year[ids1[i]]))))
      ids2 = which(miss.proj + max.proj == 2)
      yrs[ids2] <- lapply(1:length(yrs[ids2]), 
                          function(i) unique(c(yrs[[ids2[i]]], 
                                               seq(max(yrs[[ids2[i]]]), proj.year[ids2[i]], by= int), 
                                               proj.year[ids2[i]]))) 
      #yrs <- unique(c(years, project.years))
      # if (prev.year < min(yrs))
      #   yrs <- unique(c(seq(prev.year, min(yrs), by= int), yrs))
      # if (proj.year > max(yrs))
      #   yrs <- unique(c(yrs, seq(max(yrs), proj.year, by= int)))
      # 
      # if(length(project.years) == 1 & (max(project.years) - max(years[!years %in% project.years]) >= int))
      #   yrs = c(yrs[1:which.max(years[!years %in% project.years])],
      #           seq(max(years[!years %in% project.years]) + int, project.years - int, by = int),
      #           project.years)
      
    }  
  }
  
  if(class(x[,1]) %in% c("factor", "character")) {
    
    names(x)[-1] <- gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    pop_data <- split(x[ ,grepl("[0-9]", names(x))], f = x[,1])
    
  } else {
    
    names(x) <- gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    nomes <- paste("species", 1:dim(x)[1])
    pop_data <- split(x, f = nomes)
    
  }

  best.models <- NULL
  miss.years <- lapply(1:length(yrs),
                      function(i)
                        ! yrs[[i]] %in% names(pop_data[[i]]))
  
  if(any(sapply(miss.years, any))) {   # Predictions based on population trends
    
    if(length(x) < 3) 
      stop("Too few year intervals to fit a model to population trends")
    
    best.models <- as.list(rep(NA, length(pop_data)))
    
    which.pred <- which(sapply(miss.years, any))
    
    ## Renato: Gilles, il faut peut-etre mettre ici la boucle en dplyr et/ou en paralell
    ## Renato: Done in 14 Aug 2020
    
    cat("Computing the predictions based on population trends...", sep= "\n")
    
    if (parallel) {
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      
      message('Parallel running with ',
              NbeCores, ' cores')
      
      `%d%` <- foreach::`%dopar%`
      
    } else {
      `%d%` <- foreach::`%do%`
    }
    
    if (show_progress) {
      pb <- txtProgressBar(min = 0,
                           max = length(which.pred),
                           style = 3)
      
      progress <- function(n)
        setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
    } else {
      opts <- NULL
    }
    
    x <- NULL
    models.fit <- foreach::foreach(
      x = which.pred,
      .options.snow = opts
    ) %d% {
      #source("C://Users//renato//Documents//raflima//R_packages//ConR//R//pop.decline.fit.R")
      
      if (!parallel & show_progress) setTxtProgressBar(pb, x)
      
      pred.sp <- pop.decline.fit(pop.size = pop_data[[x]], 
                             years = years, 
                             models = models,
                             project.years = yrs[[x]],
                             plot.fit = FALSE
                             #)
                             ,...)
      
      pred.pop.size <- pred.sp$data$Observed
      pred.pop.size[is.na(pred.sp$data$Observed)] <- 
        pred.sp$data$Predicted[is.na(pred.sp$data$Observed)]
      names(pred.pop.size) <- pred.sp$data$Year
      res <- list(pred.pop.size, attributes(pred.sp$best.model)$best.model.name)
      res
    }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    for (i in 1:length(which.pred)) { 
      pop_data[[which.pred[i]]] <- models.fit[[i]][[1]]
      best.models[[which.pred[i]]] <- models.fit[[i]][[2]]
    }  
    
    # for (i in 1:length(which.pred)) {
    #   
    #   pred.sp <- pop.decline.fit(pop.size = pop_data[[which.pred[i]]], 
    #                              years = years, 
    #                              models = models,
    #                              project.years = yrs[[which.pred[i]]],
    #                              plot.fit = FALSE,
    #                              ...)
    #   pred.pop.size <- pred.sp$data$Observed
    #   pred.pop.size[is.na(pred.sp$data$Observed)] <- 
    #     pred.sp$data$Predicted[is.na(pred.sp$data$Observed)]
    #   names(pred.pop.size) <- pred.sp$data$Year
    #   pop_data[[which.pred[i]]] <- pred.pop.size
    #   best.models[[which.pred[i]]] <-  attributes(pred.sp$best.model)$best.model.name
    #   
    # }
  }
  
  assess.period <- lapply(1:length(pop_data), 
                      function(i) paste(unique(sort(c(prev.year[i], assess.year, proj.year[i]))), collapse="-"))
  pop_data1 <- lapply(1:length(pop_data), function(i) pop_data[[i]][names(pop_data[[i]]) %in% yrs[[i]]])
  ps.interval <- sapply(1:length(pop_data1), function(i) paste(
    unique(c(as.character(pop_data1[[i]][which(yrs[[i]] %in% prev.year[i])]),
             as.character(pop_data1[[i]][which(yrs[[i]] %in% assess.year)]),
             as.character(pop_data1[[i]][which(yrs[[i]] %in% proj.year[i])]))), collapse = "-"))
  
  ## Population reduction using IUCN criteria
  ## Gilles: I always read it is better to pre allocate space in data frame instead of appending it, for memory and speed
  Results = data.frame(
    species = names(pop_data),
    assessment.year = assess.year,
    assessment.period = as.character(unlist(assess.period)),
    assessment.pop.sizes = as.character(unlist(ps.interval)),
    predictive.model = NA,
    reduction_A12 = NA,
    reduction_A3 = NA,
    reduction_A4 = NA,
    category_A = NA,
    category_A_code = NA,
    stringsAsFactors = FALSE
  )
  row.names(Results) = NULL
  
  if(!is.null(best.models))
    Results$predictive.model = as.character(unlist(best.models))
  
  # criteria A1/A2
  if("A1" %in% subcriteria | "A2" %in% subcriteria) {
    
    # if(length(pop_data) == 1 ) {
    #   
    #   Results$reduction_A12 <-
    #     100 *(1 - (tail(as.numeric(pop_data[[1]]), 1) /        
    #            head(as.numeric(pop_data[[1]]), 1)))
    # 
    # } else {
    
      Results$reduction_A12 <-
        100 * sapply(1:length(pop_data), function(y) 
          1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)]) /
                 as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% prev.year[y])])))
      
          
    # }
  }
  
  # criteria A3
  if("A3" %in% subcriteria) {
    
    # Results$reduction_A3 <-
    #   100 * sapply(pop_data, function(y)
    #     1 - (as.numeric(y[which(names(y) %in% max(yrs))]) /
    #            as.numeric(y[which(names(y) == assess.year)])))
    Results$reduction_A3 <-
      100 * sapply(1:length(pop_data), function(y)
        1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% proj.year[y])]) /
               as.numeric(pop_data[[y]][which(names(pop_data[[y]]) == assess.year)])))
    
  }
  
  # criteria A4
  if("A4" %in% subcriteria) {
    
    ### Gilles: cannot be computed if no generation.time provided : any assumed value?
    ### Renato: I think now it should be alright now
    #anos1 <- years[(1 + which(years == min(yrs))):(which(years == assess.year) - 1)]
    #anos1 <- yrs[(1 + which(yrs == min(yrs))):(which(yrs == assess.year) - 1)]
    #anos1 <- lapply(yrs, function(x) x[(1 + which(x == min(x))):(which(x == assess.year) - 1)])
    anos1 <- lapply(1:length(yrs),
                    function(y)
                      yrs[[y]][(1 + which(yrs[[y]] == min(prev.year[y]))):(which(yrs[[y]] == assess.year) - 1)])
    
    if(is.null(generation.time)) {
      
      ids <- lapply(1:length(yrs), function(y) which(yrs[[y]] %in% (anos1[[y]] + 10)))
      #ids <- which(yrs %in% (anos1 + 10))
      #anos2 <- yrs[which(yrs %in% (anos1 + 10))]
      
    } else {
      
      ids <- lapply(1:length(yrs), 
                          function(y) {
                            try.yrs = anos1[[y]] + 3 * generation.time[y]
                            sapply(1:length(try.yrs), 
                                   function(j) which.min(abs(yrs[[y]] - try.yrs[j])))
                          })
      #ids <- lapply(1:length(yrs), 
      #              function(y) which(yrs[[y]] %in% (anos1[[y]] + 3 * generation.time[y])))
      #ids <- which(yrs %in% (anos1 + 3 * generation.time))
      #anos2 <- yrs[which(yrs %in% (anos1 + 3 * generation.time))]
      
      if(any((3 * generation.time) < 10)) {
        
        ids1 = which((3 * generation.time) < 10) 
        ids[ids1] <- lapply(1:length(yrs[ids1]), 
                    function(y) which(yrs[ids1][[y]] %in% (anos1[ids1][[y]] + 10)))
      
      }
        
    #   if(any(sapply(ids, length) == 0))  {
    #     
    #     ids1 = which(sapply(ids, length) == 0) 
    #     ids[ids1] <- lapply(1:length(yrs[ids1]), 
    #                         function(y) {
    #                             try.yrs = 0.5 + anos1[ids1][[y]] + 3 * generation.time[ids1][y]
    #                             sapply(1:length(try.yrs), 
    #                                             function(j) which.min(abs(yrs[ids1][[y]] - try.yrs[j])))
    #                             })
    #     
    #   }
    }
    
    if(any(sapply(ids, length)) > 0) {
      
      anos2 <- lapply(1:length(yrs), function(y) yrs[[y]][ids[[y]]])
      dup.yrs <- sapply(anos2, duplicated)
      
      Results$reduction_A4 <-
        100 * sapply(1:length(pop_data), function(y)
          max(1 - (as.numeric(pop_data[[y]][names(pop_data[[y]]) %in% as.character(anos2[[y]])][!dup.yrs[[y]]]) /
                     as.numeric(pop_data[[y]][names(pop_data[[y]]) %in% as.character(anos1[[y]])][!dup.yrs[[y]]])), na.rm = TRUE))

      #anos2 <- yrs[ids]
      # Results$reduction_A4 <-
      #   100 * sapply(pop_data, function(y)
      #     max(1 - (as.numeric(y[names(y) %in% as.character(anos2)]) /
      #                as.numeric(y[names(y) %in% as.character(anos1)])), na.rm = TRUE))
      
    }
  }  
  
  ## specific function to categorize taxa based on reductions values
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
  
  Results$category_A <- all_ranks$ranks_A
  Results$category_A_code <- all_ranks$cats_code
  
  Results <-
    Results[, apply(Results, MARGIN = 2, function(x)
      ! all(is.na(x)))]
  
  return(Results)
}    
