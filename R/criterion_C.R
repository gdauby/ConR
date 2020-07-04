#' @title Assess IUCN Criterion C
#' 
#' @description Preliminary assessment of species conservation status following
#'   IUCN Criterion C, which is based on small population size (e.g. <10,000
#'   mature individuals) that are declining or may decline in the near future.
#' 
#' 
#' @param x a vector (one species) or a data frame (multiple species/
#'   subpopulations) containing the population sizes (e.g. number of mature
#'   individuals) per year, from the oldest to the most recent estimate.
#' @param years a vector containing the years for which the population sizes is
#'   available (i.e. time series). Can be NULL if x contain years as names.
#' @param assess.year the year for which the assessment should be performed.
#' @param project.years a vector containing the years for which population sizes
#'   were or should be projected.
#' @param generation.time a value or vector of generation lengths, i.e. the
#'   average age of parents of the current cohort (IUCN 2019).
#' @param subpop.size a named list containing the vector of number of mature
#'   individuals per subpopulation. The length of the list must match the length
#'   and order of the taxa being assessed.
#' @param models a vector containing the names of the models to be fitted to
#'   species population data to perform projections.
#' @param subcriteria a vector containing the sub-criteria that should be
#'   included in the assessment (i.e. C1 and/or C2).
#' @param C.threshold numeric vector with the criterion C thresholds to define
#'   small population sizes (e.g. number of mature individuals). Default values
#'   are the thresholds recommended by the IUCN.
#' @param C1.threshold numeric vector with the C1 thresholds of continuing
#'   decline. Default values are the thresholds recommended by the IUCN.
#' @param C2ai.threshold numeric vector with the C2a i thresholds for the
#'   population size of the largest subpopulation. Default are the values
#'   recommended by the IUCN.
#' @param C2aii.threshold numeric vector with the C2a ii thresholds for the
#'   percentage of the population size in the same subpopulation. Default are
#'   the values recommended by the IUCN.
#' @param mag.fluct numerical. Threshold of mean order of magnitude of the
#'   differences between population minima and maxima to classify populations
#'   with extreme fluctuations. Default to 10 as recommended by IUCN (2019).
#' @param high.alter numerical. Threshold of proportion of changes that are
#'   followed by a change in the opposite direction. Default to 80%. Currently NOT implemented.
#' @param all.cats logical. Should the categories from all criteria be returned
#'   and not just the consensus categories?
#' @param ... other parameters to be passed as arguments for function `pop.decline.fit`
#' 
#' @return A data frame containing, for each of taxon, the year of assessment,
#'   the time interval of the assessment (include past and future estimates, if
#'   any), the population sizes in the interval of assessment, the model used to
#'   obtain the projections of population size, the population decline and
#'   subpopulation descritors related to sub-criteria C1 and C2, the IUCN
#'   categories associated with these sub-criteria and the consensus category
#'   for criterion C.
#'   
#' @details The function `criterion_C` is similar to another `ConR` function:
#'   `criterion_A`. The main difference between these functions relies on the
#'   differences between criteria A and C as described by IUCN (2019, p.70):
#'   "criterion C applies only to small populations, the time frame over which
#'   the decline is measured is shorter (...) and the decline rate thresholds
#'   are lower, because the populations are already small".
#'   
#'   Two basic tests are performed for each taxon for the assessment of criterion
#'   C. First, we test if the population is small. By default, we use the
#'   maximum value of the thresholds recommended by IUCN (2019): 10,000 mature
#'   individuals. If the taxon is not below this threshold, the assessment is
#'   not performed. IUCN (2019) does not specify at what time the population
#'   size should be below the threshold. Here, we consider the year of the
#'   assessment.
#'   
#'   Next, we test if population size is actually declining. IUCN (2019, p.43)
#'   defines: "A continuing decline is a recent, current or projected future
#'   decline (...) which is liable to continue unless remedial measures are
#'   taken. (...). Continuing declines at any rate can be used to qualify taxa
#'   under criteria B or C2. Estimated continuing decline (under criterion C1)
#'   has quantitative thresholds, and requires a quantitative estimate, which
#'   can be calculated using the same methods as for population reduction" (i.e.
#'   criterion A). Therefore, function `criterion_C` consider two types of
#'   decline: (i) continuing decline at any rate (sub-criteria C2) and (ii)
#'   estimated continuing decline (sub-criteria C1).
#'   
#'   The first type of decline is defined based on the mean change of population
#'   size between observations (no statistical fit); if the mean change from
#'   the first population size suggests a decline in the population size, then
#'   the population is classified as declining. Although IUCN (2019) considers
#'   declines at any rate, here we consider populations in decline those with an
#'   average decline of 0.1% or more, in order to incorporate small fluctuations
#'   in stable populations. Moreover, although (IUCN 2019, p.43) states that
#'   under criteria C2, "continuing declines can be observed, estimated,
#'   inferred or projected", here we consider only observed, estimated, inferred
#'   before the years of assessment.
#'   
#'   The second type of decline is defined on the statistical models
#'   fitted to the observed and/or projected population data. Once the best model is
#'   selected, the confidence interval of the parameters is computed. If the
#'   parameter estimates indicate a declining trend, then the population is
#'   classified as declining (e.g. the slope parameter of the linear model is
#'   negative, as well as the confidence interval around the slope estimate).
#'   For this type of decline, we consider observed, estimated or projected
#'   (IUCN 2019).
#'   
#'   In the case of taxa with population size per subpopulation, there are two
#'   ways to entering subpopulation information. The first is to provide a named
#'   list with a vector of population sizes of each species at the year of
#'   assessment. The other is to provide population sizes for each subpopulation
#'   in `x`, and repeat the name of the taxon in the first column of `x`. In the
#'   case of subpopulations, the overall reduction in population size is
#'   obtained as recommended by IUCN (2019, p.38) which is average reduction across
#'   all subpopulation, weighted by their initial size.
#'   
#'   As defined by IUCN (2019, p. 44), extreme fluctuations are variations in
#'   population size or area typically greater than one order of magnitude. In
#'   addition, "Fluctuations must be inferred only where there is reasonable
#'   certainty that a population change will be followed by a change in the
#'   reverse direction within a generation or two" IUCN (2019).
#'   
#'   
#'   @author Lima, R.A.F. & Dauby, G.
#'   @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'   
#'   
#'   @export criterion_C
#'   
#'   @examples
#'   
#'   
#'   ## No subpopulations
#'   data(example_criterionC)
#'   
#'   criterion_C(x = example_criterionC,
#'   years = NULL, 
#'   assess.year = 2000,
#'   project.years = NULL,
#'   generation.time = 10,
#'   subpop.size = NULL,
#'   models = c("linear", "quadratic", "exponential", "logistic", "general_logistic")
#'   subcriteria = c("C1", "C2")
#'   )
criterion_C = function(x,
                       years = NULL, 
                       assess.year = NULL, 
                       project.years = NULL,
                       generation.time = NULL,
                       subpop.size = NULL,
                       models = c("linear", "quadratic", "exponential", "logistic", "general_logistic","piecewise"),
                       subcriteria = c("C1", "C2"),
                       #data.type = NULL,
                       #nature.evidence = NULL,   
                       C.threshold = c(10000, 2500, 250),
                       C1.threshold = c(10, 20, 25),
                       C2ai.threshold = c(1000, 250, 50),
                       C2aii.threshold = c(90, 95, 100),
                       mag.fluct = 10,
                       high.alter = 80, 
                       all.cats = TRUE,
                       ...) {
  
  if(is.null(x))
    stop("Please provide at least two estimates of population sizes")
  
  if(!any(subcriteria %in% c("C1", "C2")))
    stop("Please provide at least one sub-criteria for the assessment: C1 and/or C2")

  if(is.null(years)) {
    
    anos <- as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    
    if(is.null(anos)) 
      stop("Please provide at least two years with estimates of population sizes") 
    
    years <- anos
    warning("The years of the population sizes were not given and were taken from the input population data", call. = FALSE)
  }
  
  if(is.vector(x)) {
    
    if(is.null(names(x))) {
      
      x <-
        as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, years)),
                      stringsAsFactors = FALSE)
      
    } else {
      
      x <-
        as.data.frame(matrix(x, ncol = length(x), dimnames = list(NULL, names(x))),
                      stringsAsFactors = FALSE)
      
    }
    
    x = cbind.data.frame(data.frame(species = "species 1"), x)
    
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
  
  if("C2" %in% subcriteria) {
    
    if(is.null(subpop.size)) {
      
      if(any(duplicated(x[,1]))) {
      
        numb.subpop <- stats::setNames(as.vector(table(x$species)), nm = unique(x$species))
        
        subpop.size <- split(x[ , which(names(x) == assess.year)], f = x[,1])
        
        x <-
          cbind.data.frame(data.frame (species = unique(x$species)),
                           rowsum(x[, which(names(x) %in% years)], x$species, reorder = FALSE, row.names = FALSE))
        
        row.names(x) <- NULL 
    
      } else {
      
        stop("Please provide the number of individuals for each subpopulation to assess sub-criterion C2 or select only subcriterion C1")
       
      }
      
    } else {
      
      if(!is.list(subpop.size))
        stop("Subpop.size must be a list with length matching the number of taxa of the input data frame 'x'")
      
      if(length(subpop.size) != dim(x)[1])
        stop("The length of subpop.size does not match the number of taxa in the input data frame 'x'")
      
      numb.subpop <- stats::setNames(sapply(subpop.size, length), nm = unique(x$species))
      if (all(numb.subpop == 1))
        stop(
          "Please provide the number of individuals for each subpopulation to assess sub-criterion C2 or select only subcriterion C1"
        )
      if (any(x[, which(names(x) == assess.year)] != sapply(subpop.size, sum)))
        stop(
          "The overall population size provided in 'x' does not match the sum of the subpopulation sizes for one or more taxa. Please, double-check the input data"
        )
      if (is.null(names(subpop.size)) &
          class(x[, 1]) %in% c("factor", "character"))
        names(subpop.size) = unique(x$species)
      warning(
        "Taxon(a) name(s) of 'subpop.size' were not given and were taken from the input population data"
      )
      if (is.null(names(subpop.size)) &
          !class(x[, 1]) %in% c("factor", "character"))
        names(subpop.size) = paste("species", 1:dim(x)[1])
      warning("Taxon(a) name(s) of 'subpop.size' were not given and were created by 'ConR'")
      
    }  
  }
  
  if(is.null(generation.time)) {
    
    prev.year1 <- assess.year - 3
    prev.year2 <- assess.year - 5
    prev.year3 <- assess.year - 10
    proj.year1 <- assess.year + 3
    proj.year2 <- assess.year + 5
    proj.year3 <- assess.year + 10
    warning(
      "Generation length not given: assuming the IUCN defaults (3, 5 and 10 years). Please, check if this is accurate for your species"
    )
    
  } else {
    
    if(dim(x)[1] != length(generation.time)) {
      
      if (length(unique(generation.time)) > 1)
        stop(
          "Number of generation lengths is different from the number of taxa in the assessment. Please provide one value for all taxa or one value for each taxa"
        )
      
      if (length(unique(generation.time)) == 1) {
        generation.time = rep(generation.time, dim(x)[1])
        warning(
          "Only one generation length provided for two or more taxa: assuming the same generation length for all taxa"
        )
        
      }
    }
    
    prev.year1 <- assess.year - 1 * generation.time
    prev.year2 <- assess.year - 2 * generation.time
    prev.year3 <- assess.year - 3 * generation.time
    proj.year1 <- assess.year + 1 * generation.time
    proj.year2 <- assess.year + 2 * generation.time
    proj.year3 <- assess.year + 3 * generation.time
    
    if (any((3 * generation.time) < 10)) {
      prev.year3[(3 * generation.time) < 10] <- assess.year - 10
      proj.year3[(3 * generation.time) < 10] <- assess.year + 10
      warning(
        "Three times the generation length was smaller than 10 years for one or more species: assuming 10 years"
      )
      
    }
    
    if (any((2 * generation.time) < 5)) {
      prev.year2[(2 * generation.time) < 5] <- assess.year - 5
      proj.year2[(2 * generation.time) < 5] <- assess.year + 5
      warning(
        "Two times the generation length was smaller than 5 years for one or more species: assuming 5 years"
      )
      
    }
    
    if (any(generation.time < 3)) {
      prev.year1[generation.time < 3] <- assess.year - 3
      proj.year1[generation.time < 3] <- assess.year + 3
      warning("Generation length was smaller than 3 years for one or more species: assuming 3 years")
      
    }
    
    if (any((proj.year1 - assess.year) > 100) |
        any((proj.year2 - assess.year) > 100) |
        any((proj.year3 - assess.year) > 100))
      warning(
        "Maximum year to project population sizes is more than 100 years into the future: assuming 100 years after the year of assessment"
      )
    
    if (any((proj.year1 - assess.year) > 100))
      proj.year1[(proj.year1 - assess.year) > 100] <-
      assess.year + 100
    
    if (any((proj.year2 - assess.year) > 100))
      proj.year2[(proj.year2 - assess.year) > 100] <-
      assess.year + 100
    
    if (any((proj.year3 - assess.year) > 100))
      proj.year3[(proj.year3 - assess.year) > 100] <-
      assess.year + 100
    
  }
  
  if ("C1" %in% subcriteria) {
    
    yrs <- lapply(1:length(prev.year1), 
                      function(i) c(prev.year3[i], prev.year2[i], prev.year1[i], 
                                          assess.year, 
                                          proj.year1[i], proj.year2[i], proj.year3[i]))

  } else {
    
    yrs <- rep(list(unique(c(years, project.years))), length(prev.year1))

  }
  
  if(class(x[,1]) %in% c("factor", "character")) {
    
    names(x)[-1] <- gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    pop_data <- split(x[ ,grepl("[0-9]", names(x))], f = x[,1])
    
  } else {
    
    names(x) <- gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    nomes <- paste("species", 1:dim(x)[1])
    pop_data <- split(x, f = nomes)
    
  }
  
  ## Continuing decline at any rate
  any.decline <- sapply(pop_data, 
                        function(x) {
                          x1 = as.numeric(x[1:which(names(x) == assess.year)])
                          x1 = x1[!is.na(x1)]
                          mean(diff(x1), na.rm=TRUE) / 
                            head(x1, 1)})
  
  ## Estimated continuing decline
  if("C1" %in% subcriteria) {
  
    if(length(x) < 3) 
      stop("Too few year intervals to fit a model to population trends")
    
    models.fit <- as.list(rep(NA, length(pop_data)))
  
    ## Renato: Gilles, il faut peut-etre mettre ici la boucle en dplyr et/ou en paralell
    ## Je propose d'utiliser foreach et parallel comme dans criterion_b.
    for(i in 1:length(pop_data)) {
    
    models.fit[[i]] <- pop.decline.fit(pop.size = pop_data[[i]], 
                                 years = years, 
                                 models = models,
                                 project.years = yrs[[i]],
                                 plot.fit = FALSE,
                                 ...)
    }
    
    cont.decline <- sapply(models.fit, pop.decline.test, assess.year = assess.year)
  
  }
  
  miss.years = lapply(1:length(yrs), 
                      function(i) !yrs[[i]] %in% names(pop_data[[i]]))
  
  best.models = NULL
  
  if(any(sapply(miss.years, any))) {   # Predictions based on the best model fit to population trends
    
    best.models = as.list(rep(NA, length(pop_data)))
    which.pred = which(sapply(miss.years, any))
    
    for(i in 1:length(which.pred)) {
      
      pred.sp <- models.fit[[which.pred[i]]]
      pred.pop.size <- pred.sp$data$Observed
      pred.pop.size[is.na(pred.sp$data$Observed)] <- 
        pred.sp$data$Predicted[is.na(pred.sp$data$Observed)]
      names(pred.pop.size) <- pred.sp$data$Year
      pop_data[[which.pred[i]]] <- pred.pop.size
      best.models[[which.pred[i]]] <-  attributes(pred.sp$best.model)$best.model.name
      
    }
  }
  
  assess.period <- lapply(1:length(pop_data),
                          function(i)
                            paste(unique(sort(
                              c(prev.year3[i],
                                assess.year,
                                proj.year3[i])
                            )), collapse = "-"))
  pop_data1 <-
    lapply(1:length(pop_data), function(i)
      pop_data[[i]][names(pop_data[[i]]) %in% yrs[[i]]])
  ps.interval <- sapply(1:length(pop_data1), function(i)
    paste(unique(
      c(
        as.character(pop_data1[[i]][which(names(pop_data1[[i]]) %in% prev.year3[i])]),
        as.character(pop_data1[[i]][which(names(pop_data1[[i]]) %in% assess.year)]),
        as.character(pop_data1[[i]][which(names(pop_data1[[i]]) %in% proj.year3[i])])
      )
    ), collapse = "-"))
  
  ## Small population size and continuing decline using IUCN criteria
  Results <- data.frame(
    species = names(pop_data),
    assessment.year = assess.year,
    assessment.period = as.character(unlist(assess.period)),
    assessment.pop.sizes = as.character(unlist(ps.interval)),
    predictive.model = NA,
    stringsAsFactors = FALSE
  )
  row.names(Results) <- NULL
  
  if(!is.null(best.models))
    Results$predictive.model[which.pred] <- as.character(unlist(best.models))[which.pred]

  ## Population size at the assessmente
  Results$assess.pop.size <- sapply(1:length(pop_data), 
                                    function(i) as.numeric(pop_data[[i]][which(names(pop_data[[i]]) %in% assess.year)]))
  
  ## Are population declining at any rate?
  Results$any.decline <- 
    sapply(1:length(any.decline), function(y) if(as.numeric(any.decline[[y]]) < -0.001) "decreasing" else "not.decreasing")
  
  ## Estimated continuing decline?
  Results$cont.decline <- 
    sapply(1:length(cont.decline), function(y) if(cont.decline[[y]] %in% c("signif.decline")) { 
          "decreasing" 
        } else {
          if (cont.decline[[y]] %in% c("non.signif.decline","decrease")) "probably.decreasing" else "not.decreasing"
        })  

  
  ## Criteria C1: under criterion C1, the decline must be observed or estimated (thus removing projections of future decline)
  if("C1" %in% subcriteria) {
    
    Results$reduction_3gen <- 100 * sapply(1:length(pop_data), function(y) 
        1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)]) /
               as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% prev.year3[y])])))
    Results$reduction_2gen <- 100 * sapply(1:length(pop_data), function(y) 
        1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)]) /
               as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% prev.year2[y])])))
    Results$reduction_1gen <- 100 * sapply(1:length(pop_data), function(y) 
        1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)]) /
               as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% prev.year1[y])])))
    # Results$reduction_gen1 <- 100 * sapply(1:length(pop_data), function(y) 
    #     1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% proj.year1[y])]) /
    #            as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)])))
    # Results$reduction_gen2 <- 100 * sapply(1:length(pop_data), function(y) 
    #     1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% proj.year2[y])]) /
    #            as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)])))
    # Results$reduction_gen3 <- 100 * sapply(1:length(pop_data), function(y) 
    #     1 - (as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% proj.year3[y])]) /
    #            as.numeric(pop_data[[y]][which(names(pop_data[[y]]) %in% assess.year)])))
      
  }  
  
  ## Criteria C2: Under criteria B1b, B2b, and C2, continuing declines can be observed, estimated, inferred or projected
  if("C2" %in% subcriteria) {
   
    Results$max.subpop.size <- 
      sapply(subpop.size, max, na.rm = TRUE)
    Results$prop.subpop.size <- 
      sapply(subpop.size, function(x) 100*max(x, na.rm = TRUE)/sum(x, na.rm = TRUE) )

    fluctuations <- 
      t(sapply(1:length(pop_data), function(i) pop.fluctuation(x = pop_data[[i]], years = years, plot.test = FALSE)))
    Results$mean.fluctuation <- as.numeric(fluctuations[,"Magnitude.fluctuation"])
    Results$alternance <- as.numeric(fluctuations[,"Alternance.prop"])
    # Results$extreme.fluctuation <- 
    #   sapply(1:length(pop_data), function(y) if(as.numeric(fluctuations[,"Magnitude.fluctuation"][y]) >= mag.fluct) "yes" else "no")
    # Results$high.alternance <-
    #   sapply(1:length(pop_data), function(y) if(as.numeric(fluctuations[,"Alternance.prop"][y]) >= high.alter) "yes" else "no")
        
  }  
  
  ## specific function to categorize taxa based on reductions values
  if("C1" %in% subcriteria) C1 <- cbind.data.frame(Results[,c("assess.pop.size", "cont.decline")],
                                                   Results[,grepl("reduction", names(Results))],
                                                   stringsAsFactors = FALSE)
  if("C2" %in% subcriteria) C2 <- cbind.data.frame(Results[,c("assess.pop.size", "any.decline","max.subpop.size","prop.subpop.size","mean.fluctuation","alternance")],
                                                   stringsAsFactors = FALSE)
  
  all_ranks <- cat_criterion_c(
    C1_df = if("C1" %in% subcriteria) C1 else NULL,
    C2_df = if("C2" %in% subcriteria) C2 else NULL,
    C.threshold = C.threshold,
    C1.threshold = C1.threshold, 
    C2ai.threshold = C2ai.threshold,
    C2aii.threshold = C2aii.threshold,
    mag.fluct = mag.fluct,
    high.alter = high.alter,
    all.cats = all.cats
  )
  
  if(all.cats & !is.null(all_ranks$all_cats))
    Results = cbind.data.frame(Results, all_ranks$all_cats,
                               deparse.level = 0,
                               stringsAsFactors = FALSE)
  Results$category_C <- all_ranks$ranks_C
  Results$category_C_code <- all_ranks$cats_code
  
  return(Results)
  
}
