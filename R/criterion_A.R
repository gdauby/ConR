#' @title Assess IUCN Criterion A
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion A, which is based on population size reductions (Criteria A1, 
#'  A2, A3 and A4)
#'
#' @param x a vector or a data frame containing the (estimated) number of mature individuals of the species
#' @param years a vector containing the years for which the number of mature individuals was estimated 
#' @param assess.year  
#' @param project.years a vector containing the years for which the number of mature individuals should be predicted
#' @param generation.time 
#' @param models a vector containing the names of the models to be fitted to species population data 
#' @param A1.threshold 
#' @param A234.threshold
#'
#' @return TO BE COMPLETED
#' 
#' @details TO BE COMPLETED 
#' 
#' @author Lima, R.A.F. & Dauby, G.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and Criteria. Version 14. Standards and Petitions Committee. Downloadable from: http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @export criterion_A
#'
#' @examples
#' data(example_criterionA)
#' criterion_A(example_criterionA, years = NULL, assess.year = 2000, project.years = seq(2002, 2030, by = 2), generation.time = 10)
#' 
criterion_A = function(x, 
                       years = NULL, 
                       assess.year = NULL, 
                       project.years = NULL,
                       generation.time = NULL,
                       models = c("linear", "accelerating", "exponential", "logistic", "general_logistic"),
                       A1.threshold = c(50, 70, 90),
                       A234.threshold = c(30, 50, 80)) {

  if(is.null(years)) { 
    anos = as.numeric(gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))]))
    if(is.null(anos)) { 
      stop("please provide at least two years with estimates of population sizes") 
    }
  }
  
  if(!is.null(years) & any(!anos %in% years)) {
    years = sort(unique(c(years, anos[!anos %in% years]))) 
  } else {
    years = anos  
  }
  
  if(dim(x)[2] < 2 | length(years) < 2) {
    stop("At least two years are needed to estimate population reduction")
  }  
  
  if(dim(x)[2] < 3 & !is.null(project.years)) {
    stop("Too few year intervals to fit a model to population trends")
  }  
 
  if(is.null(assess.year)) { 
    assess.year = max(years, na.rm = TRUE)
    warning("Year of assessment not given: assuming to be the most recent year of the time series")
  }
  
  if(max(project.years) - max(years)>100){
    stop("Projected population sizes are larer then 100 years from last year of observation")
  }  
  
  if(is.null(generation.time)) {
    closest.year = years[which.min(abs(years - (assess.year - 10)))]
    yrs = years[which(years == closest.year):which(years == assess.year)]
    warning("generation time not given: assuming 10 years previous to the year of assessment. Please, check if this is accurate for your species")
  } else {
    all.yrs = (assess.year - 3*generation.time):assess.year
    yrs = all.yrs[all.yrs %in% years]
    if(!min(all.yrs) %in% years) yrs = c(min(all.yrs), yrs)
    if(!max(all.yrs) %in% years) yrs = c(yrs, max(all.yrs))
  }   
  
  if(class(x[,1]) %in% c("factor", "character")) {
    names(x)[-1] = gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    pop_data <- split(x[ ,grepl("[0-9]", names(x))], f = x$species)
  } else {
    names(x) = gsub("[^0-9]", "", names(x)[grepl("[0-9]", names(x))])
    nomes = paste("species", 1:dim(x)[1])
    pop_data <- split(x[ ,grepl("[0-9]", names(x))], f = nomes)
  }

  if(!is.null(project.years) & any(!project.years %in% years)) {   # Predictions based on population trends
    for(i in 1:length(pop_data)) {
      pred.sp = pop.decline.fit(pop.size = pop_data[[i]], 
                      years = years, 
                      models = models,
                      project.years = project.years,
                      plot = FALSE)
      pred.pop.size = pred.sp$data$Observed
      pred.pop.size[is.na(pred.sp$data$Observed)] = pred.sp$data$Predicted[is.na(pred.sp$data$Observed)]
      names(pred.pop.size) = pred.sp$data$Year
      pop_data[[i]] = pred.pop.size
    }
  }

  ## Population reduction using IUCN criteria
    # criteria A1/A2
      reduction_A12 = sapply(pop_data, function(y)
          1 - (y[which(names(y) == assess.year)]/
                 y[which(names(y) == min(years))]))
  
    # criteria A3 
      reduction_A3 = sapply(pop_data, function(y)
        1 - (y[which(names(y) == max(project.years))])/
                y[which(names(x) == assess.year)])
      
    # criteria A4
      anos1 = years[(1 + which(years == min(years))):(which(years == assess.year) - 1)]
      anos2 = project.years[which(project.years %in% (anos1 + 3*generation.time))]
      reduction_A4 = sapply(pop_data, function(y)
        max(1 - (y[names(y) %in% as.character(anos2)]/
              y[names(y) %in% as.character(anos1)]), na.rm=TRUE))
      
    Results = cbind.data.frame(species = names(pop_data),
                           reduction_A12 = 100*as.numeric(unlist(reduction_A12)),
                           reduction_A3 = 100*as.numeric(unlist(reduction_A3)),
                           reduction_A4 = 100*as.numeric(unlist(reduction_A4)),
                           deparse.level = 0, stringsAsFactors = FALSE)
    
    Rank_A1 = rep(4, dim(Results)[1])
    Rank_A1[Results$reduction_A12 >= A1.threshold[1]] <-3
    Rank_A1[Results$reduction_A12 >= A1.threshold[2]] <-2
    Rank_A1[Results$reduction_A12 >= A1.threshold[3]] <-1
    
    Rank_A2 = rep(4, dim(Results)[1])
    Rank_A2[Results$reduction_A12 >= A234.threshold[1]] <-3
    Rank_A2[Results$reduction_A12 >= A234.threshold[2]] <-2
    Rank_A2[Results$reduction_A12 >= A234.threshold[3]] <-1

    Rank_A3 = rep(4, dim(Results)[1])
    Rank_A3[Results$reduction_A3 >= A234.threshold[1]] <-3
    Rank_A3[Results$reduction_A3 >= A234.threshold[2]] <-2
    Rank_A3[Results$reduction_A3 >= A234.threshold[3]] <-1
    
    Rank_A4 = rep(4, dim(Results)[1])
    Rank_A4[Results$reduction_A4 >= A234.threshold[1]] <-3
    Rank_A4[Results$reduction_A4 >= A234.threshold[2]] <-2
    Rank_A4[Results$reduction_A4 >= A234.threshold[3]] <-1

    ## Implemente a user defined choice of A1 or A2!!!  
    Rank_A1234 <- apply(rbind(Rank_A1, Rank_A2, Rank_A3, Rank_A4), 2, min, na.rm=TRUE)
    
    Cat = rep("LC or NT", dim(Results)[1])
    Cat[Rank_A1234 == 1] <- "CR"
    Cat[Rank_A1234 == 2] <- "EN"
    Cat[Rank_A1234 == 3] <- "VU"

    #### CHECK HERE: NOT SURE IF ALL CONDITIONS ARE EXAUSTIVE ####    
    Cat_Code <- rep(NA, dim(Results)[1])
    if(any(Rank_A1 > Rank_A2))    
      Cat_Code[Rank_A1 > Rank_A2] <- paste(Cat[Rank_A1 > Rank_A2], "A2")
    if(any(Rank_A2 > Rank_A3)) 
      Cat_Code[Rank_A2 > Rank_A3] <- paste(Cat[Rank_A2 > Rank_A3], "A3")
    if(any(Rank_A3 > Rank_A4)) 
      Cat_Code[Rank_A3 > Rank_A4] <- paste(Cat[Rank_A3 > Rank_A4], "A4")
    
    if (any(Rank_A2 == Rank_A3)) 
      Cat_Code[Rank_A2 == Rank_A3 & Rank_A2 != 4] <- paste(Cat[Rank_A2 == Rank_A3 & Rank_A2 != 4], "A2+A3") 
    if (any(Rank_A2 == Rank_A4)) 
      Cat_Code[Rank_A2 == Rank_A4 & Rank_A2 != 4] <- paste(Cat[Rank_A2 == Rank_A4 & Rank_A2 != 4], "A2+A4") 
    if (any((Rank_A2 == Rank_A3) & (Rank_A2 == Rank_A4))) 
      Cat_Code[(Rank_A2 == Rank_A3) & (Rank_A2 == Rank_A4) & Rank_A2 != 4] <- paste(Cat[(Rank_A2 == Rank_A3) & (Rank_A2 == Rank_A4) & Rank_A2 != 4], "A2+A3+A4") 
    
    Results$Category_A1 = NA
    if (any(Rank_A1 == 1)) 
      Results$Category_A1[Rank_A1 == 1] <- "CR"
    if (any(Rank_A1 == 2)) 
      Results$Category_A1[Rank_A1 == 2] <- "EN"
    if (any(Rank_A1 == 3)) 
      Results$Category_A1[Rank_A1 == 3] <- "VU"
    if (any(Rank_A1 > 3)) 
      Results$Category_A1[Rank_A1 > 3] <- "LC or NT"
  
    Results$Category_A2 = NA
    if (any(Rank_A2 == 1)) 
      Results$Category_A2[Rank_A2 == 1] <- "CR"
    if (any(Rank_A2 == 2)) 
      Results$Category_A2[Rank_A2 == 2] <- "EN"
    if (any(Rank_A2 == 3)) 
      Results$Category_A2[Rank_A2 == 3] <- "VU"
    if (any(Rank_A2 > 3)) 
      Results$Category_A2[Rank_A2 > 3] <- "LC or NT"
    
    Results$Category_A3 = NA
    if (any(Rank_A3 == 1)) 
      Results$Category_A3[Rank_A3 == 1] <- "CR"
    if (any(Rank_A3 == 2)) 
      Results$Category_A3[Rank_A3 == 2] <- "EN"
    if (any(Rank_A3 == 3)) 
      Results$Category_A3[Rank_A3 == 3] <- "VU"
    if (any(Rank_A3 > 3)) 
      Results$Category_A3[Rank_A3 > 3] <- "LC or NT"
    
    Results$Category_A4 = NA
    if (any(Rank_A4 == 1)) 
      Results$Category_A4[Rank_A4 == 1] <- "CR"
    if (any(Rank_A4 == 2)) 
      Results$Category_A4[Rank_A4 == 2] <- "EN"
    if (any(Rank_A4 == 3)) 
      Results$Category_A4[Rank_A4 == 3] <- "VU"
    if (any(Rank_A4 > 3)) 
      Results$Category_A4[Rank_A4 > 3] <- "LC or NT"
    
    Results$Category_CriteriaA <- Cat
    Results$Category_code <- Cat_Code
    
    return(Results)  
}    
