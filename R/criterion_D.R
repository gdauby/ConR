#' @title Assess IUCN Criterion D
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion D, which is based on very small or restricted populations 
#'  (Subcriteria D1 and D2).
#'
#' @param pop.size a value (one species) or a vector (multiple species/
#'   subpopulations) containing the population sizes (e.g. number of mature
#'   individuals).
#' @param name_sp a vector containing the name of the taxa to be assessed.
#'   Default to "Species 1", "Species 2", ..., "Species n", where n is the
#'   number of taxa.
#' @param AOO a vector containing the Area of Occupancy of each taxon (in km2).  
#' @param n.Locs a vector containing the number of locations for each taxon.
#' @param prop.mature a value or vector of the proportion of mature individuals in the 
#'   total population (IUCN 2019). Default to 1.
#' @param subcriteria a vector containing the sub-criteria that should be
#'   included in the assessment (i.e. D and/or D2).
#' @param D.threshold numeric vector with the criterion D thresholds of very
#'   small population sizes. Default values are the thresholds recommended by
#'   the IUCN.
#' @param AOO.threshold numeric vector containing the threshold of the Area of
#'   Occupancy (AOO) in km2. Default to NULL.
#' @param Loc.threshold numeric vector containing the threshold of the Number of
#'   Locations. Default to NULL.
#' @param all.cats logical. Should the categories from all criteria be returned
#'   and not just the consensus categories?
#'
#' @return A data frame containing the name of each taxon, the population sizes,
#'   the AOO and Number of Locations, the IUCN categories associated with each
#'   sub-criterion and the consensus category for criterion D.
#' 
#' @details This is a simple and fast function to assess IUCN criterion D. The
#'   assessment based solely on population sizes is done automatically, if a
#'   vector of population sizes is provided. Caution should be taken while
#'   assessing sub-criterion D2. IUCN (2019, p.71) emphasizes that there is no
#'   strict thresholds for D2 and that this sub-criterion should only be
#'   assessed if the "population is prone to the effects of human activities or
#'   stochastic events in an uncertain future, and is thus capable of becoming
#'   Critically Endangered or even Extinct in a very short time period. (...).
#'   So, simply meeting the suggested (or any other) threshold for AOO or number
#'   of locations is not sufficient". Therefore, `criterion_D` assumes no default
#'   thresholds for D2.
#'   
#'   The argument `prop.mature` can be used if the population data provided are not 
#'   already the number of mature individuals (i.e. population size sensu IUCN, 2019). 
#'   By default, the proportion of mature individuals in the total population proportion 
#'   is taken as 1, but the user can provide one proportion for all species or species-
#'   specific proportions.

#'
#'   Currenlty, the function does not supports data separated by subpopulation.
#'   
#' @author Renato A. Ferreira de Lima & Gilles Dauby
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#'
#' @examples
#' 
#' ## All sub-criteria
#' criterion_D(pop.size = c(25, 50, 200, 250, 500, 1000, 1500), 
#'             AOO = c(1,2,5,10,15,20,25), 
#'             n.Locs = c(1,1,2,3,4,5,6),
#'             subcriteria = c("D", "D2"),
#'             AOO.threshold = 20,
#'             Loc.threshold = 5)
#'             
#' ## Only sub-criterion D
#' criterion_D(pop.size = c(25, 50, 200, 250, 500, 1000, 1500), 
#'             subcriteria = c("D"))
#'             
#' ## Olny sub-criterion D2 (Please read the function details before its use)            
#' criterion_D(pop.size = NULL, 
#'             AOO = c(1,2,5,10,15,20,25), 
#'             n.Locs = c(1,1,2,3,4,5,6),
#'             subcriteria = c("D2"),
#'             AOO.threshold = 20,
#'             Loc.threshold = 5)
#'             
#' @importFrom stringr str_replace_all
#'
#' @export criterion_D
criterion_D = function(pop.size = NULL,
                       name_sp = NULL,
                       AOO = NULL,
                       n.Locs = NULL,
                       prop.mature = NULL,
                       subcriteria = c("D", "D2"),
                       D.threshold = c(1000, 250, 50),
                       AOO.threshold = NULL,
                       Loc.threshold = NULL,
                       all.cats = TRUE) {
  
  if(!any(subcriteria %in% c("D", "D2")))
    stop("Please provide at least one sub-criteria for the assessment: D and/or D2")

  if("D" %in% subcriteria) {
    
    if(is.null(pop.size))
      stop("Please provide at least one estimate of population size")
    
  }
  
  if("D2" %in% subcriteria & (!is.null(AOO) | !is.null(n.Locs))) {
    
    if(is.null(AOO.threshold) & is.null(Loc.threshold))
      stop("Please, provide plausible thresholds and make sure that the taxa evaluated are prone to become Critically Endangered or even Extinct in a very short time period to assess sub-criterion D2")
    
  }

  if(is.null(pop.size) & all(subcriteria == "D2"))
    pop.size = rep(NA, max(length(AOO), length(n.Locs)))
  
  # if(is.null(name_sp)) {
    
    x <- as.data.frame(matrix(pop.size, nrow = length(pop.size), 
                             dimnames = list(paste("Species",1:length(pop.size)), "pop.size")),
                      stringsAsFactors = FALSE)
    
  # } else {
  #   
  #   x <- as.data.frame(matrix(pop.size, nrow = length(pop.size), dimnames = list(name_sp, "pop.size")),
  #                     stringsAsFactors = FALSE)
  #   
  # }
  
  if("D2" %in% subcriteria & (!is.null(AOO) | !is.null(n.Locs))) {
    
    if(length(AOO) != length(n.Locs))
      stop("The number of values in AOO and n.Locs is not the same. Please check entry data")
    
    if(!is.null(AOO))  
      x <-  cbind.data.frame(x, AOO = AOO,
                           stringsAsFactors = FALSE, row.names = NULL)
    if(!is.null(n.Locs))  
      x  <- cbind.data.frame(x, n.Locs = n.Locs,
                           stringsAsFactors = FALSE, row.names = NULL)
  }
  
  if(is.null(prop.mature)) {
    
    prop.mature <- rep(1, dim(x)[1])
    
  } else {
    
    if(any(prop.mature>1) | any(prop.mature<0))
      warning("The proportion of mature individuals normally ranges between 0 and 1")
    
    if(dim(x)[1] != length(prop.mature)) {
      
      if(length(unique(prop.mature)) > 1)
        stop("Number of proportions of mature individuals in the population is different from the number of taxa in the assessment. Please provide one value for all taxa or one value for each taxa")
      
      if(length(unique(prop.mature)) == 1) {
        
        prop.mature <- rep(prop.mature, dim(x)[1])
        warning("Only one proportion of mature individuals provided for two or more taxa: assuming the same proportion for all taxa")
        
      }
    }
  } 
  
  x$pop.size <- as.double(x$pop.size) * prop.mature

    ## Very small or restricted population using IUCN criteria
  Results <- data.frame(
    tax = if (!is.null(name_sp)) name_sp else paste("Species",1:length(pop.size)),
    stringsAsFactors = FALSE
  )
  
  rpl.cds <- c("CR", "EN", "VU", "LC or NT")
  names(rpl.cds) <- c("0", "1", "2", "3")
  
  assess_D =  list(D = NULL,
                  D2.AOO = NULL,
                  D2.Loc = NULL)

  ## Criteria D: very small population sizes
  if("D" %in% subcriteria & !is.null(D.threshold)) {

    Results$pop.size <- x$pop.size
    obj = list(ranks = NULL, cats = NULL, codes = rep("D", length(x$pop.size)))
    obj[[1]] <- findInterval(x$pop.size, sort(D.threshold))
    obj[[2]] <- stringr::str_replace_all(obj[[1]], rpl.cds)
    obj[[3]][obj[[2]] == "LC or NT"] <- ""
    assess_D[["D"]] <- obj
    
  } 
  
  ##Criteria D2
  if("D2" %in% subcriteria) { 
    
    if(!is.null(AOO) & !is.null(AOO.threshold)) {
      
      Results$AOO <- x$AOO
      obj = list(ranks = rep(3, length(x$AOO)), cats = NULL, codes = rep("", length(x$AOO)))
      obj[[1]][x$AOO < AOO.threshold] <- 2
      obj[[2]] <- stringr::str_replace_all(obj[[1]], rpl.cds)
      obj[[3]][obj[[2]] == "VU"] <- "D2"
      assess_D[["D.AOO"]] <- obj
      
    }
    
    if(!is.null(n.Locs) & !is.null(Loc.threshold)) {
      
      Results$Numb.Locations <- x$n.Locs
      obj = list(ranks = rep(3, length(x$n.Locs)), cats = NULL, codes = rep("", length(x$n.Locs)))
      obj[[1]][x$AOO < AOO.threshold] <- 2
      obj[[2]] <- stringr::str_replace_all(obj[[1]], rpl.cds)
      obj[[3]][obj[[2]] == "VU"] <- "D2"
      assess_D[["D2.Loc"]] <- obj

    }
  } 
  
  assess_D <- assess_D[!unlist(lapply(assess_D, is.null))]
  ranks <- do.call(cbind.data.frame, lapply(assess_D, function(x) x[[1]]))
  ranks_D <- as.character(apply(ranks, 1, min, na.rm=TRUE))
  ranks_D <- stringr::str_replace_all(ranks_D, rpl.cds)
  
  cats <- do.call(cbind.data.frame, lapply(assess_D, function(x) x[[2]]))
  cats_code <- suppressWarnings(apply(ranks, 1,
                                      FUN = function(x) {
                                        y = names(x[x == min(x[x<3], na.rm = T)])
                                        paste(y, collapse = "+")
                                      }))

  if(all.cats &  inherits(ranks, "data.frame"))
    Results = cbind.data.frame(Results, cats,
                               deparse.level = 0,
                               stringsAsFactors = FALSE)
  Results$category_D <- ranks_D
  Results$category_D_code <- cats_code
  
  return(Results)
}