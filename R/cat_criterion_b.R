#' @title Categorize taxa according to IUCN criterion B
#'
#' @description Provide IUCN category based on the four sub-criteria and provided thresholds
#'
#' @param EOO numeric vector
#' @param AOO numeric vector
#' @param locations numeric vector
#' @param protected numeric vector providing estimated percentage of species distribution in protected areas
#' @param decline string vector providing sub-populations decline status. If different of 'Decreasing', sub-criterion (b) of criterion B will not be met
#' @param protected.threshold numeric, one value indicating the threshold for protected value above which a taxa would not be threatened whatever the others parameters, by default is 100
#' @param EOO.threshold numeric vector
#' @param AOO.threshold numeric vector
#' @param Loc.threshold numeric vector
#'
#' @return list
#' 
#' @details The function categorizes taxa following criterion B and categories of the IUCN. 
#' 
#' @author Dauby, G. & Lima, R.A.F.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @export cat_criterion_b
#'
#' @examples
#' 
#' EOO <- c(34000, 5000)
#' AOO <- c(300, 25)
#' locations <- c(9, 12)
#' cat_criterion_b(EOO = EOO,AOO = AOO, locations = locations)
#' 
#' EOO <- c(34000, 80)
#' AOO <- c(300, 25)
#' locations <- c(1, 1)
#' cat_criterion_b(EOO = EOO, AOO = AOO, locations = locations)
#' 
#' EOO <- c(50, 5000)
#' AOO <- c(5, 25)
#' locations <- c(1, 10)
#' protected <- c(80, 50)
#' decline <- c("Stable", "Decreasing")
#' cat_criterion_b(EOO = EOO, AOO = AOO, locations = locations, protected = protected, decline = decline)
#' 
#' EOO <- c(34000, 5000)
#' AOO <- c(300, 25)
#' locations <- c(9, 12)
#' protected <- c(100, 80)
#' decline <- c("Decreasing", "Decreasing")
#' cat_criterion_b(EOO = EOO, AOO = AOO, locations = locations, protected = protected, decline = decline)
#' 
#' @importFrom utils tail
#' 
cat_criterion_b <- function(EOO = NULL,
                            AOO = NULL,
                            locations = NULL,
                            protected = NULL, ## to include yet
                            decline = NULL,
                            protected.threshold = 100,
                            EOO.threshold = c(20000, 5000, 100), 
                            AOO.threshold = c(2000, 500, 10), 
                            Loc.threshold = c(10, 5, 1)) {
  
  all.identical <-
    function(l)
      all(mapply(identical, head(l, 1), tail(l,-1)))
  
  if (!all.identical(c(
    length(EOO),
    length(AOO),
    length(locations),
    ifelse(is.null(protected), length(EOO), length(protected)),
    ifelse(is.null(decline), length(EOO), length(decline))
  )))
    stop("Numbers of values provided for each parameters should be identical")
  
  if(protected.threshold > 100 | protected.threshold <= 0)
    stop("protected.threshold must be higher than 0 and lower or equal to 100")
  
  # ,
  # length(protected)
  
  rank_eoo <- findInterval(EOO, sort(EOO.threshold))
  
  rank_aoo <- findInterval(AOO, sort(AOO.threshold))
  
  rank_loc <- 
    findInterval(locations, sort(Loc.threshold), left.open = T)

  all_ranks <-  cbind.data.frame(B1a = rank_eoo,
                                 B2a = rank_aoo,
                                 Ba = rank_loc,
                                 deparse.level = 0, 
                                 stringsAsFactors = FALSE)
  
  ranks_B12a <- 
    as.character(apply(all_ranks, 
                       1, FUN = function(x) {
                         
                         min_b <- 
                           y <- min(x[1:2], na.rm = T)
                         
                         y <- 
                           max(c(min_b, x[3]), na.rm = T)
                         
                         # if(x[3] < 3) {
                         #   
                         #  y <- 
                         #    max(c(min_b, x[3]), na.rm = T)
                         #  
                         # } else {
                         #   
                         #   y <-
                         #     min(c(min_b, x[3]), na.rm = T)
                         #   
                         # }
                         
                         return(y)
                       }
                       ))
  
  if(!is.null(protected)) {
    
    if(any(protected >= protected.threshold) & any(ranks_B12a != "3")) {
      
      if(any(ranks_B12a[which(protected >= protected.threshold)] != '3')) {
        
        message("Some taxa categorized as Threatened finally assessed as Not Threatened because percentage of their area in protected areas above the protected.threshold")
        
        ranks_B12a[which(protected >= protected.threshold)] <- 
          "3"
        
      }
    }
  }
  
  if(!is.null(decline)) {
    
    if(any(decline != "Decreasing") & any(ranks_B12a != "3")) {
      
      
      if(any(ranks_B12a[which(decline != "Decreasing")] != '3')) {
        
        message("Some taxa categorized as Threatened based on EOO/AOO/locations finally assessed as Not Threatened because no Decline detected")
        
        ranks_B12a[which(decline != "Decreasing")] <- 
          "3"
        
      }
      
      ranks_B12a[which(decline != "Decreasing")] <- 
        "3"
      
    }
  }
  
  ranks_B12a <- 
    gsub("0", "CR", ranks_B12a)
  ranks_B12a <- 
    gsub("1", "EN", ranks_B12a)
  ranks_B12a <- 
    gsub("2", "VU", ranks_B12a)
  ranks_B12a <- 
    gsub("3", "LC or NT", ranks_B12a)
  # ranks_B12a <- 
  #   gsub("4", "LC or NT (protected)", ranks_B12a)
  
  cat_codes <- 
    apply(
    all_ranks,
    1,
    FUN = function(x) {
      y <- names(x[x == min(x, na.rm = T)])
      paste(y[!is.na(y)], collapse = "+")
    }
  )
  
  cat_codes[ranks_B12a == 'LC or NT'] <- NA
  
  return(list(ranks_B12a = ranks_B12a, cat_codes = cat_codes))
}

