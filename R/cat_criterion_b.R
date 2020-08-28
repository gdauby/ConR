#' @title Categorize taxa according to IUCN criterion B
#'
#' @description Provide IUCN threat categories based on B sub-criteria, conditions and thresholds.
#'
#' @param EOO numeric vector with species extent of occurrence - EOO (i.e. sub-criterion B1)
#' @param AOO numeric vector with species area of occupancy - AOO (i.e. sub-criterion B2)
#' @param locations numeric vector with the number of locations where the species occur (i.e. condition 'a')
#' @param protected numeric vector providing estimated percentage of species distribution in protected areas
#' @param decline string vector providing the status of the species continuing decline in EOO, AOO, habitat,
#'  locations or subpopulations or population size (i.e. condition 'b'). If different of 'Decreasing', 
#'  the species is classified as condition 'b' of criterion B will not be met.
#' @param ext.fluct numeric. vector with the mean order of magnitude of the
#'   differences between population minima and maxima. Currently not implemented.
#' @param EOO.threshold numeric vector with the EOO thresholds to convert estimates into threat categories. 
#'  Default is the thresholds recommended by IUCN.
#' @param AOO.threshold numeric vector with the AOO thresholds to convert estimates into threat categories. 
#'  Default is the thresholds recommended by IUCN.
#' @param Loc.threshold numeric vector with the thresholds of number of locations (condition 'a'). 
#'  Default is the thresholds recommended by IUCN.
#' @param protected.threshold numeric, one value indicating the threshold for protected value above which 
#'  a taxa would not be threatened whatever the others parameters, by default is 100
#' @param fluct.threshold numeric. Threshold of the order of magnitude of the
#'   differences between population minima and maxima to classify extreme fluctuations. 
#'   Default to 10 as recommended by IUCN.
#' @param all.cats logical. Should the categories from all criteria be returned and not just the consensus categories? Default to TRUE.
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
                            protected = NULL,
                            decline = NULL,
                            ext.fluct = NULL,
                            EOO.threshold = c(20000, 5000, 100), 
                            AOO.threshold = c(2000, 500, 10), 
                            Loc.threshold = c(10, 5, 1),
                            protected.threshold = 100,
                            fluct.threshold = 10,
                            all.cats = TRUE
                            ) {
  
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
  
  rank_eoo <- findInterval(EOO, sort(EOO.threshold))
  
  rank_aoo <- findInterval(AOO, sort(AOO.threshold))
  
  rank_loc <- 
    findInterval(locations, sort(Loc.threshold), left.open = T)

  

    all_ranks <-  cbind.data.frame(B1a = rank_eoo,
                                   B2a = rank_aoo,
                                   Ba = rank_loc,
                                   deparse.level = 0, 
                                   stringsAsFactors = FALSE)    

  
  if(!is.null(decline)) {
    
    names(all_ranks) <- 
      paste0(names(all_ranks), "b")
    
  } else {
    
    message("No information on decline range provided, continuing decline is assumed to be true")
    
  }
  
  
  all_missing <- 
    apply(all_ranks, 1, function(x) ifelse(all(is.na(x)), TRUE, 
                                           ifelse(all(is.na(x[1:2])), TRUE, FALSE)))
  
  if(sum(all_missing) > 0) {
    warning(paste(sum(all_missing), "taxa are not categorized because EOO & AOO or EOO & AOO & locations are missing"))
  }
  
  ranks_B12a <- 
    as.character(apply(all_ranks[!all_missing, ], 
                       1, FUN = function(x) {
                         
                         min_b <-
                           min(x[1:2], na.rm = T)
                         
                         y <- 
                           max(c(min_b, x[3]), na.rm = T)
                         
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
        
        message("Some taxa categorized as threatened based on EOO/AOO/locations were finally assessed as not threatened because no Decline detected")
        
        ranks_B12a[which(decline != "Decreasing")] <- 
          "3"
        
      }
      
      ranks_B12a[which(decline != "Decreasing")] <- 
        "3"
      
    }
  }
  
  replace_code <- 
    data.frame(code = c("0", "1", "2", "3"),
             cat = c("CR", "EN", "VU", "LC or NT"))
  
  for (i in 1:nrow(replace_code))
    ranks_B12a <-
    gsub(replace_code[i, 1], replace_code[i, 2], ranks_B12a)
  
  # ranks_B12a <- 
  #   gsub("0", "CR", ranks_B12a)
  # ranks_B12a <- 
  #   gsub("1", "EN", ranks_B12a)
  # ranks_B12a <- 
  #   gsub("2", "VU", ranks_B12a)
  # ranks_B12a <- 
  #   gsub("3", "LC or NT", ranks_B12a)
  
  if(all.cats) {
    ranks_B1a <- 
      as.character(apply(all_ranks, 
                         1, FUN = function(x) {
                           
                           min_b <- min(x[1])
                           
                           y <- 
                             max(c(min_b, x[3]))
                           
                           return(y)
                         }
      ))
    
    ranks_B2a <- 
      as.character(apply(all_ranks, 
                         1, FUN = function(x) {
                           
                           min_b <- min(x[2])
                           
                           y <- 
                             max(c(min_b, x[3]))
                           
                           return(y)
                         }
      ))
    
    for (i in 1:nrow(replace_code))
      ranks_B1a <-
      gsub(replace_code[i, 1], replace_code[i, 2], ranks_B1a)
    
    for (i in 1:nrow(replace_code))
      ranks_B2a <-
      gsub(replace_code[i, 1], replace_code[i, 2], ranks_B2a)
  }
  
  cat_codes <- 
    apply(
      all_ranks[!all_missing,][,1:2],
      1,
      FUN = function(x) {
        y <- names(x[x == min(x, na.rm = T)])
        paste(y[!is.na(y)], collapse = "+")
      }
    )
  
  ranks_B12a_final <- cat_codes_final <- 
    vector(mode = "character", length = nrow(all_ranks))
  
  ranks_B12a_final[!all_missing] <- 
    ranks_B12a
  ranks_B12a_final[all_missing] <- 
    NA
  
  cat_codes_final[!all_missing] <- 
    cat_codes
  cat_codes_final[all_missing] <- 
    NA
  
  cat_codes_final[ranks_B12a_final == 'LC or NT'] <- NA
  
  if(!all.cats) return(list(cat_cb = ranks_B12a_final, cat_codes_cb = cat_codes_final))
  
  if(all.cats) return(list(cat_cb = ranks_B12a_final, cat_codes_cb = cat_codes_final, ranks_B1a = ranks_B2a, ranks_B2a = ranks_B2a ))
  
}

