#' @title Categorize taxa according to IUCN criterion B
#'
#' @description Provide IUCN category based on the four sub-criteria and provided thresholds
#'
#' @param EOO numeric vector
#' @param AOO numeric vector
#' @param locations numeric vector
#' @param protected numeric vector
#' @param EOO.threshold numeric vector
#' @param AOO.threshold numeric vector
#' @param Loc.threshold numeric vector
#'
#' @return list
#' 
#' @details The function ... 
#' 
#' @author Dauby, G. & Lima, R.A.F.
#'
#' @references 
#'
#' @export
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
#' 
cat_criterion_b <- function(EOO = NULL,
                            AOO = NULL,
                            locations = NULL,
                            protected = NULL, ## to include yet
                            EOO.threshold = c(20000, 5000, 100), 
                            AOO.threshold = c(2000, 500, 10), 
                            Loc.threshold = c(10, 5, 1)) {
  
  all.identical <-
    function(l)
      all(mapply(identical, head(l, 1), tail(l,-1)))
  
  if (!all.identical(c(
    length(EOO),
    length(AOO),
    length(locations)
  )))
    stop("Numbers of values provided for each parameters should be identical")
  
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
                         
                         if(x[3] > 10) {
                          y <- min(c(min_b, x[3]), na.rm = T)
                         }else{
                           y <- max(c(min_b, x[3]), na.rm = T)
                         }
                         return(y)
                       }
                       ))
  
  ranks_B12a <- 
    gsub("0", "CR", ranks_B12a)
  ranks_B12a <- 
    gsub("1", "EN", ranks_B12a)
  ranks_B12a <- 
    gsub("2", "VU", ranks_B12a)
  ranks_B12a <- 
    gsub("3", "LC or NT", ranks_B12a)
  
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

