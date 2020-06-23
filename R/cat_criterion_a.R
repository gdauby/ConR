#' @title Categorize taxa according to IUCN criterion A
#'
#' @description Provide IUCN category based on the four sub-criteria and provided thresholds
#'
#' @param A1_val numeric vector
#' @param A2_val numeric vector
#' @param A3_val numeric vector
#' @param A4_val numeric vector
#' @param A1.threshold numeric vector
#' @param A234.threshold numeric vector
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
#' A1_val <- c(40, 90, 75, 10)
#' A2_val <- c(90, 95, 30, 10)
#' A3_val <- c(10, 30, 45, 15)
#' A4_val <- c(32, 18, 80, 20)

#' cat_criterion_a(
#'   A1_val = A1_val,
#'   A2_val = A2_val,
#'   A3_val = A3_val,
#'   A4_val = A4_val
#' )
#'
cat_criterion_a <- function(A1_val = NULL,
                            A2_val = NULL,
                            A3_val = NULL,
                            A4_val = NULL,
                            A1.threshold = c(50, 70, 90), 
                            A234.threshold = c(30, 50, 80)) {
  
  all.identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))
  
  if (!all.identical(c(
    length(A1_val),
    length(A2_val),
    length(A3_val),
    length(A4_val)
  )))
    stop("Numbers of values provided for each criterion should be identical")
  
  rank_a1 <- findInterval(A1_val, sort(A1.threshold))
  
  rank_a2 <- findInterval(A2_val, sort(A234.threshold))
  
  rank_a3 <- findInterval(A3_val, sort(A234.threshold))
  
  rank_a4 <- findInterval(A4_val, sort(A234.threshold))
  
  all_ranks <-  cbind.data.frame(A1 = rank_a1,
                                 A2 = rank_a2,
                                 A3 = rank_a3,
                                 A4 = rank_a4,
                                 deparse.level = 0, 
                                 stringsAsFactors = FALSE)
  
  ranks_A1234 <- 
    as.character(apply(all_ranks, 1, max, na.rm=TRUE))
  
  ranks_A1234 <- 
    gsub("3", "CR", ranks_A1234)
  ranks_A1234 <- 
    gsub("2", "EN", ranks_A1234)
  ranks_A1234 <- 
    gsub("1", "VU", ranks_A1234)
  ranks_A1234 <- 
    gsub("0", "LC or NT", ranks_A1234)
  
  cat_codes <- apply(
    all_ranks,
    1,
    FUN = function(x) {
      y = names(x[x == max(x, na.rm = T)])
      paste(y[!is.na(y)], collapse = "+")
    }
  )
  
  return(list(ranks_A1234 = ranks_A1234, cat_codes = cat_codes))
}