#' @title Categorize taxa according to IUCN criterion A
#' 
#' @description Provide the consensus IUCN category based on the sub-criteria of 
#' IUCN criterion A (A1, A2, A3 and A4) and the thresholds recommended by IUCN.
#'
#' @param A1_val numeric vector of the estimates of population decline based on IUCN sub-criteria A1.
#' @param A2_val numeric vector of the estimates of population decline based on IUCN sub-criteria A2.
#' @param A3_val numeric vector of the estimates of population decline based on IUCN sub-criteria A3.
#' @param A4_val numeric vector of the estimates of population decline based on IUCN sub-criteria A4.
#' @param A1.threshold numeric vector with the A1 thresholds to convert decline estimates into categories. Default is the thresholds recommended by IUCN.
#' @param A234.threshold numeric vector with the A2, A3 and A4 thresholds to convert decline estimate into categories. Default is the thresholds recommended by IUCN.
#' @param all.cats logical. Should the categories from all criteria be returned and not just the consensus categories? Default to TRUE.
#' 
#' @return A list containing a vector of the consensus category from all sub-criteria evaluated for each taxon (`ranks_A`) and
#' the sub-criteria used to obtain the consensus category (`cat_codes`). If `all.cats  == TRUE` the function also returns a
#' data frame containing the categories classified by each sub-criteria individually (`all.cats`).
#' 
#' @details By default, the function provides the consensus category, following the recommendations of IUCN (2019)
#' that states "Only the criteria for the highest category of threat that the taxon qualifies for should be listed".
#' Therefore, the consensus category is the highest category of threat among the sub-criteria evaluated.
#' 
#' The function assumes that the order of the values in A1_val, A2_val, A3_val and A4_val are 
#' from the same taxa (i.e. first element from A1_val until A4_val is always the same species i). 
#' Therefore, the order of the estimates of population decline for each sub-criteria *must* be the same.
#' 
#' @author Dauby, G. & Lima, R.A.F.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and Criteria. Version 14. Standards and Petitions Committee. Downloadable from: http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @importFrom stringr str_replace_all
#' 
#' @export cat_criterion_a
#'
#' @examples
#' A1_val <- c(40, 90, 75, 10)
#' A2_val <- c(90, 95, 30, 10)
#' A3_val <- c(10, 30, 45, 15)
#' A4_val <- c(32, 18, 80, 20)
#' 
#' ## All four criteria
#' cat_criterion_a(
#'   A1_val = A1_val,
#'   A2_val = A2_val,
#'   A3_val = A3_val,
#'   A4_val = A4_val
#' )
#' 
#' ## All four criteria, no categories for all criteria, just the consesus
#' cat_criterion_a(
#'   A1_val = A1_val,
#'   A2_val = A2_val,
#'   A3_val = A3_val,
#'   A4_val = A4_val,
#'   all.cats = FALSE
#' )
#' 
#' ## One or more criteria not evaluated
#' cat_criterion_a(
#'   A1_val = NULL,
#'   A2_val = A2_val,
#'   A3_val = A3_val,
#'   A4_val = A4_val
#' )
#' 
#' cat_criterion_a(
#'   A1_val = A1_val,
#'   A2_val = NULL,
#'   A3_val = NULL,
#'   A4_val = NULL
#' )
#'
#' ## One or more criteria not evaluated, one species
#' cat_criterion_a(
#'   A1_val = NULL,
#'   A2_val = 0.4,
#'   A3_val = NULL,
#'   A4_val = NULL
#' )
#'
cat_criterion_a <- function(A1_val = NULL,
                            A2_val = NULL,
                            A3_val = NULL,
                            A4_val = NULL,
                            A1.threshold = c(50, 70, 90), 
                            A234.threshold = c(30, 50, 80),
                            all.cats = TRUE) {
  
  L <- list(A1 = A1_val,
            A2 = A2_val,
            A3 = A3_val,
            A4 = A4_val)
  rank <- cat <- L[!unlist(lapply(L, is.null))]
  if(length(unique(lapply(rank, length)))>1)
    stop("Numbers of values provided for each criterion should be identical")
  
  if(any(unlist(lapply(L, is.null))))
    warning(paste("The following criteria were not provided and are not used in the assessment: ", paste(names(L)[unlist(lapply(L, is.null))], collapse = ", ")),
            call. = FALSE)
  
  # all.identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))
  # 
  # if (!all.identical(c(
  #   length(A1_val),
  #   length(A2_val),
  #   length(A3_val),
  #   length(A4_val)
  # )))
  #   stop("Numbers of values provided for each criterion should be identical")
  
  rpl.cds <- c("CR", "EN", "VU", "LC or NT")
  names(rpl.cds) <- c("3", "2", "1", "0")
  
  if(!is.null(A1_val)) {
    
    rank[["A1"]] <- findInterval(rank[["A1"]], sort(A1.threshold))
    cat[["A1"]] <- stringr::str_replace_all(rank[["A1"]], rpl.cds)
    #rank_a1 <- findInterval(A1_val, sort(A1.threshold)) 
    
  }
  
  if(!is.null(A2_val)|!is.null(A3_val)|!is.null(A4_val)) {  
    
    ids <- !names(rank) %in% "A1"
    rank[ids] <- lapply(rank[ids], findInterval, vec = sort(A234.threshold))
    cat[ids] <- lapply(rank[ids], stringr::str_replace_all, rpl.cds)
    #rank_a2 <- findInterval(A2_val, sort(A234.threshold))
    #rank_a3 <- findInterval(A3_val, sort(A234.threshold))
    #rank_a4 <- findInterval(A4_val, sort(A234.threshold))
    
  }    
  
  all_ranks <- do.call(cbind.data.frame, rank)
  all_cats  <- do.call(cbind.data.frame, cat)
  
  # all_ranks <-  cbind.data.frame(A1 = rank_a1,
  #                                A2 = rank_a2,
  #                                A3 = rank_a3,
  #                                A4 = rank_a4,
  #                                deparse.level = 0, 
  #                                stringsAsFactors = FALSE)
  
  ranks_A1234 <- 
    as.character(apply(all_ranks, 1, max, na.rm=TRUE))
  ranks_A1234 <- 
    stringr::str_replace_all(ranks_A1234,  rpl.cds)
  
  # ranks_A1234 <- 
  #   gsub("3", "CR", ranks_A1234)
  # ranks_A1234 <- 
  #   gsub("2", "EN", ranks_A1234)
  # ranks_A1234 <- 
  #   gsub("1", "VU", ranks_A1234)
  # ranks_A1234 <- 
  #   gsub("0", "LC or NT", ranks_A1234)
  
  cats_code <- apply(
    all_ranks,
    1,
    FUN = function(x) {
      y = names(x[x == max(x, na.rm = T)])
      paste(y[!is.na(y)], collapse = "+")
    }
  )
  
  if(all.cats & dim(all_cats)[2] > 1) {
    
    return(list(ranks_A = ranks_A1234, cats_code = cats_code, all_cats = all_cats))
    
  } else {
    
    return(list(ranks_A = ranks_A1234, cats_code = cats_code))
    #return(list(ranks_A1234 = ranks_A1234, cat_codes = cat_codes))
    
  }
  
}
