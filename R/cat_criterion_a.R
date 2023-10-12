#' Categorize taxa according to IUCN criterion A
#'
#' Provide the consensus IUCN category based on the sub-criteria of 
#' IUCN criterion A (A1, A2, A3 and A4) and the thresholds recommended by IUCN.
#'
#' @param A1_val a numeric vector of the estimates of population decline based on IUCN sub-criterion A1
#' @param A2_val a numeric vector of the estimates of population decline based on IUCN sub-criterion A2
#' @param A3_val a numeric vector of the estimates of population decline based on IUCN sub-criterion A3
#' @param A4_val a numeric vector of the estimates of population decline based on IUCN sub-criterion A4
#' @param A1.threshold a numeric vector of threshold values for sub-criterion A1, default is c(50, 70, 90) following IUCN guidelines
#' @param A234.threshold a numeric vector of A2, A3 and A4 threshold values for sub-criteria A2, A3 and A4, default is c(30, 50, 80) following IUCN guidelines
#' @param all.cats a logical value indicating whether to return all categories or only the consensus, default is TRUE
#' 
#' @return A list of two or three elements:
#' \itemize{
#' \item \code{ranks_A}: a character vector with the consensus category of each taxon based on subcriteria A1-A4
#' \item \code{cats_code}: a character vector with the final IUCN category for each taxon
#' \item \code{all_cats}: a data frame with the categorization of each taxon based on all sub criteria, only returned if \code{all.cats=TRUE}
#' }
#' 
#' @details By default, the function provides the consensus category, following the recommendations of IUCN (2019)
#' that states "Only the criteria for the highest category of threat that the taxon qualifies for should be listed".
#' Therefore, the consensus category is the highest category of threat among the sub-criteria evaluated.
#' 
#' The function assumes that the order of the values in A1_val, A2_val, A3_val and A4_val are 
#' from the same taxa (i.e. first element from A1_val until A4_val is always the same species i). 
#' Therefore, the order of the estimates of population decline for each sub-criterion *must* be the same.
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
#' @author Gilles Dauby & Renato A. Ferreira de Lima
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'   
#' @importFrom stringr str_replace_all
#' 
#' @export cat_criterion_a
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
    warning(paste("The following subcriteria were not used in the assessment: ", paste(names(L)[unlist(lapply(L, is.null))], collapse = ", ")),
            call. = FALSE)
  
  if(any(unlist(lapply(L, function(x) any(x<0)))))
    warning(paste("The following subcriteria had population increase for one or more species: ", paste(names(L)[unlist(lapply(L, function(x) any(x<0)))], collapse = ", ")),
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
