#' @title Downlist Threat Categories
#'
#' @description Perform the downlisting of threat categories which is often
#'   necessary for regional conservation assessments (IUCN 2012).
#'
#' @param cats character. The vector containing the IUCN threat categories.
#' @param down.by numerical. The number of steps to downlist the categories.
#'   Default to 1.
#' @param dd logical. Should the Data Deficient (DD) category be included in the
#'   downlisting? Default to FALSE.
#' @param sign logical. Should the degrees sign be indicated in the downlisted
#'   category. Default to TRUE.
#'
#' @details For regional conservation assessments, the IUCN recommends an
#'   additional step, which relates to the possible "effect of populations of
#'   the same taxon in neighbouring regions on the regional population", and
#'   thus to the possibility of a possible rescue effect (IUCN 2012).
#'
#'   Although the IUCN (2012) considers the exceptional possibility of
#'   uplisting, the function only performs downlisting of the categories of
#'   threat. Downlisting normally is a one-step change in the category (i.e.
#'   from EN to VU), but this can be controlled by the argument ```down.by```
#'   (default to 1).
#'
#'   By default, the Data Deficient category is excluded from the downlisting.
#'   The Least Concern category remains unaltered as well.
#'
#'   Note that "if it is unknown whether or not extra-regional populations
#'   influence the extinction risk of the regional population, the category
#'   (...) should be kept unaltered" (IUCN 2012).
#'
#' @author Renato A. Ferreira de Lima
#'
#' @references IUCN (2012). Guidelines for Application of IUCN Red List Criteria
#'   at Regional and National Levels (Version 4.0). IUCN. Gland, Switzerland and
#'   Cambridge, UK. 41pp.
#'
#'
#' @examples
#' cats <- c("EX","CR","EN","VU","NT","DD","LC","NA","NE")
#' cat_downlist(cats)
#' cat_downlist(cats, down.by = 2)
#' cat_downlist(cats, dd = TRUE)
#' cat_downlist(cats, sign = FALSE)
#' cat_downlist(cats, down.by = 2, dd = TRUE, sign = FALSE)
#' 
#' 
#'
#' @importFrom stringr str_replace_all
#' 
#' @export cat_downlist
cat_downlist <- function(cats = NULL, down.by = 1, dd = FALSE, sign = TRUE){

  rpl.cds <- 
    c("5", "5", "5", "4", "3", "2", "1", "0.5", "0", "0")
  names(rpl.cds) <- 
    c("EW","EX","RE","CR", "EN", "VU", "NT", "DD", "LC", "LC or NT")
  
  cats1 <- 
    stringr::str_replace_all(cats, rpl.cds)
  cats1 <- 
    suppressWarnings(as.double(cats1)) - down.by
  
  cats1[!is.na(cats1) & cats1 < 0] <- 0
  
  rpl.cds1 <- c("CRo", "ENo", "VUo", "NTo", "LCo")
  names(rpl.cds1) <- c("4", "3", "2", "1", "0")
  cats1 <- stringr::str_replace_all(cats1, rpl.cds1)

  if (!dd) 
    cats1[cats %in% "DD"] <- "DD" 
  
  if (any(is.na(cats1)))
    cats1[is.na(cats1)] <- cats[is.na(cats1)]

  if (!sign)
    cats1 <- gsub("o$", "", cats1)

  return(cats1)
}
