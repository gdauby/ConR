#' @title Internal function
#'
#' @description Compute prop and nbr taxa per cell
#'
#' @param Cell_count data.frame
#' @param threshold integer
#' 
.prop_threat <- function(Cell_count, threshold) {
  NbeRec <- nrow(Cell_count)
  if(NbeRec >= threshold) {
    NbeEsp <- length(unique(Cell_count$tax))
    NbeThreatened <- length(unique(Cell_count[which(Cell_count$Category_CriteriaB %in% c("CR","EN","VU")),"tax"]))
    PropThreatened <- round(NbeThreatened/NbeEsp*100,1)
  }else{
    NbeEsp <- NbeThreatened <- PropThreatened <- NA
  }
  c(NbeRec, NbeEsp, NbeThreatened, PropThreatened)
}