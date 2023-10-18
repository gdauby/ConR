#' @title Generate dummy distribution data
#'
#' @description 
#'  Generates dummy geographic distribution data for multiple species
#'  
#' 
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @return A data.frame with three columns
#'   
#' @examples 
#' dummy_dist(n = 5) ## five species dataset
#' dummy_dist(xmin = -20, xmax = 30, n = 30)
#'
#' @export dummy_dist
dummy_dist <- function(n = 10, xmin = 0, xmax = 5, ymin = 0, ymax = 5, step = 0.1) {
  nbe_occ <- sample(seq(1, n, 1), n, replace = TRUE)
  ddlat <- lapply(nbe_occ, function(x) sample(seq(xmin, xmax, step), x))
  ddlon <- lapply(nbe_occ, function(x) sample(seq(ymin, ymax, step), x))
  names(ddlon) <- names(ddlat) <- paste0("tax", seq(1,length(ddlat), 1))
  ddlat <- lapply(seq_along(ddlat), function(i) data.frame(ddlat =  ddlat[[i]]))
  ddlon <- lapply(seq_along(ddlon), function(i) data.frame(ddlon =  ddlon[[i]], taxa = names(ddlon)[i]))
  test_data <- 
    data.frame(ddlat = do.call('rbind', ddlat), ddlon = do.call('rbind', ddlon)[,1], taxa = do.call('rbind', ddlon)[,2])
  return(test_data)
}