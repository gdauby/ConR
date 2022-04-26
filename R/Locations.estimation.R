#' @title Internal function
#'
#' @description Number of location estimation
#' 
#' @param coordEAC data.frame
#' @param cell_size integer
#' @param nbe_rep integer
#' @param export_shp logical
#' @param proj_type character string
#' @param Rel_cell_size numeric, if \code{method_locations="sliding scale"}, \code{Cell_size_locations} is ignored and the resolution is given by the maximum distance separating two occurrences multiplied by \code{Rel_cell_size}. By default, it is 0.05
#' 
Locations.estimation <- function(coordEAC,
                                 cell_size = 10,
                                 nbe_rep = 0,
                                 # poly_borders = NULL,
                                 export_shp = FALSE,
                                 proj_type = proj_type,
                                 method = "fixed_grid",
                                 Rel_cell_size = 0.05
                                 
) {
  
  
  if (any(method == "sliding scale")) {
    
    if (nrow(coordEAC) > 1) {
      
      pairwise_dist <- stats::dist(coordEAC[, 1:2],  upper = F)
      
      cell_size <- max(pairwise_dist) * Rel_cell_size / 1000
      
    } else{
      
      cell_size <- 10
      
    }
  }
  
  res <-
    cell.occupied(
      nbe_rep = nbe_rep,
      size = cell_size,
      coord = coordEAC[,c(2, 1)],
      export_shp = export_shp,
      proj_type = proj_type
    )
  
  locations <- res[[2]]
  
  if (export_shp)
    return(list(locations = locations, poly_locations = res[[1]]))
  
  if (!export_shp)
    return(locations)
  
}