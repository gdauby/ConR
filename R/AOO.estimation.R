
#' @title Internal function
#'
#' @description AOO estimation
#' 
#' @param coordEAC data.frame
#' @param cell_size integer
#' @param nbe_rep integer
#' @param export_shp logical
#' @inheritParams proj_crs
#' 
#' 
#' @keywords internal
#' @noRd
AOO.estimation <- function(coordEAC,
                           cell_size = 2,
                           nbe_rep = 0,
                           export_shp = FALSE,
                           proj_type = proj_type
) {
  
  res <-
    cell.occupied(
      nbe_rep = nbe_rep,
      size = cell_size,
      coord = coordEAC[,c(2, 1)],
      export_shp = export_shp,
      proj_type = proj_type
    )
  
  
  
  AOO <- res[[2]] * cell_size * cell_size
  
  if (export_shp)
    return(list(AOO = AOO, poly_AOO = res[[1]]))
  
  if (!export_shp)
    return(AOO)
  
}