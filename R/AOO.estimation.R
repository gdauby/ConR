
#' @title Internal function
#'
#' @description AOO estimation
#' 
#' @param coordEAC data.frame
#' @param cell_size integer
#' @param nbe_rep integer
#' @param export_shp logical
#' @param proj_type character string
#' 
AOO.estimation <- function(coordEAC,
                           cell_size = 2,
                           nbe_rep = 0,
                           # poly_borders = NULL,
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
  
  # mbm <- microbenchmark("1" = { res <-
  #   cell.occupied(
  #     nbe_rep = nbe_rep,
  #     size = cell_size,
  #     coord = coordEAC[,c(2, 1)],
  #     export_shp = export_shp,
  #     proj_type = proj_type
  #   )
  #                       },
  #                       "2" = {
  #                         res <-
  #                           .cell.occupied.stars(
  #                             nbe_rep = nbe_rep,
  #                             size = cell_size,
  #                             coord = coordEAC[,c(2, 1)],
  #                             export_shp = export_shp,
  #                             proj_type = proj_type
  #                           )
  #                       },
  #                       check = NULL)
  
  
  
  # Corners <- rbind(c(min(coordEAC[, 1]),
  #                    max(coordEAC[, 1])),
  #                  c(min(coordEAC[, 2]),
  #                    max(coordEAC[, 2])))
  
  
  AOO <- res[[2]] * cell_size * cell_size  ### AOO
  
  if (export_shp)
    return(list(AOO = AOO, poly_AOO = res[[1]]))
  
  if (!export_shp)
    return(AOO)
  
}