
#' @title Internal function
#'
#' @description AOO estimation
#' 
#' @param coordEAC data.frame
#' @param cell_size integer
#' @param nbe_rep integer
#' @param export_shp logical
#' 
AOO.estimation <- function(coordEAC, 
                            cell_size = 2, 
                            nbe_rep = 0, 
                            # poly_borders = NULL, 
                            export_shp=FALSE) {
  
  crs_proj <- 
    proj_crs()
  
  res <-
    cell.occupied(
      nbe_rep = nbe_rep,
      size = cell_size,
      coord = coordEAC,
      export_shp = export_shp
    )
  
  Corners <- rbind(c(min(coordEAC[, 1]),
                     max(coordEAC[, 1])),
                   c(min(coordEAC[, 2]),
                     max(coordEAC[, 2])))
  
  # if (nbe_rep == 0) {
  #   
  #   Occupied_cells <- vector(mode = "numeric", length = 4)
  #   decal <- c(0, 1, 2, 3)
  #   
  #   for (h in decal) {
  #     ext <-
  #       raster::extent(
  #         floor(Corners[1, 1]) - h * (cell_size * 1000 / 4) - 2 * cell_size * 1000,
  #         floor(Corners[1, 2]) + h * (cell_size * 1000 /
  #                                       4) + 2 * cell_size * 1000,
  #         floor(Corners[2, 1]) - h * (cell_size * 1000 /
  #                                       4) - 2 * cell_size * 1000,
  #         floor(Corners[2, 2]) + h * (cell_size * 1000 /
  #                                       4) + 2 * cell_size * 1000
  #       )
  #     
  #     r <-
  #       raster::raster(ext, resolution = cell_size * 1000, crs = crs_proj)
  #     
  #     r2_AOO <-
  #       raster::rasterize(coordEAC[, 1:2], r)
  #     
  #     OCC <-
  #       length(which(!is.na(raster::values(r2_AOO))))
  #     
  #     Occupied_cells[h + 1] <- OCC
  #     
  #     ### If only one occupied cell, stop the production of raster
  #     if (OCC == 1)
  #       break
  #   }
  #   # h <- decal[which.min(Occupied_cells)]
  #   # Occupied_cells <- min(Occupied_cells)
  # }
  # 
  # if (nbe_rep > 0) {
  #   Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
  #   
  #   for (h in 1:nbe_rep) {
  #     rd.1 <- runif(1) * cell_size * 1000
  #     rd.2 <- runif(1) * cell_size * 1000
  #     
  #     ext = raster::extent(
  #       floor(Corners[1, 1]) - rd.1 - 2 * cell_size * 1000,
  #       floor(Corners[1, 2]) + rd.1 + 2 * cell_size * 1000,
  #       floor(Corners[2, 1]) - rd.2 - 2 * cell_size *
  #         1000,
  #       floor(Corners[2, 2]) + rd.2 + 2 * cell_size * 1000
  #     )
  #     r = raster::raster(ext, resolution = cell_size * 1000, crs = crs_proj)
  #     # r
  #     r2_AOO <- raster::rasterize(coordEAC[, 1:2], r)
  #     OCC <- length(which(!is.na(raster::values(r2_AOO))))
  #     Occupied_cells[h] <- OCC
  #     # rd.1.vec <- c(rd.1.vec, rd.1)
  #     # rd.2.vec <- c(rd.2.vec, rd.2)
  #     if (OCC == 1)
  #       break
  #   }
  #   
  # }
  
  # Occupied_cells <- Occupied_cells[Occupied_cells>0]
  # Occupied_cells <- min(Occupied_cells)
  
  AOO <- res[[2]] * cell_size * cell_size  ### AOO
  if (export_shp)
    return(list(AOO, res[[1]]))
  if (!export_shp)
    return(AOO)
  
}