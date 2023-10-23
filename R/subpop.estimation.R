#' @title Internal function
#'
#' @description Subpopulations estimation
#' 
#' @param XY data.frame with first two columns are projected coordinates (numeric)
#' @param Resol_sub_pop numeric. Defines the radius of the circles
#'   around each occurrence, in kilometres.
#' @param export_shp logical. Whether the resulting shapefiles should be
#'   exported. FALSE by default.
#' @param proj_type character string or numeric or object of CRS class, by
#'   default is "cea"
#'
#' @import sf
#' @keywords internal
#' @noRd
subpop.estimation <- function(XY,
                         Resol_sub_pop, 
                         export_shp = FALSE,
                         proj_type = "cea") {
  
  points_sf <- st_as_sf(XY, coords = c(2, 1))
  buff_sf <- st_buffer(points_sf, dist = Resol_sub_pop * 1000)
  buff_sf <- st_union(buff_sf) 
  buff_sf <- st_cast(buff_sf, "POLYGON")
  SubPopPoly <-
    st_as_sf(data.frame(buff_sf))

  st_crs(SubPopPoly) <- proj_type
  
  NbeSubPop <- nrow(SubPopPoly)

  if (export_shp) {
    OUTPUT <- list(NbeSubPop, SubPopPoly)
    names(OUTPUT) <- c("number_subpop", "poly_subpop")
  
  } else {
    OUTPUT <- NbeSubPop
    names(OUTPUT) <- c("number_subpop")
  }
  
  return(OUTPUT)
}
