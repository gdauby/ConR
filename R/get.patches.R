#' Internal function
#'
#' Get Occupied Cells and Patches
#'
#' @param points XY data frame
#' @param Cell_size_AOO numeric
#' @param EOO.poly Spatial polygon
#' @param proj_user projected coordinate system (in meters)
#' @param Resol_sub_pop numeric. Defines in kilometers the radius of the circles around each occurrence
#' 
#' @details The ...
#'    
#' @examples
#' 
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-28.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-17.2),
#'                    tax = rep("a", 7),
#'                    valid = c(c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE)),
#'                    stringsAsFactors = FALSE)
#' mydf$classes = as.double(mydf$valid)
#' shp <- SpatialPolygonsDataFrame(EOO.computing(mydf[mydf$valid,], export_shp = TRUE)[[2]], data.frame(tax = "a"))
#' plot(mydf[,2:1])
#' plot(shp, add=TRUE)
#' .over.valid.poly(shp, mydf)
#' .over.valid.poly(shp, mydf, proj_user = 5641)
#' .over.valid.poly(shp, mydf, proj_user = 5641, value = "flag")  
#' 
#' @importFrom sf st_intersects st_as_sf st_crs st_distance
#' @importFrom fields rdist
#' @importFrom dplyr full_join
#' 
#' @export over.valid.poly
#' 
get.patches <- function(points, 
                        Cell_size_AOO = 2,
                        EOO.poly, 
                        proj_user,
                        Resol_sub_pop = 5) {
  
  if (is.null(proj_user)) {
    proj_user <- 3857
    warning("no projected coordinate reference system provided by the user: assuming WGS 84 Pseudo-Mercator (see https://epsg.io)"
    )
  }
  
  poly_sf <- sf::st_as_sf(EOO.poly)
  if (is.na(sf::st_crs(poly_sf)[[1]]))
    sf::st_crs(poly_sf) <- proj_user
  
  points_sf <- sf::st_as_sf(points, coords = c("ddlon", "ddlat"))
  if (is.na(sf::st_crs(points_sf)[[1]]))
    sf::st_crs(points_sf) <- sf::st_crs(poly_sf)
  
  points_sf <- sf::st_transform(points_sf, crs = proj_user)
  poly_sf <- sf::st_transform(poly_sf, crs = proj_user)
  

}