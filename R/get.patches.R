#' Internal function
#'
#' Get Occupied Cells and Patches
#'
#' @param points XY data frame
#' @param cell_size numeric
#' @param nbe_rep numeric
#' @param EOO.poly Spatial polygon
#' @param proj_type projected coordinate system (in meters)
#' @param Resol_sub_pop numeric. Defines in kilometers the radius of the circles around each occurrence
#' 
#' @details The ...
#' 
#' "A taxon can be considered to be severely fragmented if most (>50%) of its total area of occupancy 
#' is in habitat patches that are (...) (2) separated from other habitat patches by a large distance."
#'   
#' @examples
#' 
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-48.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-37.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#' 
#' @importFrom sf st_intersects st_as_sf st_crs st_distance
#' 
#' 
get.patches <- function(points, 
                        cell_size = 2,
                        nbe_rep = 0,
                        EOO.poly = NULL, 
                        proj_type = "3857",
                        Resol_sub_pop = 50) {
  
  #proj_type <- proj_crs(proj_type = proj_type)

  ## Constructing the sf object with the XY data	
  points_sf <- sf::st_as_sf(points, coords = c("ddlon", "ddlat"))
  if (is.na(sf::st_crs(points_sf)[[1]]))
    sf::st_crs(points_sf) <- sf::st_crs("epsg:4326") # assuming WSG84 if crs is missing
  
  # poly_sf <- sf::st_as_sf(EOO.poly)
  # if (is.na(sf::st_crs(poly_sf)[[1]]))
  #   sf::st_crs(poly_sf) <- proj_type
  
  ## Transform the sf object to the desired projection 
  points_sf <- sf::st_transform(points_sf, crs = "epsg:3857") ## Gilles, here I could not use the argument itself to provide the crs, due to the projection 'bordel'. Needs checking 
  # poly_sf <- sf::st_transform(poly_sf, crs = proj_type)

  ## Creating and merging the buffers around the points  
  points_buff <- sf::st_buffer(points_sf, dist = Resol_sub_pop*1000)
  points_buff <- sf::st_union(points_buff)
  
  #Creating the raster at the desired resolution (again I had problems with the projection...)
  r <- sf::st_as_sf(raster::raster(raster::extent(points_sf),
                 resolution = cell_size * 1000,
                 crs = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

  # Getting the cells by the species 
  sf::st_intersects(points_sf, r)

  # Getting the cells within the buffers (probably occupied cells)  
  
  # obtenir le numero de habitat patches (gid cells que se touche)

  # calculer numero de patches/numero de pixels occupés (

}