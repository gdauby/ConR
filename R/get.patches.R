#' Internal function
#'
#' Get Occupied Cells and Patches
#'
#' @param points XY data frame
#' @param cell_size numeric
#' @param nbe_rep numeric
#' @param proj_type projected coordinate system (in meters)
#' @param Resol_sub_pop numeric. Defines in kilometres the radius of the circles around each occurrence
#' @param dist_isolated numeric. Distance in kilometres for identifying subpopulations considered to be isolated
#' 
#' @details 
#' This function evaluates the proportion of the total area of occupancy 
#' in habitat patches separated from others by a large distance.
#' Based on IUCN guidelines page(48):
#' "A taxon can be considered to be severely fragmented if most (>50%) of its total area of occupancy 
#' is in habitat patches that are (...) (2) separated from other habitat patches by a large distance."
#' 
#' This function interpret subpopulations as obtained using ```subpop.estimation``` function as habitat patches. 
#' First, subpopulations that are isolated i.e. distant to all others subpopulations by at ```least dist_isolated``` kilometres, are identified.
#' Then, the proportion of the AOO concerned by those "isolated" subpopulations is calculated.
#'   
#' @examples
#' 
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-48.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-37.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#' 
#' get.patches(XY = mydf, dist_isolated = 200)
#' 
#' @import sf
#' 
#' 
get.patches <- function(XY, 
                        cell_size = 2,
                        nbe_rep = 0,
                        proj_type = "cea",
                        Resol_sub_pop = 10,
                        dist_isolated) {
  
  proj_type <- 
    proj_crs(proj_type = proj_type)
  
  XY_proj <- 
    st_as_sf(XY, coords = c("ddlon", "ddlat"))
  st_crs(XY_proj) <- 4326
  XY_proj <- st_transform(XY_proj, crs = proj_type)
  XY_proj_coord <- st_coordinates(XY_proj)
  XY_proj_coord <- as.data.frame(XY_proj_coord)
  
  res_aoo <- 
    AOO.estimation(
    coordEAC = XY_proj_coord[, c(2, 1)],
    cell_size = cell_size,
    nbe_rep = nbe_rep,
    export_shp = TRUE,
    proj_type = proj_type
  )
  
  res_aoo_poly <- 
    res_aoo$poly_AOO
  
  res_subpop <-
    subpop.estimation(
      XY = XY_proj_coord,
      Resol_sub_pop = Resol_sub_pop,
      proj_type = proj_type,
      export_shp = TRUE
    )
  
  res_subpop_poly <- 
    res_subpop$poly_subpop
  
  mapview::mapview(res_subpop_poly) + mapview::mapview(res_aoo_poly, col.regions = "red")
  
  dist_btw_poly <- st_distance(res_subpop_poly)
  dist_btw_poly <- matrix(dist_btw_poly, nrow = nrow(dist_btw_poly), ncol =nrow(dist_btw_poly))
  
  # dist_btw_poly[which(dist_btw_poly > 100000)]
  
  
  dist_btw_poly[upper.tri(dist_btw_poly, diag = T)] <- NA
  
  mat_above_thres <- dist_btw_poly > dist_isolated*1000
  
  isolated_subpop1 <- apply(mat_above_thres, 1, FUN = function(x) all(x, na.rm = T))
  isolated_subpop2 <- apply(mat_above_thres, 2, FUN = function(x) all(x, na.rm = T))
  
  isolated_subpop_poly <- res_subpop_poly[isolated_subpop1 & isolated_subpop2,]
  nrow(isolated_subpop_poly)
  connected_subpop_poly <- res_subpop_poly[!(isolated_subpop1 & isolated_subpop2),]
  nrow(connected_subpop_poly)
  # plot(connected_subpop_poly)
  # plot(isolated_subpop_poly, add = T, col = "red")
  
  
  intersect_aoo_isolated <- 
    st_intersects(res_aoo_poly, isolated_subpop_poly)
  
  fraction_aoo <- 
    length(unlist(intersect_aoo_isolated))/nrow(res_aoo_poly)*100
  
  return(list(fraction_aoo = fraction_aoo, isolated_subpop_poly = isolated_subpop_poly))
  
  # mapview::mapview(connected_subpop_poly) + mapview::mapview(isolated_subpop_poly, col.regions = "red")
  # 
  # ## Constructing the sf object with the XY data	
  # points_sf <- sf::st_as_sf(XY, coords = c("ddlon", "ddlat"))
  # if (is.na(sf::st_crs(points_sf)[[1]]))
  #   sf::st_crs(points_sf) <- sf::st_crs("epsg:4326") # assuming WSG84 if crs is missing
  # 
  # # poly_sf <- sf::st_as_sf(EOO.poly)
  # # if (is.na(sf::st_crs(poly_sf)[[1]]))
  # #   sf::st_crs(poly_sf) <- proj_type
  # 
  # ## Transform the sf object to the desired projection 
  # points_sf <- sf::st_transform(points_sf, crs = "epsg:3857") ## Gilles, here I could not use the argument itself to provide the crs, due to the projection 'bordel'. Needs checking 
  # # poly_sf <- sf::st_transform(poly_sf, crs = proj_type)
  # 
  # ## Creating and merging the buffers around the points  
  # points_buff <- sf::st_buffer(points_sf, dist = Resol_sub_pop*1000)
  # points_buff <- sf::st_union(points_buff)
  # 
  # #Creating the raster at the desired resolution (again I had problems with the projection...)
  # r <- sf::st_as_sf(raster::raster(raster::extent(points_sf),
  #                resolution = cell_size * 1000,
  #                crs = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  # 
  # # Getting the cells by the species 
  # sf::st_intersects(points_sf, r)

  # Getting the cells within the buffers (probably occupied cells)  
  
  # obtenir le numero de habitat patches (gid cells que se touche)

  # calculer numero de patches/numero de pixels occupés (

}