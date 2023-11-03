#' @title Internal function
#'
#' @description 
#'  `r lifecycle::badge("experimental")`
#'  Get Occupied Cells and Patches
#'
#' @param points XY data frame
#' @param cell_size numeric
#' @param nbe_rep numeric
#' @param AOO numeric. AOO of the species in squared kilometres.
#' @param Resol_sub_pop numeric. Defines in kilometres the radius of the circles
#'   around each occurrence
#' @param subpop_poly Simple feature collection, output of
#'   ```subpop.estimation``` function
#' @param dist_isolated numeric. Distance in kilometres for identifying
#'   subpopulations considered to be isolated
#' @param proj_type projected coordinate system (in meters)
#' @param export_shp logical
#' 
#' @details This function evaluates the proportion of the total area of
#' occupancy in habitat patches separated from others by a large distance. Based
#' on IUCN guidelines page(48): "A taxon can be considered to be severely
#' fragmented if most (>50%) of its total area of occupancy is in habitat
#' patches that are (...) (2) separated from other habitat patches by a large
#' distance."
#'
#' This function interpret subpopulations as obtained using
#' ```subpop.estimation``` function as habitat patches. First, subpopulations
#' that are isolated i.e. distant to all others subpopulations by at least
#' ```dist_isolated``` kilometres, are identified. Then, the percentage of the
#' AOO concerned by those "isolated" subpopulations is calculated.
#'   
#' @examples
#' 
#' \donttest{
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-48.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-37.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#' 
#' get.patches(XY = mydf, dist_isolated = 200)
#' }
#' 
#' 
#' @keywords internal
#' 
#' @importFrom sf st_as_sf st_transform st_coordinates st_distance st_geometry
#' 
#' @noRd
get.patches <- function(XY, 
                        cell_size = NULL,
                        nbe_rep = 0,
                        AOO = NULL,
                        Resol_sub_pop = NULL,
                        subpop_poly = NULL,
                        dist_isolated = NULL,
                        proj_type = "cea",
                        export_shp = FALSE) {
  
  if (!is.null(AOO) & is.null(cell_size))
    stop("Please provide the size (in km) of the grid cells used for AOO")
  
  proj_type <- 
    proj_crs(proj_type = proj_type)

  ## Creating and projecting the points  
  XY_proj <-
    sf::st_as_sf(XY, coords = c("ddlon", "ddlat"))
  sf::st_crs(XY_proj) <- 4326   
  XY_proj <- sf::st_transform(XY_proj, crs = proj_type)
  XY_proj_coord <- sf::st_coordinates(XY_proj)
  XY_proj_coord <- as.data.frame(XY_proj_coord)
  
  
  ## if no AOO is provided, getting the AOO raster for the points
  if (is.null(AOO)) {
    res_aoo <-
      #For Casearia sylvestris (AOO = 3104): 24.5 secs; For Annona neosalicifoli (AOO = 604): 1.11 secs
      AOO.estimation(
        coordEAC = XY_proj_coord[, c(2, 1)],
        cell_size = cell_size,
        nbe_rep = nbe_rep,
        export_shp = TRUE,
        proj_type = proj_type
      )
    
    res_aoo_poly <-
      res_aoo$poly_AOO
    cells <- nrow(res_aoo_poly)
    
  } else {
    cells <- round(AOO / (cell_size * cell_size), 0)
  }  
  ## if no subpopulation sf provided, getting it
  if (is.null(subpop_poly)) {
      
    res_subpop <-          #For Casearia sylvestris: 0.30 secs; For Annona neosalicifolia: 0.07 secs
      subpop.estimation(
        XY = XY_proj_coord,
        Resol_sub_pop = Resol_sub_pop,
        proj_type = proj_type,
        export_shp = TRUE
      )

    res_subpop_poly <-
      res_subpop$poly_subpop
    
  } else {
    
    res_subpop_poly <- 
      subpop_poly
    
  }
  
  # mapview::mapview(res_subpop_poly) + mapview::mapview(res_aoo_poly, col.regions = "red")
  
  ##Obtaining the distances between subpopulations above the dist_isolated threshold
  dist_btw_poly <- sf::st_distance(res_subpop_poly) #For Casearia sylvestris (11 subpops): 0.85 secs; For Annona neosalicifolia (12 subpops): 0.45 secs 

  dist_btw_poly <- matrix(dist_btw_poly,
            nrow = nrow(dist_btw_poly), ncol = nrow(dist_btw_poly))
  dist_btw_poly[upper.tri(dist_btw_poly, diag = TRUE)] <- NA
  
  mat_above_thres <- dist_btw_poly > dist_isolated * 1000
  
  isolated_subpop1 <- 
    apply(mat_above_thres, 1, FUN = function(x) all(x, na.rm = T))
  isolated_subpop2 <- 
    apply(mat_above_thres, 2, FUN = function(x) all(x, na.rm = T))
  
  isolated_subpop_poly <-
    res_subpop_poly[isolated_subpop1 & isolated_subpop2, ]
  
  if (nrow(isolated_subpop_poly) > 0) {
    df <- data.frame(tax = as.character(unique(XY$tax)),
                     frag = "isolated", stringsAsFactors = FALSE)
    isolated_subpop_poly <-
      sf::st_as_sf(data.frame(sf::st_geometry(isolated_subpop_poly), df))
  }

  connected_subpop_poly <- res_subpop_poly[!(isolated_subpop1 & isolated_subpop2),]
  
  if (nrow(connected_subpop_poly) > 0) {
    df <- data.frame(tax = as.character(unique(XY$tax)),
                     frag = "connected")
    connected_subpop_poly <- 
      sf::st_as_sf(data.frame(sf::st_geometry(connected_subpop_poly), df))
  }

  # nrow(connected_subpop_poly)
  # plot(sf::st_geometry(connected_subpop_poly))
  # plot(isolated_subpop_poly[,1], add = T, col = "red")
  
  #### GILLES: WHY THIS STEP IS NECESSARY? CAN'T YOU JUST DO: ####
  fraction_aoo <- 100 * nrow(isolated_subpop_poly)/cells 
  
  # intersect_aoo_isolated <- 
  #   sf::st_intersects(res_aoo_poly, isolated_subpop_poly)
  # 
  # fraction_aoo <- 
  #   length(unlist(intersect_aoo_isolated))/nrow(res_aoo_poly)*100

  # mapview::mapview(connected_subpop_poly) + mapview::mapview(isolated_subpop_poly, col.regions = "red")
  # 

  if(export_shp) {
    
    polys <- rbind(isolated_subpop_poly, connected_subpop_poly)
    
    return(list(fraction_aoo = fraction_aoo, subpop_poly = polys))
    
  } else {
    names(fraction_aoo) <- as.character(unique(XY$tax))
    return(fraction_aoo)
  }
  
  
 
  # 
  # #Creating the raster at the desired resolution (again I had problems with the projection...)
  # r <- sf::st_as_sf(raster::raster(raster::extent(XY_proj),
  #                resolution = cell_size * 1000,
  #                crs = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  # 
  # # Getting the cells by the species 
  # sf::st_intersects(XY_proj, r)
  # 
  # # Getting the cells within the buffers (probably occupied cells)  
  # sf::st_intersects(res_subpop_poly, r)
  
    
  # obtenir le numero de habitat patches (gid cells que se touche)

  # calculer numero de patches isolé/numero de pixels occupés (

}