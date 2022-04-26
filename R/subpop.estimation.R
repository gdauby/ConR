#' @title Internal function
#'
#' @description Subpopulations estimation
#' 
#' @param XY data.frame with first two columns are projected coordinates (numeric)
#' @param Resol_sub_pop numeric. Defines the radius of the circles
#'   around each occurrence, in kilometres.
#' @param export_shp logical
#' @param proj_type
#'
#' @importFrom utils packageVersion
#' @import sf
#' 
subpop.estimation <- function(XY,
                         Resol_sub_pop, 
                         export_shp = FALSE,
                         proj_type) {
  
  if(utils::packageVersion("sf") < '0.9.0')
    stop('A more recent version of the "sf" package is needed, please update this package')
  
  # sp::coordinates(XY) <- c(2, 1)
  # 
  # XY_sp_proj_buff <-
  #   rgeos::gBuffer(XY, 
  #                  width = Resol_sub_pop * 1000, 
  #                  id = 1)
  
  #### GILLES: SOME CODES FROM OTHER FUNCTIONS TO USE sf INSTEAD OF sp and rgeos
  
  # centroid <- st_as_sf(data.frame(mean_ddlon = mean(XY$ddlon), 
  #                     mean_ddlat = mean(XY$ddlat)), 
  #          coords = c("mean_ddlon", "mean_ddlat"))
  # st_crs(centroid) <- 4326
  # centroid_proj <- st_transform(centroid, crs = proj_crs(3786))
  # centroid_proj <- st_coordinates(centroid_proj)
  # centroid_proj <- 
  #   rbind(centroid_proj, centroid_proj + Resol_sub_pop * 1000)
  # centroid_proj <- st_as_sf(as.data.frame(centroid_proj), 
  #                      coords = c("X", "Y"))
  # st_crs(centroid_proj) <- proj_crs(3786)
  # centroid_added <- st_transform(centroid_proj, crs = 4326)
  
  # st_distance(centroid_added)
  # 
  # dist(st_coordinates(centroid_added))
  
  
  points_sf <- sf::st_as_sf(XY, coords = c(1, 2))
  buff_sf <- sf::st_buffer(points_sf, dist = Resol_sub_pop * 1000)
  buff_sf <- sf::st_union(buff_sf) ### GILLES, WE NEED TO JOIN THE OVERLAPPING CIRCLES, RIGHT? PLEASE CHECK!
  # buff_sf_multipoly <- sf::st_cast(buff_sf, "MULTIPOLYGON")
  buff_sf <- sf::st_cast(buff_sf, "POLYGON")
  SubPopPoly <-
    sf::st_as_sf(data.frame(buff_sf))
  # SubPopPoly_multipoly <-
  #   sf::st_as_sf(data.frame(buff_sf_multipoly, tax = ifelse(any(colnames(XY) == "tax"), unique(XY$tax), "taxa1")))
  
  # XY_sf <- 
  #   st_as_sf(XY, coords = c(1, 2))
  # XY_sp_proj_buff <- 
  #   st_buffer(XY_sf, Resol_sub_pop * 1000)
  
  # SubPopPoly <- 
  #   sf::st_cast(as(XY_sp_proj_buff, "sf"), "POLYGON")

  sf::st_crs(SubPopPoly) <- proj_type
  
  NbeSubPop <- nrow(SubPopPoly)

  if (export_shp) { ## IF/ELSE ADDED BY RENATO
    #SubPopPoly <- as(SubPopPoly, "Spatial") # not transforming to sp anymore
    OUTPUT <- list(NbeSubPop, SubPopPoly)
    names(OUTPUT) <- c("number_subpop", "poly_subpop")
  
  } else {
    OUTPUT <- NbeSubPop
    names(OUTPUT) <- c("number_subpop")
  }
  
  return(OUTPUT)
}
