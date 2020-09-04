#' @title Internal function
#'
#' @description Subpopulations estimation
#' 
#' @param XY data.frame
#' @param Resol_sub_pop integer
#' @param export_shp logical
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
  points_sf <- sf::st_as_sf(XY[,1:3], coords = c("ddlon", "ddlat"))
  buff_sf <- sf::st_buffer(points_sf, dist = Resol_sub_pop * 1000)
  buff_sf <- sf::st_union(buff_sf) ### GILLES, WE NEED TO JOIN THE OVERLAPPING CIRCLES, RIGHT? PLEASE CHECK!
  buff_sf <- sf::st_cast(buff_sf, "POLYGON")
  SubPopPoly <-
    sf::st_as_sf(data.frame(buff_sf, tax = unique(XY$tax)))
  
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
