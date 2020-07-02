#' @title Internal function
#'
#' @description Subpopulations estimation
#' 
#' @param XY data.frame
#' @param Resol_sub_pop integer
#' 
#'
#' @importFrom rgeos gBuffer
#' @importFrom utils packageVersion
#' @importFrom sf st_transform sf_project st_crs st_cast
#' 
subpop.estimation <- function(XY,
                         Resol_sub_pop, 
                        proj_type) {
  
  if(utils::packageVersion("sf") < '0.9.0')
    stop('A more recent version of the "sf" package is needed, please update this package')
  
  sp::coordinates(XY) <- c(2, 1)
  
  XY_sp_proj_buff <-
    rgeos::gBuffer(XY, 
                   width = Resol_sub_pop * 1000, 
                   id = 1)
  
  
  # XY_sf <- 
  #   st_as_sf(XY, coords = c(1, 2))
  # XY_sp_proj_buff <- 
  #   st_buffer(XY_sf, Resol_sub_pop * 1000)
  
  SubPopPoly <- 
    sf::st_cast(as(XY_sp_proj_buff, "sf"), "POLYGON")

  sf::st_crs(SubPopPoly) <- proj_type
    
  # SubPopPoly <- 
  #   sf::st_transform(SubPopPoly, crs = 4326)
  
  
  NbeSubPop <- nrow(SubPopPoly)
  
  SubPopPoly <- as(SubPopPoly, "Spatial")
  
  OUTPUT <- list(NbeSubPop, SubPopPoly)
  names(OUTPUT) <- c("number_subpop", "poly_subpop")
  return(OUTPUT)
}