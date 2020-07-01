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
subpop.comp <- function(XY,
                         Resol_sub_pop) {
  
  if(utils::packageVersion("sf") < '0.9.0')
    stop('A more recent version of the "sf" package is needed, please update this package')
  
  projEAC <- .proj_crs()
  
  XY <- XY[, c(2, 1)]
  
  # coordEAC <-
  #   as.data.frame(matrix(unlist(
  #     rgdal::project(as.matrix(XY), proj = as.character(projEAC), inv = FALSE)
  #   ), ncol = 2))
  # rownames(coordEAC) <- seq(1, nrow(coordEAC), 1)
  
  XY_sf_proj <- 
    sf::sf_project(from = sf::st_crs(4326), 
                   to = 
                     sf::st_crs(projEAC), 
                   pts = XY)
  
  XY_sp <- as.data.frame(XY_sf_proj)
  
  sp::coordinates(XY_sp) <- c(1, 2)
  
  XY_sp_proj_buff <-
    rgeos::gBuffer(XY_sp, 
                   width = Resol_sub_pop * 1000, 
                   id = 1)
  
  XY_sf <- as(XY_sp_proj_buff, "sf")
  sf::st_crs(XY_sf) <- projEAC
  
  # if (rgdal::rgdal_extSoftVersion()[1] >= "3.0.0") 
  #   sf::st_crs(XY_sf) <- 4326
  #   
  #   sp::proj4string(XY_sp) <- sp::CRS(SRS_string = "EPSG:4326", doCheckCRSArgs = FALSE)
  # if (rgdal::rgdal_extSoftVersion()[1] < "3.0.0") 
  #   sp::proj4string(XY_sp) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs", doCheckCRSArgs = FALSE)
  
  #   XY_sp <- as(XY_sf, "Spatial")
  #   
  #   XY_sp_proj <- 
  #     sp::spTransform(XY_sp, CRSobj = projEAC)
  #   
  #   coordEAC <- as.data.frame(XY_sp_proj)
  # 
  # XY_sp_proj_no_proj <- XY_sp_proj
  # raster::crs(XY_sp_proj_no_proj) <- NA
  # 
  # 
  # 
  # XY_sf <- as(XY_sp_proj_buff, "sf")
  
  SubPopPoly <- sf::st_cast(XY_sf, "POLYGON")
  SubPopPoly <- 
    sf::st_transform(SubPopPoly, crs = 4326)
  
  NbeSubPop <- nrow(SubPopPoly)
  
  SubPopPoly <- as(SubPopPoly, "Spatial")
  
  
  # sp::proj4string(SubPopPoly) <- projEAC
  # if (rgdal::rgdal_extSoftVersion()[1] >= "3.0.0")
  #   SubPopPoly <- 
  #   sp::spTransform(x = SubPopPoly, sp::CRS(SRS_string='epsg:4326', doCheckCRSArgs=TRUE))
  # if (rgdal::rgdal_extSoftVersion()[1] < "3.0.0") 
  #   SubPopPoly <- sp::spTransform(x = SubPopPoly, sp::CRS("+proj=longlat +datum=WGS84 +no_defs", doCheckCRSArgs=TRUE))
  
  # if (utils::packageVersion("sp") >= "1.3.3") poly_masked@proj4string <- sp::CRS(SRS_string='epsg:4326')
  # if (utils::packageVersion("sp") < "1.3.3") poly_masked@proj4string <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  # p2 <-
  #   rgeos::readWKT(paste("POINT(", 
  #                        mean(unique(coordEAC)[1, 1]), " ", 
  #                        mean(unique(coordEAC)[1, 2]), ")", sep =
  #                          ""))
  
  
  
  # p2_Buffered1 <-
  #   rgeos::gBuffer(p2, width = Resol_sub_pop * 1000, id = 1)
  # if (nrow(unique(coordEAC)) > 1) {
  #   for (LL in 2:nrow(unique(coordEAC))) {
  #     p2 <-
  #       rgeos::readWKT(paste("POINT(", mean(unique(coordEAC)[LL, 1]), " ", 
  #                            mean(unique(coordEAC)[LL, 2]), ")", sep =
  #                              ""))
  #     p2_Buffered <-
  #       rgeos::gBuffer(p2, width = Resol_sub_pop * 1000, id = LL)
  #     p2_Buffered1 <- rgeos::gUnion(p2_Buffered1, p2_Buffered)
  #   }
  # }
  # 
  # splited_pol <-
  #   lapply(p2_Buffered1@polygons, slot, "Polygons")[[1]]
  # 
  # NbeSubPop <- length(splited_pol)
  # 
  # SubPopPoly <-
  #   sp::SpatialPolygons(
  #     Srl = list(p2_Buffered1@polygons[[1]]),
  #     pO = as.integer(1),
  #     proj4string = projEAC
  #   )
  
  # SubPopPoly <-
  #   sp::spTransform(SubPopPoly,
  #                   sp::CRS("+proj=longlat +datum=WGS84"))
  
  
  OUTPUT <- list(NbeSubPop, SubPopPoly)
  names(OUTPUT) <- c("Number of subpopulation", "subpop.poly")
  return(OUTPUT)
}