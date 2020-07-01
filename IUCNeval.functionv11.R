






# .crop.poly <- function(poly, crop) {
#   
#   # @importFrom spatstat as.owin area.owin union.owin setminus.owin 
#   
#   raster::crs(poly) <- NA
#   raster::crs(crop) <- NA
#   
#   # crs_ <- sp::CRS(SRS_string='epsg:4326')
#   # sf::st_crs(poly_sf) <- 4326
#   
#   
#   poly_sf <- methods::as(poly, "sf")
#   # crop_sf <- sf::st_combine(as(crop, "sf"))
#   crop_sf <- methods::as(crop, "sf")
#   
#   suppressMessages(suppressWarnings(diff_croped <- 
#                                       sf::st_intersection(crop_sf, poly_sf)))
#   
#   poly_masked <- methods::as(sf::st_union(diff_croped), "Spatial")
#   
#   if (rgdal::rgdal_extSoftVersion()[1] >= "3.0.0") 
#     poly_masked@proj4string <- sp::CRS(SRS_string='epsg:4326', doCheckCRSArgs = FALSE)
#   if (rgdal::rgdal_extSoftVersion()[1] < "3.0.0") 
#     poly_masked@proj4string <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs", doCheckCRSArgs = FALSE)
#   
#   # crs_crop <- raster::crs(crop)
#   # 
#   # raster::crs(poly) <- NA
#   # raster::crs(crop) <- NA
#   # 
#   # p1_owin <- spatstat::as.owin(poly)
#   # africa_owin <- spatstat::as.owin(crop)
#   # 
#   # if (round(spatstat::area.owin(spatstat::union.owin(p1_owin, africa_owin)), 3) != round(spatstat::area.owin(africa_owin), 3)) {
#   #   w <- spatstat::setminus.owin(p1_owin, africa_owin)
#   #   w2 <- spatstat::setminus.owin(p1_owin, w)
#   #   poly_masked <- as(w2, "SpatialPolygons")
#   #   
#   #   raster::crs(poly_masked) <- crs_crop
#   #   
#   # } else{
#   #   poly_masked <- poly
#   #   
#   # }
#   
#   EOO <-
#     round(suppressWarnings(geosphere::areaPolygon(poly_masked)) / 1000000, 1)
#   
#   return(list(EOO, poly_masked))
# }
































