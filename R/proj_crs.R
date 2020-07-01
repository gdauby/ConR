
#' @title Get Coordinate Reference System
#'
#' @description get proj CRS
#' 
#' @param proj_type string or numeric
#' 
#' @importFrom utils packageVersion
#' @importFrom rgdal showWKT rgdal_extSoftVersion make_EPSG
#' 
#' @export
proj_crs <- function(proj_type = "cea") {
  
  ## https://epsg.io/54032
  # World Azimuthal Equidistant
  # proj <-
  #   "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  ## https://epsg.io/54002
  # World Azimuthal Equidistant
  # proj <-
  #   "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  
  # https://spatialreference.org/ref/sr-org/world-cylindrical-equal-area/
  # Equal Area Cylindrical
  # proj <-
  #   "+proj=cea +lat_ts=0.0  +lon_0=0.0 +ellps=GRS80 +k_0=1.0 +x_0=0.0 +y_0=0.0"
  # if(proj_type == "cea")
  #   proj <-
  #     "+proj=cea +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  # https://epsg.io/6933
  
  if(!is.numeric(proj_type)) {
    
    if(proj_type == "cea")
      proj <-
        "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    
    if(proj_type == "Antarctic")
      proj <-
        "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
  } else {
      
    all_epsg <-
        rgdal::make_EPSG()
    
    proj <- 
      all_epsg[which(all_epsg$code == proj_type), "prj4"]
    
    if(length(proj) == 0)
      stop("No projection found given proj_type")
  }
   
  # if(proj_type == "cea")

  
  
  # https://spatialreference.org/ref/epsg/3031/
  # Equal Area Cylindrical
  # proj <-
  #   "+proj=cea +lat_ts=0.0  +lon_0=0.0 +ellps=GRS80 +k_0=1.0 +x_0=0.0 +y_0=0.0"
  
  
  
  ## https://epsg.io/102022
  ## Africa Albers Equal Area Conic
  # proj <- 
  #   "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  # "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  if (rgdal::rgdal_extSoftVersion()[1] >= "3.0.0")  {
    
    # if(proj_type == "cea")
    #   wkt_crs <- 
    #     'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["standard_parallel_1",0.0],UNIT["Meter",1.0]]'

    if(!is.numeric(proj_type)) {
      if(proj_type == "cea")
        wkt_crs <- 
          'EPSG:6933'
      
      if(proj_type == "Antarctic")
        wkt_crs <- 
          'EPSG:3031'
    } else {
      
      wkt_crs <- paste0("EPSG:", proj_type)
      
    }
    
    # wkt_crs <-
    #   rgdal::showWKT(
    #     proj
    #   )
    
    # rgdal::make_EPSG() %>% 
    #   dplyr::as_tibble() %>% 
    #   dplyr::filter(code == 3031)
    # 
    # 
    # rgdal::make_EPSG() %>% 
    #   dplyr::as_tibble() %>% 
    #   dplyr::filter(grepl("+proj=cea", prj4)) %>% 
    #   dplyr::pull(prj4)
    
    crs_proj <- 
      sp::CRS(projargs = proj, 
              SRS_string = wkt_crs, doCheckCRSArgs = TRUE)
    
    # crs_proj <- 
    #   sp::CRS(projargs = proj, 
    #                     SRS_string = wkt_crs, doCheckCRSArgs = TRUE)
  }
  
  if (rgdal::rgdal_extSoftVersion()[1] < "3.0.0")
    crs_proj <-
      sp::CRS(projargs = proj, doCheckCRSArgs = FALSE)
  
  return(crs_proj)
}
