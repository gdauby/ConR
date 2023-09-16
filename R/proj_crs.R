
#' @title Get the projection system
#'
#' @description get class "CRS" of coordinate reference system
#' 
#' @author Gilles Dauby & Renato A. Ferreira de Lima
#' 
#' @param proj_type string or numeric
#' 
#' @import sf
#' 
#' @return `crs` class object
#' 
#' 
#' @examples 
#' ## By default, the output projection is Equal-Area Scalable Earth (EASE) Global
#' ## Global cylindrical equal-area projection
#' ## See https://epsg.io/54034
#' proj_crs(proj_type = "cea")
#' 
#' ## Another readily available projection is Antarctic Polar Stereographic
#' ## https://epsg.io/3031
#' proj_crs(proj_type = "Antarctic")
#' 
#' ## Otherwise, many projection can be retrieved from EPSG code
#' 
#' 
#' @export
#' 
proj_crs <- function(proj_type = "cea") {
  
  if (length(proj_type) > 1)
    stop("proj_type should be of length 1")
  
  if (is.character(proj_type)) {
    
    if (any(proj_type == c("cea", "Antarctic", "Africa_eac"))) {
      
      ## https://epsg.io/54034
      # Cylindrical Equal Area
      if(proj_type == "cea")
        proj_type <-
          "ESRI:54034"
      
      # https://epsg.io/3031
      if(proj_type == "Antarctic")
        proj_type <-
          'EPSG:3031'
      
      # https://epsg.io/102022
      if(proj_type == "Africa_eac")
        proj_type <-
          "ESRI:102022"
    }
    
    
    crs_proj <- sf::st_crs(proj_type)
    
  }
  
  if (is.numeric(proj_type))
    crs_proj <- sf::st_crs(proj_type)
  
  return(crs_proj)
  
  
}


# proj_crs <- function(proj_type = "cea", wkt = FALSE) {
#   
#   if (length(proj_type) > 1)
#     stop("proj_type should be of length 1")
#   
#   ## https://epsg.io/54032
#   # World Azimuthal Equidistant
#   # proj <-
#   #   "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#   
#   ## https://epsg.io/54002
#   # World Azimuthal Equidistant
#   # proj <-
#   #   "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#   
#   if(is.character(proj_type)) {
#     
#     match.arg(proj_type, choices = c("cea", "Antarctic", "Africa_eac"))
#     
#     # https://epsg.io/6933
#     if(proj_type == "cea")
#       proj <-
#         "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
#     
#     # https://epsg.io/3031
#     if(proj_type == "Antarctic")
#       proj <-
#         "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# 
#     # https://epsg.io/102022
#     if(proj_type == "Africa_eac")
#       proj <-
#         "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#     
#   } 
#   
#   if (is.numeric(proj_type)) {
#       
#     proj <- 
#       all_epsg[which(all_epsg$code == proj_type), "prj4"]
#     
#     if(length(proj) == 0)
#       stop("No projection found given proj_type")
#     
#     print(all_epsg[which(all_epsg$code == proj_type), ])
#     
#   }
# 
#   
#   if (rgdal::rgdal_extSoftVersion()[1] >= "3.0.0")  {
#     
#     # if(proj_type == "cea")
#     #   wkt_crs <- 
#     #     'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["standard_parallel_1",0.0],UNIT["Meter",1.0]]'
# 
#     if(!is.numeric(proj_type)) {
#       if(proj_type == "cea")
#         # wkt_crs <- 
#         #   'EPSG:6933'
#         wkt_crs <- 
#           'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["standard_parallel_1",0],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["ESRI","54034"]]'
#       
#       if(proj_type == "Antarctic")
#         wkt_crs <- 
#           'EPSG:3031'
#       
#       if(proj_type == "Africa_eac")
#         wkt_crs <- 
#           'PROJCS["Africa_Albers_Equal_Area_Conic",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["longitude_of_center",25],PARAMETER["Standard_Parallel_1",20],PARAMETER["Standard_Parallel_2",-23],PARAMETER["latitude_of_center",0],UNIT["Meter",1],AUTHORITY["EPSG","102022"]]'
#       
#     } else {
#       
#       wkt_crs <- paste0("EPSG:", proj_type)
#       
#     }
#     
#     crs_proj <- 
#       sp::CRS(projargs = proj, 
#               SRS_string = wkt_crs, 
#               doCheckCRSArgs = TRUE)
#     
#   }
#   
#   if (rgdal::rgdal_extSoftVersion()[1] < "3.0.0")
#     crs_proj <-
#       sp::CRS(projargs = proj, doCheckCRSArgs = FALSE)
#   
#   if (wkt) {
#     
#     if (rgdal::rgdal_extSoftVersion()[1] >="3.0.0") {
#       
#       return(wkt_crs)
#     } else {
#       
#       stop("Rgdal version >= 3.0.0 mandatory, please update your version")
#       
#     }
#     
#   } else {
#     
#     return(crs_proj)
#     
#   }
#   
#   
# }


