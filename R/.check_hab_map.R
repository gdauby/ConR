#' @title  Internal function
#' 
#' @description Internal function to verify and process hab.map used in AOH.estimation and AOO.decline
#' 
#' @param hab.map raster, raster layer/stack or spatial polygons containing the
#'   habitat spatial information
#' @param hab.map.type logical, vector of same lenght of hab.map, if TRUE means hab.map is suitable for species, if FALSE unsuitable
#' 
#' @importFrom  terra rast
#' 
.check_hab_map <- function(hab.map, hab.map.type = NULL) {
  
  if (!any(class(hab.map) == "list"))
    hab.map <- list(hab.map)
  
  classes.hab.map <-
    data.frame(
      class_hab = unlist(lapply(hab.map, function(x)
        class(x)[1])),
      rast = grepl("Raster", unlist(lapply(hab.map, function(x)
        class(x)[1]))),
      sf = grepl("sf", unlist(lapply(hab.map, function(x)
        class(x)[1]))),
      nbe.layers = NA
    )
  
  if (!any(apply(classes.hab.map[, 2:ncol(classes.hab.map)], MARGIN = 2, any))) {
    stop("hab.map must be class of raster, RasterBrick, RasterStack, SpatRaster or sf")
    
  }
  
  if (any(classes.hab.map$rast))  {
    
    hab.map[[which(classes.hab.map$rast)]] <-
      terra::rast(hab.map[[which(classes.hab.map$rast)]])
    
    classes.hab.map[which(classes.hab.map$rast), "nbe.layers"] <- 
      terra::nlyr(hab.map[[which(classes.hab.map$rast)]])
    
  }
  
  classes.hab.map[which(!classes.hab.map$rast), "nbe.layers"] <- 
    1
  
  if (length(hab.map.type) != length(hab.map)) {
    
    stop("hab.map.type should be of equal size of the number of hab.map provided")
    
  }
  
  
  if (is.null(hab.map.type)) {
    
    message("hab.map.type not provided, all hab.map provided are considered as suitable")
    hab.map.type <- rep(TRUE, length(hab.map))
    
  }
  
  return(
    list(
      hab.map = hab.map,
      classes.hab.map = classes.hab.map,
      hab.map.type = hab.map.type
    )
  )
}