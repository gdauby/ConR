#' Internal function
#' @noRd
check_hab_map <- function(hab.map, hab.map.type = NULL) {
  
  if (!any(class(hab.map) == "list"))
    hab.map <- list(hab.map)
  
  # classes.hab.map <-
  #   data.frame(
  #     class_hab = unlist(lapply(hab.map, function(x)
  #       class(x)[1])),
  #     rast = grepl("Raster", unlist(lapply(hab.map, function(x)
  #       class(x)[1]))),
  #     sf = grepl("sf", unlist(lapply(hab.map, function(x)
  #       class(x)[1]))),
  #     nbe.layers = NA
  #   )
  
  classes.hab.map <-
    data.frame(
      SpatRaster = unlist(lapply(hab.map, function(x) inherits(x, "SpatRaster"))),
      RasterLayer = unlist(lapply(hab.map, function(x) inherits(x, "RasterLayer"))),
      RasterStack = unlist(lapply(hab.map, function(x) inherits(x, "RasterStack"))),
      sf = unlist(lapply(hab.map, function(x) inherits(x, "sf"))),
      nbe.layers = as.numeric(0)
    )
  
  if (!any(apply(classes.hab.map[, 1:(ncol(classes.hab.map) - 1)], MARGIN = 2, any))) {
    stop("hab.map must be class of raster, RasterBrick, RasterStack, SpatRaster or sf")
    
  }
  
  if (any(classes.hab.map$RasterLayer)) {

    hab.map[[which(classes.hab.map$RasterLayer)]] <-
      terra::rast(hab.map[[which(classes.hab.map$RasterLayer)]])
    
    classes.hab.map$SpatRaster[which(classes.hab.map$RasterLayer)] <- TRUE
    
  }
    
  if (any(classes.hab.map$SpatRaster))
    classes.hab.map[which(classes.hab.map$SpatRaster), "nbe.layers"] <- 
    terra::nlyr(hab.map[[which(classes.hab.map$SpatRaster)]])
  
  if (any(classes.hab.map$sf))
    classes.hab.map[which(classes.hab.map$sf), "nbe.layers"] <-  1
  
  # classes.hab.map[which(!classes.hab.map$rast), "nbe.layers"] <- 
  #   1
  
  if (is.null(hab.map.type)) {
    
    message("hab.map.type not provided, all hab.map provided are considered as suitable")
    hab.map.type <- rep(TRUE, length(hab.map))
    
  }
  
  if (length(hab.map.type) != length(hab.map)) {
    
    stop("hab.map.type should be of equal size of the number of hab.map provided")
    
  }
  
  return(
    list(
      hab.map = hab.map,
      classes.hab.map = classes.hab.map,
      hab.map.type = hab.map.type
    )
  )
}