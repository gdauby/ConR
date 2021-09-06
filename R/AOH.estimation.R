#' @title  Terrestrial Area of Habitat
#' 
#' @description The function compute the amount of 'habitat' within the EOO of
#'   each species, which may represent habitat itself or any other proxy, value or
#'   category the user may find relevant (e.g. area covered by protected areas,
#'   altitude, etc...). If two or more time intervals are provided, the function
#'   also returns the habitat change between intervals. Currently, the type of
#'   habitat maps objects currently accepted are sp, sf and Raster*. In
#'   addition, habitat change is only computed for rasters (e.g. RasterStack
#'   object).
#' 
#' @param XY \code{"dataframe"} see Details
#' @param hab.map raster, raster layer/stack or spatial polygons containing the
#'   habitat spatial information
#' @param hab.class classes of values in ```hab.map``` to be considered as 'habitat'.
#' @param years numeric. Time interval between the first and last ```hab.map```
#'   if more than one raster is provided (e.g. if ```hab.map``` is a RasterStack object).
#' @param output_raster the output format in the case raster are continuous
#'   variables and not categories. Default to "summary".
#' @param ID_shape string. Indicate the field name in ```hab.map``` data frame
#'   containing the ID or name of each polygon.
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param proj_type character string, numeric or object of CRS class. Default to "cea".
#' @param parallel logical. Should computing run in parallel? Default to FALSE.
#' @param NbeCores integer. The number of cores for parallel execution. Default to 2.
#' @param show_progress logical. whether a bar showing progress in computation
#'   should be shown. Default to TRUE.
#' @param exclude.area a logical, if TRUE, areas outside of \code{country_map}
#' are cropped of \code{SpatialPolygons} used for calculating EOO. By default
#' is FALSE
#' @param country_map a \code{SpatialPolygonsDataFrame} or
#' \code{SpatialPolygons} showing for example countries or continent borders.
#' This shapefile will be used for cropping the \code{SpatialPolygons}l if
#' exclude.area is TRUE
#' 
#' @details If multiple rasters are provided, the last raster is taken as the
#'   most recent one to calculate current habitat proportion and area.
#'   Similarly, the first and last rasters are used to calculate habitat change.
#'   In this case the user should provide the interval in years between the two
#'   rasters. If this interval is not provided the function assumes the number
#'   of rasters provided minus one as the 'default' time interval.
#'
#'   For faster processing of maps, the function automatically crops the EOO
#'   polygons to the extent of the ```hab.map``` provided. Then, for each taxon,
#'   the ```hab.map``` is cropped to the extent of the taxon EOO polygon. If the
#'   taxon EOO and the raster provided don't intersect the functions returns NA.
#'
#'   If the size of the EOO is too small in respect to the raster resolution
#'   (small enough to contain an entire pixel) or if the EOO polygon is
#'   elongated (perimeter/area relationship > 2), the function returns the
#'   estimates for all pixels touching the species EOO boundaries and not only
#'   the pixels within the EOO. However, if the EOO polygon is small enough to
#'   fit entirely within a pixel, metrics will refer to this single pixel. So,
#'   make sure that the resolution of the rasters provided make sense to the
#'   range of the taxa being assessed.
#'
#' @return If ```hab.map``` is a *Raster object with categorical variables, the
#'   function returns a data frame containing the taxon name (tax), the
#'   percentage of habitat within the taxon EOO (prop.EOO) and the approximate
#'   area that this percentage represents. If habitat change is also estimated
#'   the function returns additional columns, namely: 
#'   - the total percentage loss and recover of habitat within the EOO (loss,
#'   recover);
#'   - the percentage loss and recover relative to the amount of habitat and
#'   non-habitat within the EOO at the beginning of the time interval (rel.loss,
#'   rel.recover);
#'   - the annual rate of habitat loss and recover (rate.loss, rate.recover,
#'   i.e. rel.loss and rel.recover divided by ```years```);
#'   - A habitat quality index which is obtained by averaging values of quality
#'   to pixels according to their transition between the first and last year of
#'   the times series, namely habitat to habitat = 2, non-habitat to habitat =
#'   1, non-habitat to non-habitat = 0, and habitat to non-habitat = -1
#'   (hab.quality, hab.quality.lo, hab.quality.hi). The index thus vary from -1
#'   (i.e., all pixels were habitat in time one and all was lost in time two) until 2
#'   (i.e., all pixels were habitat in time one and they remained so in time 2).
#'    
#'   If ```hab.map``` is a *Raster object with continuous variables, the
#'   function returns a data frame containing the taxon name (tax), and, by
#'   default, the summary statistics of all pixels within the taxon EOO. If
#'   argument ```output_raster``` is "prop.table" the function return the proportion
#'   of each pixel value within the taxon EOO and if it is "area.table" it returns 
#'   the approximate area of each pixel value within the taxon EOO.
#'   
#'   If ```hab.map``` is a 'sp' or 'sf' object with multiple polygons, the
#'   function returns a data frame containing the taxon name (tax), the number
#'   of unique polygons (numb.polys), the percentage of area (prop.EOO) and the
#'   total area of the polygons within the taxon EOO.
#'   
#'   If ```export_shp``` is TRUE, the function also returns the transition rasters
#'   used to calculated 'hab.quality' or the cropped polygons in the case of 
#'   polygons. But note that for rasters, this include additional steps that can
#'   be time and memory consuming.
#'   
#' @examples No examples available yet.
#' 
#' @importFrom  raster extent crop crs res extract getValues mask area cellStats rasterFromXYZ coordinates extend brick
#' @import sf
#' 
#' @export AOH.estimation
AOH.estimation <- function(XY,
                     show_progress = TRUE,
                     proj_type = "cea",
                     mode = "spheroid",
                     hab.map = NULL,
                     hab.class = NULL,
                     years = NULL,
                     country_map = NULL,
                     exclude.area = T,
                     ID_shape = NULL,
                     output_raster = NULL,
                     ...) {
  
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  } else {
    
    if(any(grepl('sf', class(country_map))))
      country_map <- 
        suppressWarnings(as(country_map, "Spatial"))
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
    country_map <- 
      as(country_map, "sf")
  }
  
  proj_type_ <- proj_crs(proj_type = proj_type)
  
  EOO.res <- EOO.computing(XY = XY, export_shp = TRUE, mode = mode, exclude.area = FALSE)
  
  EOO.shp <- EOO.res$spatial
  
  EOO.shp.proj <- 
    st_transform(EOO.shp, proj_type_)
  country_map.proj <- 
    st_transform(country_map, proj_type_)
  
  if (exclude.area) {
    
    EOO.shp.proj <-
      suppressWarnings(suppressMessages(
        sf::st_intersection(EOO.shp.proj, st_union(country_map.proj))
      ))
    
    # sf::st_crs(p1) <-
    #   4326
    # 
    # if(length(p1) == 0) {
    #   warning("After excluding areas, the convex hull is empty. EOO is set as NA.")
    #   
    # }
  }
  
  ### checking arguments validity
  if (is.null(hab.map))
    stop("Please provide hab.map (a sf polygon object or a raster)")
  
  if (!is.null(hab.map)) {
    
    if (!any(grepl("Raster", class(hab.map))) & !any(grepl("sf", class(hab.map)))) {
      
      stop("hab.map must be class of raster, RasterBrick, RasterStack or sf")
      
    }
  }
  
  if (!is.null(hab.map) & !is.null(hab.class)) {
    
    if (any(grepl("Raster", class(hab.map)))) {
      
      if (dim(hab.map)[3] != length(hab.class)) {
        
        stop ("The number of layers of hab.map is different to the number hab.class values. It should be identical")
        
      }
      
      if (any(hab.map@data@isfactor) & is.factor(hab.class)) {
        
        stop ("hab.map values are factors while hab.class are not. They must be of same types")
        
      }
    }
  }
  
  if (!any(class(hab.map) == "list"))
    hab.map <- list(hab.map)
  
  
  EOO.res.all <- data.frame(species = row.names(EOO.res$results),
                            EOO = EOO.res$results$EOO)
  
  # areas_hab.map <- vector('list', length(hab.map))
  AOH.poly <- vector('list', length(hab.map))
  for (i in 1:length(hab.map)) {
    
    hab.map.selected <- hab.map[[i]]
    
    if (any(grepl("Raster", class(hab.map.selected)))) {
      
      if (as.character(crs(hab.map.selected)) != proj_type_@projargs) {
        
        hab.map.selected <- raster::projectRaster(hab.map.selected, crs = proj_type_)
        
      }
      
      ## masking input raster given hab.class
      if (!is.null(hab.class)) {
        
        hab.map_masked <- vector('list', dim(hab.map.selected)[3])
        for (i in 1:dim(hab.map.selected)[3]) {
          
          hab.map_masked[[i]] <-
            raster::calc(hab.map.selected[[i]], function(x) ifelse(x < hab.class[i], NA, 1))
          
        }
        
        hab.map_masked <- raster::stack(hab.map_masked)
        
      } else {
        
        hab.map_masked <-
          raster::calc(hab.map.selected, function(x) ifelse(!is.na(x), 1, NA))
        
      }
      
      hab.map_masked_crop <-
        raster::crop(hab.map_masked, EOO.shp.proj) # cropping raster to the extent of the polygon
      
      # values(hab.map_masked_crop) <- 
      #   1:length(hab.map_masked_crop@data@values)
      
      # intersect_rast_EOO <- raster::extract(hab.map_masked_crop, EOO.shp.proj)
      
      # hab.map_masked_crop[[1]]
      
      hab.map_poly <- 
        st_as_sf(raster::rasterToPolygons(hab.map_masked_crop, dissolve = TRUE))
      
      EOO.shp.proj.hab.map <- st_intersection(EOO.shp.proj, hab.map_poly)
      
      AOH.poly[[i]] <- 
        EOO.shp.proj.hab.map
      
      res. <- data.frame(species = EOO.shp.proj.hab.map$tax,
                         AOH = as.numeric(st_area(EOO.shp.proj.hab.map)) / 1000000)
      
      colnames(res.)[2] <-
        paste0("hab.map.", i)
      
      EOO.res.all <- 
        merge(EOO.res.all, res.,
              by = "species", all.x = T)
      
      # EOO.res.all <- data.frame(species = row.names(EOO.res$results),
      #                           EOO = EOO.res$results$EOO,
      #                           AOH = as.numeric(st_area(EOO.shp.proj.hab.map)) / 1000000
      # )
      
      
    }
    
  }
  
  return(list(AOH = EOO.res.all,
              AOH.poly = AOH.poly))
  
}

# AOH.estimation <- function(EOO.poly,
#                         hab.map = NULL,
#                         hab.class = NULL,
#                         years = NULL,
#                         export_shp = FALSE,
#                         plot = FALSE,
#                         output_raster = "summary",
#                         ID_shape = NULL,
#                         mode = "spheroid",
#                         proj_type = "cea",
#                         parallel = FALSE,
#                         NbeCores = 2,
#                         show_progress = TRUE
# ) {
# 
#   # Initial checks
#   if (!any(grepl("SpatialPolygon|sf|Raster", class(hab.map))))
#     stop("Please provide habitat maps as sp, sf or Raster objects")
#   
#   if (!"tax" %in% names(EOO.poly))
#     stop("Please provide taxa polygons with a data frame containing taxa names in a column 'tax'")
#     
#   if (any(grepl("SpatialPolygon|sf", class(hab.map))) & is.null(ID_shape))
#     stop("Please provide the field name containing the ID of the polygons")
#   
#   if (any(grepl("Raster", class(hab.map))) & is.null(hab.class) & !output_raster %in% c("summary", "prop.table", "area.table"))
#     stop("Please provide one of the following output types: summary, prop.table, area.table")
# 
#   # Defining, extracting and transforming projections
#   proj_type_ <- proj_crs(proj_type = proj_type)
# 
#   if (is.na(sf::st_crs(hab.map)) ) {
#     stop("Please provide a valid coordinate system for the habitat map")
#   }
#   
#   # proj_hab <- raster::crs(hab.map)
#   
#   if (any(grepl("Raster", class(hab.map)))) {
#     
#     if (st_crs(crs(hab.map))$input != proj_type_@projargs) {
#       
#       hab.map <- raster::projectRaster(hab.map, crs = proj_type_)
#       
#     }
#     
#     # ext <- raster::extent(sf::st_bbox(EOO.poly)[c(1, 3, 2, 4)])
#     # hab.map_masked_crop <-
#     #   raster::crop(hab.map, ext) # cropping raster to the extent of the polygon
#     
#   }
#   
#   if(any(grepl("SpatialPolygon", class(hab.map))))
#     hab.map <- sf::st_as_sf(hab.map)
#   
#   
#   if (any(grepl("sf", class(hab.map)))) {
#     
#     if (st_crs(hab.map)$input != proj_type_@projargs) {
#       
#       hab.map <- st_transform(hab.map, crs = proj_type_)
#       
#     }
#     
#   }
#   
#   # proj_hab <- sf::st_crs(hab.map)
#   
#   #### GILLES: EOO.poly is supposed to be the output of EOO.comp and if so should never have empty crs. Useful? 
#   ### RESPONSE: that's safer
#   if (is.na(sf::st_crs(EOO.poly))) {
#     sf::st_crs(EOO.poly) <- 4326
#     warning("No coordinate system for EOO polygons: assuming WSG84")
#   }
#   
#   
#   EOO.poly <- 
#     sf::st_transform(EOO.poly, crs = proj_type_)
#   
#   # Croping EOO shapefiles to the habitat map extent (to speed up computational time)
#   EOO.poly.crop <- suppressMessages(suppressWarnings(
#     sf::st_intersection(EOO.poly, sf::st_as_sfc(sf::st_bbox(hab.map)))))
#   
#   ### GILLES: it was EOO.poly$geometry, but the output of EOO.comp is 'geom' and not 'geometry'. Dont know why, perhaps different classes
#   if (length(EOO.poly$geom) > length(EOO.poly.crop$geom))
#     warning(paste0("The EOO of ",
#                    length(EOO.poly$geom) - length(EOO.poly.crop$geom),
#                    " species is empty after croping them to the habitat map extent"))
#   
#   # Getting main descriptors for the EOO polygons
#   EOO.poly.area <- as.numeric(sf::st_area(EOO.poly.crop)) / 1000000 # area in km2
#   EOO.poly.peri <- as.numeric(lwgeom::st_perimeter(EOO.poly.crop)) / 1000 # perimeter in km
#   
#   # EOO.poly.proj <- sf::st_transform(EOO.poly.crop, crs = proj_crs)
#   # EOO.poly.area <- as.numeric(sf::st_area(EOO.poly.proj)) / 1000000 # area in km2
# 
#   # EOO.poly.peri <- as.numeric(lwgeom::st_perimeter(EOO.poly.proj)) / 1000 # perimeter in km
#   per.area <- EOO.poly.peri / # perimeter-area relationship (to detect alongated polygons)
#                 EOO.poly.area
#   
#   EOO.poly.name <- EOO.poly.crop
#   sf::st_geometry(EOO.poly.name) <- NULL
#   EOO.poly.name <- EOO.poly.name[,"tax"]
#   
#   # Getting raster values/categories if habitat classes are not provided
#   if (is.null(hab.class) & any(grepl("Raster", class(hab.map))))
#     classes <- seq(hab.map[[1]]@data@min, hab.map[[1]]@data@max, 1)
# 
#   ### EXTRACTION PER SE ###
#   #Preparing to parallel computing
#   if (parallel) {
#     cl <- snow::makeSOCKcluster(NbeCores)
#     doSNOW::registerDoSNOW(cl)
# 
#     message('Parallel running with ',
#             NbeCores, ' cores')
#     `%d%` <- foreach::`%dopar%`
#   } else{
#     `%d%` <- foreach::`%do%`
#   }
# 
#   x <- NULL
#   if (show_progress) {
#     pb <-
#       txtProgressBar(
#         min = 0,
#         max = nrow(EOO.poly.crop),
#         style = 3
#       )
# 
#     progress <- function(n)
#       setTxtProgressBar(pb, n)
#     opts <- list(progress = progress)
#   } else{
#     opts <- NULL
#   }
# 
#   #Starting the foreach loops
#   output <-
#     foreach::foreach(
#       
#       x = 1:nrow(EOO.poly.crop),
#       #.packages=c("raster","sf"),
#       .options.snow = opts
#       
#     ) %d% {
#       
#       if (!parallel & show_progress)
#         setTxtProgressBar(pb, x)
#       
#       # Methods is hab.map are sp or sf objects
#       if (any(grepl("sf", class(hab.map)))) {
#         
#         crop1 <-
#           sf::st_intersection(hab.map, EOO.poly.crop[x, ]) # cropping to the extent of the polygon
#         
#         crop.proj <- sf::st_transform(crop1, crs = proj_crs_)
#         
#         area.hab <-
#           as.numeric(sum(sf::st_area(crop.proj), na.rm = TRUE)) / 1000000
#         
#         area.poly <-
#           as.numeric(sf::st_area(sf::st_transform(EOO.poly.crop[x, ], crs = proj_crs_))) / 1000000
#         
#         tmp <-
#           c(area.poly - area.hab, area.hab)
#         
#         tmp1 <- crop.proj
#         sf::st_geometry(tmp1) <- NULL
#         tmp1 <- as.data.frame(tmp1)
#         n.polys <-
#           length(unique(as.character(tmp1[, ID_shape])))
#         hab.mat <- matrix(
#           c(NA,  n.polys, 100 * tmp / area.poly, tmp),
#           ncol = 3,
#           nrow = 2,
#           byrow = FALSE
#         )
#         row.names(hab.mat) <- c("non_habitat", "habitat")
#         colnames(hab.mat) <- c("numb.polys", "prop.EOO", "area.EOO")
#         res <- cbind.data.frame(tax = EOO.poly.name[x],
#                                 hab.mat,
#                                 stringsAsFactors = FALSE)["habitat", ]
#         
#         
#         if (export_shp) {
#           ## Tried to combine the two sf objects but could not make it with the df info
#           # toto <- c(st_geometry(crop1), st_geometry(EOO.poly.crop[x, ]))
#           # merge(toto, crop1)
#           # toto <- sf::st_multipolygon(st_geometry(eoo1), st_geometry(crop1))
#           # crop1 <- st_cast(crop1, "POLYGON")
#           # eoo1 <- st_cast(EOO.poly.crop[x,], "POLYGON")
#           
#           #### CHECK HERE ####
#           #crop1 <- sf::st_cast(crop1, "POLYGON")
#           res <- list(data.frame = res, poly = crop1)
#           
#           if (plot) {
#             par(
#               mar = c(3, 3, 2, 2),
#               las = 1,
#               tcl = -0.25,
#               mgp = c(2.5, 0.5, 0)
#             )
#             plot(sf::st_geometry(EOO.poly.crop[x,]),
#                  main = EOO.poly.name[x],
#                  bg = 0)
#             plot(sf::st_geometry(crop1), add = TRUE, col = 1)
#           }
#           
#         }
#       }
# 
#       # Methods if hab.map is a Raster* object
#       if (any(grepl("Raster", class(hab.map)))) {
#         
#         ext <- raster::extent(sf::st_bbox(EOO.poly.crop[x,])[c(1, 3, 2, 4)])
#         crop1 <-
#           raster::crop(hab.map, ext) # cropping raster to the extent of the polygon
#         
#         raster.res <- raster::res(crop1)
#         if (grepl("longlat", raster::crs(crop1))) {
#           
#           raster.res.km <- ((raster.res * 0.11) / 0.000001) / 1000
#           aprox.cell.size <- raster.res.km[1] * raster.res.km[2]
#           
#         } else {
#           
#           aprox.cell.size <- raster.res[1] * raster.res[2] / 1000000 #### area of pixel in km^2
#           
#         }
#         
#         if (EOO.poly.area[x] < aprox.cell.size |
#             (per.area[x] > 2 &
#              EOO.poly.area[x] < 10 * aprox.cell.size)) {
#           
#           transect <- sf::st_zm(sf::st_linestring(sf::st_coordinates(EOO.poly.crop[x,])))
#           transect <-
#             sf::st_sfc(transect, crs = sf::st_crs(EOO.poly.crop))
#           transect <- sf::st_sf(transect)
#           ext.transect <- raster::extract(crop1, transect,
#                                           along = TRUE, cellnumbers = TRUE)
#           tmp <-
#             ext.transect[[1]][!duplicated(ext.transect[[1]][, 1]),-1, drop = FALSE]
#           #Setting an empty object of the raster area (to ease code flow)
#           mask1 <- NULL
#           r.area <- NULL
#           
#         } else {
#           
#           mask1 <- raster::mask(crop1, mask = EOO.poly.crop[x,])
#           tmp  <-
#             raster::getValues(mask1) # much faster than raster::extract
#           
#           #Getting the approximated area of the masked raster (in km2)
#           mask.area <- raster::area(mask1[[1]], na.rm = TRUE)
#           r.area <- raster::cellStats(mask.area, 'sum')
#           
#         }
#         
#         #Removing NA pixels for calculations
#         if (class(tmp) %in% c("numeric", "integer", "logical"))
#           tmp <- tmp[!is.na(tmp)] # excluding pixels outside the EOO
#         
#         if (class(tmp) %in% c("matrix", "data.frame"))
#           tmp  <-
#           tmp[!is.na(tmp[, 1]), , drop = FALSE] # excluding pixels outside the EOO
#         
#         #Getting approximate raster area touching the EOO polygon (very small or elongated EOOs)
#         #This way of getting the area overestimates the real EOO area...
#         if (is.null(r.area) & class(tmp) %in% c("numeric", "integer", "logical"))
#           r.area <- length(tmp) * raster.res.km[1] * raster.res.km[2]
#         
#         if (is.null(r.area) &
#             class(tmp) %in% c("matrix", "data.frame"))
#           r.area <- dim(tmp)[1] * raster.res.km[1] * raster.res.km[2]
#         
#         #Methods if raster contains categorical variables (e.g. land-use classes)
#         if (!is.null(hab.class)) {
#           # Getting habitat quantities (land-use maps)
#           if (r.area == 0) {
#             
#             res <- matrix(
#               NA,
#               nrow = 1,
#               ncol = 2,
#               dimnames = list("habitat", c("prop.EOO", "area.EOO"))
#             )
#             
#           } else {
#             
#             if (is.logical(hab.class)) {
#               
#               if (!is.logical(tmp)) 
#                 stop("hab.class is logical but the provided raster values are not logical")
#               
#               hab.mat <- 
#                 data.frame(matrix(ncol = 3, nrow = 2))
#               rownames(hab.mat) <- c("non_habitat", "habitat")
#               colnames(hab.mat) <- c("n.pixs", "prop.EOO", "area.EOO")
#               hab.mat$n.pixs <- table(tmp)
#               hab.mat$prop.EOO <- hab.mat$n.pixs / length(tmp)
#               hab.mat$area.EOO <- (hab.mat$n.pixs * r.area) / 100
#               
#             } else {
#               
#               mask.dim <- dim(tmp)
#               pixs <- mask.dim[1]
#               last <- mask.dim[2]
#               hab.mat <-
#                 as.matrix(table(factor(
#                   tmp[, last] %in% hab.class,
#                   levels = c(FALSE, TRUE)
#                 )))
#               row.names(hab.mat) <- c("non_habitat", "habitat")
#               hab.mat <- cbind(hab.mat, 100 * hab.mat / pixs)
#               hab.mat <- cbind(hab.mat, (hab.mat[, 2] * r.area) / 100)
#               colnames(hab.mat) <- c("n.pixs", "prop.EOO", "area.EOO")
#               
#             }
#           }
#           
#           # Extra codes if more than one raster is provided
#           if (class(crop1) == "RasterBrick") {
#             if (r.area == 0) {
#               res <- cbind.data.frame(
#                 tax = EOO.poly.name[x],
#                 res,
#                 matrix(
#                   NA,
#                   nrow = 1,
#                   ncol = 9,
#                   dimnames = list(
#                     "habitat",
#                     c(
#                       "loss",
#                       "recover",
#                       "rel.loss",
#                       "rel.recover",
#                       "rate.loss",
#                       "rate.recover",
#                       "hab.quality",
#                       "hab.quality.lo",
#                       "hab.quality.hi"
#                     )
#                   )
#                 ),
#                 stringsAsFactors = FALSE
#               )
#               
#             } else {
#               # Getting the overal percentage of each transition
#               transit.mat <-
#                 100 * as.matrix(table(
#                   factor(tmp[, 1] %in% hab.class, levels = c(FALSE, TRUE)),
#                   factor(tmp[, last] %in% hab.class, levels = c(FALSE, TRUE))
#                 )) / pixs
#               colnames(transit.mat) <- c("non_habitat_t1", "habitat_t1")
#               row.names(transit.mat) <-
#                 c("non_habitat_t0", "habitat_t0")
#               hab.mat <-
#                 cbind(hab.mat, loss = c(transit.mat[1, 2], transit.mat[2, 1]))
#               hab.mat <-
#                 cbind(hab.mat, recover = c(transit.mat[2, 1], transit.mat[1, 2]))
#               
#               # % of habitat loss and recovery
#               hab_loss <-
#                 100 * transit.mat[2, 1] / sum(transit.mat[2,])
#               non_hab_loss <-
#                 100 * transit.mat[1, 2] / sum(transit.mat[1,])
#               hab_rec <- 100 * transit.mat[1, 2] / sum(transit.mat[, 2])
#               non_hab_rec <-
#                 100 * transit.mat[2, 1] / sum(transit.mat[, 1])
#               hab.mat <-
#                 cbind(hab.mat, rel.loss = c(non_hab_loss, hab_loss))
#               hab.mat <-
#                 cbind(hab.mat, rel.recover = c(non_hab_rec, hab_rec))
#               
#               # rate of loss and recovery
#               if (is.null(years)) {
#                 years <- dim(hab.map)[3] - 1
#                 warning("Time interval not provided and takes as the number of rasters minus one")
#               }
#               
#               non_hab_loss_rate <- non_hab_loss / years
#               hab_loss_rate <- hab_loss / years
#               non_hab_rec_rate <- non_hab_rec / years
#               hab_rec_rate <- hab_rec / years
#               hab.mat <-
#                 cbind(hab.mat,
#                       rate.loss = c(non_hab_loss_rate, hab_loss_rate))
#               hab.mat <-
#                 cbind(hab.mat,
#                       rate.recover = c(non_hab_rec_rate, hab_rec_rate))
#               
#               # defining habitat quality at time t+1
#               transit <- matrix(NA, ncol = 1, nrow = pixs)
#               transit[tmp[, 1] %in% hab.class &
#                         tmp[, last] %in% hab.class, 1] <- 2
#               transit[tmp[, 1] %in% hab.class &
#                         !tmp[, last] %in% hab.class, 1] <- -2
#               transit[!tmp[, 1] %in% hab.class &
#                         tmp[, last] %in% hab.class, 1] <- 1
#               transit[!tmp[, 1] %in% hab.class &
#                         !tmp[, last] %in% hab.class, 1] <- 0
#               if (dim(transit)[1] > 3) {
#                 mod <- stats::lm(base::jitter(transit, factor = 0.1) ~ 1)
#                 ci <- round(stats::confint(mod), 4)
#                 hab.mat <- cbind(
#                   hab.mat,
#                   hab.quality = round(as.numeric(stats::coef(mod)), 4),
#                   hab.quality.lo = ci[1],
#                   hab.quality.hi = ci[2]
#                 )
#                 
#               } else {
#                 hab.mat <- cbind(
#                   hab.mat,
#                   hab.quality = round(mean(transit), 4),
#                   hab.quality.lo = NA,
#                   hab.quality.hi = NA
#                 )
#                 
#               }
#               
#               res <-
#                 cbind.data.frame(tax = EOO.poly.name[x],
#                                  hab.mat[, -1],
#                                  stringsAsFactors = FALSE)["habitat", ]
#               
#               # location of decline and recover
#               #### CHECK HERE FOR A BETTER SOLUTION FOR SPECIES WIHOUT A MASK1 OBJECT ####
#               #THESE ARE SPECIES WITH VERY SMALL OR ELONGATED EOOs...
#               if (export_shp & !is.null(mask1)) {
#                 loc.loss <- data.frame(
#                   cell = 1:length(mask1[[1]]@data@values),
#                   raster::coordinates(mask1[[1]]),
#                   classes = NA_integer_
#                 )
#                 
#                 loc.loss$classes[!is.na(mask1[[1]]@data@values)] <-
#                   transit
#                 r <-
#                   raster::rasterFromXYZ(loc.loss[!is.na(loc.loss$classes), 2:4],
#                                         crs = raster::crs(mask1[[1]]))
#                 names(r) <- EOO.poly.name[x]
#                 
#                 if (plot) {
#                   par(
#                     mar = c(3, 3, 2, 2),
#                     las = 1,
#                     tcl = -0.25,
#                     mgp = c(2.5, 0.5, 0)
#                   )
#                   plot(
#                     r,
#                     col = c("red", "grey", "green", "darkgreen"),
#                     main = EOO.poly.name[x]
#                   )
#                   plot(sf::st_geometry(EOO.poly.crop[x, ]),
#                        add = TRUE,
#                        bg = 0)
#                 }
#                 
#                 res <- list(data.frame = res, raster = r)
#                 
#               }
#             }
#           }
#           
#           #Methods if raster contains continuous variables (e.g. altitude)
#         } else {
#           
#           #Getting summary stats for the variable in the EOO (quantitative maps)
#           if (class(tmp) %in% c("integer", "numeric", "logical")) {
#             
#             vals <- tmp # excluding pixels outside the EOO
#             pixs <- length(vals)
#             
#           }
#           
#           if (class(tmp) %in% c("matrix", "data.frame")) {
#             
#             last <- dim(tmp)[2]
#             vals <- tmp[, last] # excluding pixels outside the EOO
#             pixs <- length(vals)
#             
#           }
#           
#           sumario <-
#             matrix(
#               unclass(summary(vals)),
#               nrow = 1,
#               ncol = 6,
#               dimnames = list("", c(
#                 "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max"
#               ))
#             )
#           
#           ### CONFIDENCE INTERVALS WERE VERY UNINFORMATIVE (CLOSE TO THE MEAN), SO ARE NO LONGER RETURNED
#           # if(class(vals) %in% "integer") {
#           #   mod <- stats::glm(vals ~ 1, family = "poisson" )
#           #   ci <- suppressMessages(round(exp(stats::confint(mod)),4))
#           # }
#           # if(class(tmp) %in% "numeric") {
#           #   mod <- stats::lm(vals ~ 1)
#           #   ci <- suppressMessages(round(stats::confint(mod),4))
#           # }
#           
#           #Getting outputs in the format required by the user
#           if (output_raster %in% "summary") {
#             
#             res <- cbind.data.frame(tax = EOO.poly.name[x],
#                                     sumario, # ci,
#                                     stringsAsFactors = FALSE)
#             
#             
#           } else {
#             
#             #Tabulating pixels
#             tmp1 <- table(factor(vals, levels = classes))
#             prop <-
#               matrix(
#                 round(100 * unclass(tmp1) / pixs, 8),
#                 nrow = 1,
#                 ncol = length(classes),
#                 dimnames = list("", classes)
#               )
#             
#             if (output_raster == "prop.table")
#               res <- cbind.data.frame(tax = EOO.poly.name[x],
#                                       prop,
#                                       stringsAsFactors = FALSE)
#             
#             if (output_raster == "area.table") {
#               area <- round((prop * r.area) / 100, 3)
#               res <- cbind.data.frame(tax = EOO.poly.name[x],
#                                       area,
#                                       stringsAsFactors = FALSE)
#               
#             }
#             
#             
#           }
#         }
#       }
#       
#       res
#   }
#   
#   if(parallel) snow::stopCluster(cl)
#   if(show_progress) close(pb)
# 
#   #### THIS IF/ELSE NEEDS CHECKING WHEN 'export_shp' == TRUE ####
#   if (export_shp) {
#     #Getting the data frames
#     result <- 
#       sapply(output, '[[', "data.frame", simplify = FALSE, USE.NAMES = FALSE)
#     result <- 
#       do.call(rbind.data.frame, result)
#     rownames(result) <- NULL
#     
#     #Getting the rasters
#     if (any(grepl("Raster", class(hab.map)))) {
#       maps <- 
#         sapply(output, '[[', "raster", simplify = FALSE, USE.NAMES = FALSE)
#       
#       if (length(maps) > 1)
#         maps <- sapply(maps, function(x) raster::extend(x, hab.map))
#         maps <- raster::brick(maps)
#     }  
#     
#     if (grepl("sf", class(hab.map))) {
#       maps <- 
#         output[names(output) == "poly"]
#       if (length(maps) > 1)
#         maps <-
#           do.call("rbind", maps)
#         row.names(maps) <- NULL
#     }
#     
#   } else {
#     
#     result <- do.call(rbind.data.frame, output)
#     rownames(result) <- NULL
#     
#   }
#   
#   ## Merging the output with the entry data frame and returning
#   tax.df <- data.frame(order = 1:length(EOO.poly$tax), tax = EOO.poly$tax,
#                        stringsAsFactors = FALSE)
#   result1 <- merge(tax.df, result, 
#                    by = "tax", all = TRUE, sort = FALSE)
#   result1 <- result1[order(result1$order),]
#   result1 <- result1[,!names(result1) %in% "order", drop = FALSE]
#   
#   if (export_shp) {
#     return(list(data.frame = result1, raster.EOO = rasters))
#   } else {
#     return(result1)
#   }  
# }

# XY_sf <- st_transform(XY_sf, proj_type_)
# 
# data.frame(x = median(XY$ddlon), 
#            y = median(XY$ddlat))
# 
# elevation_rast <- 
#   elevatr::get_elev_raster(locations = data.frame(x = XY$ddlon, 
#                                                   y = XY$ddlat), z = 2, prj = sp::wkt(proj_type_))
# 
# elevation_rast <- crop(elevation_rast, XY_sf)
# 
# test <- calc(elevation_rast, fun=function(x){ x[x < 200 | x > 500] <- NA; return(x)})
# hab.map = test









