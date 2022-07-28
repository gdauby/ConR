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
#' 
#' @param XY `"dataframe"` see Details
#' @param show_progress logical. whether a bar showing progress in computation
#'   should be shown. Default to TRUE.
#' @param proj_type character string, numeric or object of CRS class. Default to "cea".
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param hab.map raster, raster layer/stack or spatial polygons containing the
#'   habitat spatial information
#' @param hab.map.type logical, vector of same lenght of hab.map, if TRUE means hab.map is suitable for species, if FALSE unsuitable
#' @param hab.class classes of values in ```hab.map``` to be considered as 'habitat'.
#' @param years numeric. Time interval between the first and last ```hab.map```
#'   if more than one raster is provided (e.g. if ```hab.map``` is a RasterStack object).
#' @param country_map a `SpatialPolygonsDataFrame` or
#' @param exclude.area a logical, if TRUE, areas outside of `country_map`
#' are cropped of `SpatialPolygons` used for calculating EOO. By default
#' is FALSE
#' @param parallel logical. Should computing run in parallel? Default to FALSE.
#' @param NbeCores integer. The number of cores for parallel execution. Default to 2.
#' @param simplifiy_poly logical whether the resulting polygon should be simplified using ms_simplify function of rmapshaper package
#' @param buffer logical


#' `SpatialPolygons` showing for example countries or continent borders.
#' This shapefile will be used for cropping the `SpatialPolygons`l if
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
#' @importFrom terra project crop vect rast app crs
#' @importFrom stars st_as_stars
#' @importFrom rmapshaper ms_simplify
#' @import sf
#' 
#' @export AOH.estimation
AOH.estimation <- function(XY,
                     show_progress = TRUE,
                     proj_type = "cea",
                     mode = "planar",
                     hab.map = NULL,
                     hab.map.type = NULL,
                     hab.class = NULL,
                     years = NULL,
                     country_map = NULL,
                     exclude.area = TRUE,
                     parallel = FALSE,
                     NbeCores = 2,
                     simplifiy_poly = TRUE,
                     buffer = NULL,
                     ...) {
  
  ### checking arguments validity
  if (is.null(hab.map))
    stop("Please provide hab.map (a sf polygon object or a raster or a SpatRaster)")
  
  
  ## converting to terra if raster
  # if (any(class(hab.map) == "RasterLayer")) {
  #   
  #   hab.map <- terra::rast(hab.map)[[1]]
  #   
  # }
  

  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  } else {
    
    # if(any(grepl('sf', class(country_map))))
    #   country_map <- 
    #     suppressWarnings(as(country_map, "Spatial"))
    
    
    country_map <-
      suppressWarnings(sf::st_buffer(country_map, dist = 0))
    
    # country_map <- 
    #   as(country_map, "sf")
  }
  
  proj_type_ <- proj_crs(proj_type = proj_type, wkt = T)
  
  EOO.res <- EOO.computing(XY = XY, export_shp = TRUE, mode = mode, exclude.area = FALSE) # , ...
  
  EOO.shp <- EOO.res$spatial
  
  EOO.shp.proj <- 
    st_transform(EOO.shp, proj_type_)
  
  if (!is.null(buffer))
    EOO.shp.proj <- st_buffer(EOO.shp.proj, buffer)
  
  country_map.proj <- 
    st_transform(country_map, proj_type_)
  
  if (!is.null(hab.map) & !is.null(hab.class)) {
    
    # if (any(grepl("SpatRaster", class(hab.map)))) {
    #   
    #   if (any(class(hab.class) == "data.frame")) {
    #     
    #     hab.class <- as.data.frame(hab.class)
    #     
    #     sp_with_eoo <- unique(EOO.shp.proj$tax)
    #     
    #     if (any(!sp_with_eoo %in% hab.class[,1])) {
    #       
    #       stop("Values in hab.class data.frame missing for ",sum(!sp_with_eoo %in% hab.class[,1]), " species")
    #       
    #     }
    #     
    #     if (any(as.numeric(hab.class[,3]) - as.numeric(hab.class[,2]) < 0)) {
    #       
    #       stop("The second column of hab.class should the minimum value while the third column should be the maximum. This is not the case.")
    #       
    #     }
    #     
    #     hab.class <- 
    #       hab.class[which(hab.class[,1] %in% sp_with_eoo),]
    #     
    #   }
    #   
    #   if (class(hab.class) == "numeric") {
    #     
    #     # if (dim(hab.map)[3] != length(hab.class)) {
    #     #   
    #     #   stop ("The number of layers of hab.map is different to the number hab.class values. It should be identical")
    #     #   
    #     #   if (any(hab.map@data@isfactor) & is.factor(hab.class)) {
    #     #     
    #     #     stop ("hab.map values are factors while hab.class are not. They must be of same types")
    #     #     
    #     #   }
    #     # }
    #   }
    # }
  }
  
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
  
  # if (any(class(hab.map) == "RasterStack"))
  #   hab.map <- raster::unstack(hab.map)
  
  if (!is.null(years)) {
    if (length(hab.map) != length(years))
      stop ("the number of years provided is different to the number of layers/polygon provided")
    
    names(hab.map) <- years
    
  }
  
  EOO.res.all <- EOO.res$results

  
  # loop across hab.map provided
  AOH.poly <- vector('list', length(hab.map) + 1)
  for (i in 1:length(hab.map)) {
    
    hab.map.selected <- hab.map[[i]]
    
    if (any(grepl("Raster", class(hab.map.selected)))) {
      
      if (as.character(terra::crs(hab.map.selected, describe = TRUE)$code) != unlist(strsplit(proj_type_, ":"))[2]) {
        
        # hab.map.selected <- raster::projectRaster(hab.map.selected, crs = proj_type_)
        
        hab.map.selected <- terra::project(hab.map.selected, proj_type_)
        
        # hab.map.selected.terra <- terra::rast(hab.map.selected)
        # 
        # 
        # microbenchmark::microbenchmark(raster = raster::projectRaster(hab.map.selected, crs = proj_type_),
        #                                terra = terra::project(hab.map.selected.terra, as.character(proj_type_)), 
        #                                times = 5)
        
        
      }
      
      hab.map.selected <-
        terra::crop(hab.map.selected, terra::vect(EOO.shp.proj))
        # raster::crop(hab.map.selected, st_bbox(EOO.shp.proj)) # cropping raster to the extent of the polygon
      
      ## masking input raster given hab.class
      if (!is.null(hab.class)) {
        
        if (any(class(hab.class) == "numeric")) {
          
          hab.map_masked <- vector('list', dim(hab.map.selected)[3])
          
          # system.time(
          #   for (j in 1:dim(hab.map.selected)[3]) {
              
          for (j in 1:length(hab.class)) {
            
            cropped.hab.map <- 
              terra::app(
                hab.map.selected,
                fun = function(x) {
                  x[x >=  as.numeric(hab.class[j])] <- NA
                  return(x)
                }
              )
            
            cropped.hab.map <- 
              terra::app(cropped.hab.map, function(x) {
              x[!is.na(x)] <- 1
              return(x)
            })
            
            hab.map_poly <- st_as_sf(stars::st_as_stars(cropped.hab.map), merge = TRUE)
            message("converting to valid polygons")
            hab.map_poly <- st_make_valid(hab.map_poly)
            message("combining subseted pollygons to multipolygon")
            hab.map_poly <- st_combine(hab.map_poly)
            hab.map_masked[[j]] <- st_as_sf(hab.map_poly)
            
            
            # split_vector <- rep(1:NbeCores, each = nrow(EOO.shp.proj[1:6,]) / NbeCores, length.out = nrow(EOO.shp.proj[1:6,]))
            
            # system.time(
            #   split_results <- split(EOO.shp.proj[1:6,], split_vector) %>%
            #     parallel::mclapply(function(x) st_intersection(x, hab.map_masked), mc.cores = 1)              
            # )

          }
          
              # cropped.hab.map <- 
              #   raster::calc(
              #     hab.map.selected,
              #     fun = function(x) {
              #       x[x <  min(as.numeric(hab.class))] <- NA
              #       return(x)
              #     }
              #   )
              
              # hab.map_poly <- 
              #   st_union(hab.map_masked[which(hab.map_masked$layer > as.numeric(hab.class[j])),])
              
              # hab.map_masked[[j]]  <- 
              #   st_as_sf(raster::rasterToPolygons(cropped.hab.map, dissolve = TRUE))
              
              
              
              # hab.map_masked <- st_as_sf(stars::st_as_stars(cropped.hab.map), merge = TRUE)
              # 
              # hab.map_masked <- st_make_valid(hab.map_masked)
              
              # hab.map_masked[[j]] <-
              #   raster::calc(hab.map.selected[[j]], function(x) ifelse(x < hab.class[j], NA, 1))
              
          #   }
          # )
        }
        
        if (any(class(hab.class) == "data.frame")) {
          
          hab.class <- as.data.frame(hab.class)
          
          cropped.hab.map <- 
            terra::app(hab.map.selected,
                     fun = function(x) {
                       x[x <  min(as.numeric(hab.class[, 2])) |
                           x > max(as.numeric(hab.class[, 3]))] <- NA
                       return(x)
                     })
          
          # cropped.hab.map <- 
          #   raster::calc(
          #     raster::raster(hab.map.selected),
          #     fun = function(x) {
          #       x[x <  min(as.numeric(hab.class[, 2]))|
          #           x > max(as.numeric(hab.class[, 3]))] <- NA
          #       return(x)
          #     }
          #   )
          
          # test <- terra::as.polygons(cropped.hab.map, trunc = TRUE, values = TRUE)
          
          # r2_pol <-
          #   terra::as.polygons(
          #     cropped.hab.map)
          # 
          # r2_pol_sf <- sf::st_as_sf(r2_pol)
          # 
          # r2_pol_sf <- suppressWarnings(sf::st_cast(r2_pol_sf, "POLYGON"))
          # row.names(r2_pol_sf) <- 1:nrow(r2_pol_sf)
          
          hab.map_poly <- st_as_sf(stars::st_as_stars(cropped.hab.map), merge = TRUE)
          
          if (simplifiy_poly) {
            
            message("Simplifying the hab.map multi polygon")
            hab.map_poly <- rmapshaper::ms_simplify(hab.map_poly)
            
          }
          
          message("converting to valid polygons")
          hab.map_masked <- st_make_valid(hab.map_poly)
          # message("combining subseted polygons to multipolygon")
          
          # hab.map_masked <- 
          #   st_as_sf(raster::rasterToPolygons(cropped.hab.map, dissolve = TRUE))
          
          # hab.map_masked <- vector('list', nrow(hab.class))
          # for (j in 1:nrow(hab.class)) {
          #   
          #   sp.selected <- hab.class[j,1]
          #   
          #   hab.class.selected <-
          #     hab.class[which(hab.class[, 1] == sp.selected), ]
          #   
          #   EOO.shp.selected <- 
          #     EOO.shp.proj[which(EOO.shp.proj$tax == sp.selected),]
          #   
          #   cropped.hab.map <- 
          #     calc(
          #     hab.map.selected,
          #     fun = function(x) {
          #       x[x < as.numeric(hab.class[j, 2]) |
          #           x > as.numeric(hab.class[j, 3])] <- NA
          #       return(x)
          #     }
          #   )
          #   
          #   hab.map_masked[[j]] <-
          #     crop(cropped.hab.map, EOO.shp.selected)
          #     
          # }
        }
        
        # hab.map_masked <- raster::stack(hab.map_masked)
        
      } else {
        
        cropped.hab.map <-
          terra::app(hab.map.selected, function(x) ifelse(!is.na(x), 1, NA))
        
        hab.map_masked <- st_make_valid(st_as_sf(stars::st_as_stars(cropped.hab.map), merge = TRUE))
        
      }
        
        if (any(class(hab.class) == "data.frame")) {
          
          if (parallel) {
            
            cl <- snow::makeSOCKcluster(NbeCores)
            doSNOW::registerDoSNOW(cl)
            
            message('Parallel running with ',
                    NbeCores, ' cores')
            
            `%d%` <- foreach::`%dopar%`
            
          } else {
            
            `%d%` <- foreach::`%do%`
            
          }
          
          x <- NULL
          if (show_progress) {
            pb <-
              txtProgressBar(min = 0,
                             max = nrow(hab.class),
                             style = 3)
            
            progress <- function(n)
              setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
          } else {
            opts <- NULL
          }
          
          output <-
            foreach::foreach(
              x = 1:nrow(hab.class),
              .options.snow = opts,
              .packages = 'sf'
            ) %d% {
              
              if (!parallel & show_progress)
                setTxtProgressBar(pb, x)
              
              hab.map_poly <- 
                hab.map_masked[which(hab.map_masked$values > as.numeric(hab.class[x,2]) & hab.map_masked$values < as.numeric(hab.class[x,3])),]
              
              test <- 
                st_union(st_intersection(EOO.shp.proj[which(EOO.shp.proj$tax == hab.class[x,1]),], 
                                         hab.map_poly))
              
            test
              
            }
          
          EOO.shp.proj.hab.map <- 
            do.call('rbind', lapply(output, function(x) st_as_sf(x)))
          EOO.shp.proj.hab.map <- cbind(EOO.shp.proj.hab.map, tax = EOO.shp.proj$tax, hab.map = i)
          
          names(EOO.shp.proj.hab.map)[names(EOO.shp.proj.hab.map) == attr(EOO.shp.proj.hab.map, "sf_column")] = "geometry"
          st_geometry(EOO.shp.proj.hab.map) = "geometry"
          
          # EOO.shp.hab.map_list <- vector('list', nrow(hab.class))
          # for (j in 1:nrow(hab.class)) {
          #   
          #   hab.map_poly <- 
          #     hab.map_masked[which(hab.map_masked$layer > as.numeric(hab.class[j,2]) & hab.map_masked$layer < as.numeric(hab.class[j,3])),]
          #   
          #   EOO.shp.hab.map_list[[j]] <- 
          #     st_union(st_intersection(EOO.shp.proj[which(EOO.shp.proj$tax == hab.class[j,1]),], hab.map_poly))
          #   
          # }
          # 
          # 
          # EOO.shp.proj.hab.map <- do.call('rbind', lapply(EOO.shp.hab.map_list, function(x) st_as_sf(x)))
          # EOO.shp.proj.hab.map <- cbind(EOO.shp.proj.hab.map, tax = EOO.shp.proj$tax)
        }
        
      if (any(class(hab.class) == "numeric")) {
        message("intersection with EOO polygons")
        
        EOO.shp.hab.map_list <- vector('list', length(hab.class))
        for (j in 1:length(hab.class)) {
          # EOO.shp.hab.map <- vector('list', nrow(EOO.shp.proj))
          # for (k in 1:nrow(EOO.shp.proj)) {
          #
          #   EOO.shp.hab.map[[k]] <-
          #     st_intersection(EOO.shp.proj[k,], hab.map_poly)
          #
          # }
          
          if (parallel) {
            cl <- snow::makeSOCKcluster(NbeCores)
            doSNOW::registerDoSNOW(cl)
            
            message('Parallel running with ',
                    NbeCores, ' cores')
            
            `%d%` <- foreach::`%dopar%`
          } else {
            `%d%` <- foreach::`%do%`
          }
          
          x <- NULL
          if (show_progress) {
            pb <-
              txtProgressBar(min = 0,
                             max = nrow(EOO.shp.proj),
                             style = 3)
            
            progress <- function(n)
              setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
          } else {
            opts <- NULL
          }
          
          output <-
            foreach::foreach(
              x = 1:nrow(EOO.shp.proj),
              .options.snow = opts,
              .packages = 'sf'
            ) %d% {
              if (!parallel & show_progress)
                setTxtProgressBar(pb, x)
              
              test <-
                st_intersection(EOO.shp.proj[x, ], hab.map_masked[[j]])
              
              test
              
            }
          
          EOO.shp.hab.map_list[[j]] <-
            do.call('rbind', output)
          
        }
        
        EOO.shp.proj.hab.map <- EOO.shp.hab.map_list
        
      }
        
      if (is.null(hab.class)) {
        hab.map_poly <-
          st_union(hab.map_masked)
        
        EOO.shp.proj.hab.map <-
          st_intersection(EOO.shp.proj, hab.map_poly)
        
      }
      
      # values(hab.map_masked_crop) <- 
      #   1:length(hab.map_masked_crop@data@values)
      
      # intersect_rast_EOO <- raster::extract(hab.map_masked_crop, EOO.shp.proj)
      
      # hab.map_masked_crop[[1]]
      
      AOH.poly[[i]] <- 
        EOO.shp.proj.hab.map
      
      if (any(class(EOO.shp.proj.hab.map) == "list")) {
        
        
         for (j in 1:length(EOO.shp.proj.hab.map)) {
           
           res. <- data.frame(species = EOO.shp.proj.hab.map[[j]]$tax,
                              AOH = as.numeric(st_area(EOO.shp.proj.hab.map[[j]])) / 1000000)
           
           colnames(res.)[2] <-
             paste0("hab.map_", hab.class[j],"_", ifelse(is.null(years), i, names(hab.map)[i]))
           
           colnames(EOO.res.all)[1] <- "species"
           
           EOO.res.all <- 
             merge(EOO.res.all, res.,
                   by = "species", 
                   all.x = T)
         }
        
      } else {
        
        res. <- data.frame(species = EOO.shp.proj.hab.map$tax,
                           AOH = as.numeric(st_area(EOO.shp.proj.hab.map)) / 1000000)
        
        colnames(res.)[2] <-
          paste0("hab.map_", ifelse(is.null(years), i, names(hab.map)[i]))
        
        colnames(EOO.res.all)[1] <- "species"
        
        EOO.res.all <- 
          merge(EOO.res.all, res.,
                by = "species", all.x = T)
        
      }
    }
    
    if (any(grepl("sf", class(hab.map.selected)))) {
      
      if (st_crs(hab.map.selected) != st_crs(EOO.shp.proj)) {
        
        hab.map.selected <- st_transform(hab.map.selected, proj_type_)
        
      }
      
      if (any(!st_is_valid(hab.map.selected))) {
        message("sf object provided not valid - trying to make it valid")
        hab.map.selected <- 
          sf::st_make_valid(hab.map.selected)
        
      }
      
      if (hab.map.type[i])
        EOO.shp.proj.hab.map <-
          st_intersection(EOO.shp.proj, st_union(hab.map.selected))
      
      if (!hab.map.type[i])
        EOO.shp.proj.hab.map <-
          st_difference(EOO.shp.proj, st_union(hab.map.selected))
      
      EOO.shp.proj.hab.map <- cbind(EOO.shp.proj.hab.map,  hab.map = i)
      
      names(EOO.shp.proj.hab.map)[names(EOO.shp.proj.hab.map) == attr(EOO.shp.proj.hab.map, "sf_column")] = "geometry"
      st_geometry(EOO.shp.proj.hab.map) = "geometry"
      
      res. <- data.frame(species = EOO.shp.proj.hab.map$tax,
                         AOH = round(as.numeric(st_area(EOO.shp.proj.hab.map)) / 1000000, 3))
      
      colnames(res.)[2] <-
        paste0("hab.map.", i)
      
      EOO.res.all <- 
        merge(EOO.res.all, res.,
              by = "species", all.x = T)
      
      EOO.res.all[, i + 2][which(!is.na(EOO.res.all[,2]) & is.na(EOO.res.all[, i + 2]))] <- 0
      
      AOH.poly[[i]] <- 
        EOO.shp.proj.hab.map
      
      # hab.map.selected <-
      #   sf::st_intersection(hab.map.selected, EOO.shp.proj)
      
    }
  }
  
  # all hab.map together
  
  if (any(classes.hab.map$rast) & is.data.frame(hab.class)) {
    
    all_poly <- bind_rows(AOH.poly)
    
    nbe_sp <- length(unique(all_poly$tax))
    
    all_sp <- unique(all_poly$tax)
    
    if (parallel) {
      
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      
      message('Parallel running with ',
              NbeCores, ' cores')
      
      `%d%` <- foreach::`%dopar%`
      
    } else {
      
      `%d%` <- foreach::`%do%`
      
    }
    
    x <- NULL
    if (show_progress) {
      pb <-
        txtProgressBar(min = 0,
                       max = nbe_sp,
                       style = 3)
      
      progress <- function(n)
        setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }
    
    output <-
      foreach::foreach(
        x = 1:nbe_sp,
        .options.snow = opts,
        .packages = 'sf'
      ) %d% {
        
        if (!parallel & show_progress)
          setTxtProgressBar(pb, x)
        
        selected_poly <- all_poly[which(all_poly$tax == all_sp[x]), ]
        
        intersect_poly <- st_intersection(selected_poly)
        intersect_poly <- intersect_poly[which(intersect_poly$n.overlaps > 1), ]

        intersect_poly[,c("tax", "geometry")]
        
      }
    
    EOO.shp.proj.hab.map <- 
      do.call('rbind', lapply(output, function(x) st_as_sf(x)))
    
    res. <- data.frame(species = EOO.shp.proj.hab.map$tax,
                       AOH = round(as.numeric(st_area(EOO.shp.proj.hab.map)) / 1000000, 3))
    
    colnames(res.)[2] <-
      "all.hab.map"
    
    EOO.res.all <- 
      merge(EOO.res.all, res.,
            by = "species", all.x = T)
    
    EOO.res.all[, ncol(EOO.res.all)][which(!is.na(EOO.res.all[,2]) & is.na(EOO.res.all[, ncol(EOO.res.all)]))] <- 0
    
    AOH.poly[[length(AOH.poly)]] <- 
      EOO.shp.proj.hab.map
    
    
  }
  
  return(list(EOO = EOO.shp.proj,
              AOH = EOO.res.all,
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









