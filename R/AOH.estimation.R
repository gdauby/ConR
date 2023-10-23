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
#' @param proj_type character string, numeric or object of CRS class. Default to
#'   "cea".
#' @param mode character string either 'spheroid' or 'planar'. By default
#'   'spheroid'
#' @param hab.map raster, raster layer/stack or spatial polygons containing the
#'   habitat spatial information
#' @param hab.map.type logical, vector of same lenght of hab.map, if TRUE means
#'   hab.map is suitable for species, if FALSE unsuitable
#' @param hab.class classes of values in ```hab.map``` to be considered as
#'   'habitat'.
#' @param years numeric. Time interval between the first and last ```hab.map```
#'   if more than one raster is provided (e.g. if ```hab.map``` is a RasterStack
#'   object).
#' @param country_map a `SpatialPolygonsDataFrame` or
#' @param exclude.area a logical, if TRUE, areas outside of `country_map` are
#'   cropped of `SpatialPolygons` used for calculating EOO. By default is FALSE
#' @param parallel logical. Should computing run in parallel? Default to FALSE.
#' @param NbeCores integer. The number of cores for parallel execution. Default
#'   to 2.
#' @param simplifiy_poly logical whether the resulting polygon should be
#'   simplified using ms_simplify function of rmapshaper package
#' @param buffer logical
#' @param ... additional arguments passed to `EOO.computing`
#' 
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
#' @rawNamespace import(terra, except = c(points, median, na.omit, quantile, head, tail, predict))
#' @importFrom stars st_as_stars
#' @importFrom rmapshaper ms_simplify
#' @importFrom dplyr bind_rows
#' @import sf
#' 
#' @export AOH.estimation
#' 
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
  
  proj_type_ <- proj_crs(proj_type = proj_type)
  
  EOO.res <- EOO.computing(XY = XY, export_shp = TRUE, mode = mode, exclude.area = FALSE, ...) # 
  
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
      
      hab.map.selected <-
        terra::crop(hab.map.selected, terra::vect(EOO.shp.proj))
      # raster::crop(hab.map.selected, st_bbox(EOO.shp.proj)) # cropping raster to the extent of the polygon
      
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
            
            if (!requireNamespace("rmapshaper", quietly = TRUE))
              stop(
                "The 'rmapshaper' package is required to run this function.",
                "Please install it first."
              )
            
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
          
          activate_parallel(parallel = parallel)
          
          pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_data))
          opts <- pro_res$opts
          pb <- pro_res$pb
          
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
          
          if(parallel) parallel::stopCluster(cl)
          if(show_progress) close(pb)
          
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
          
          activate_parallel(parallel = parallel)
          
          pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_data))
          opts <- pro_res$opts
          pb <- pro_res$pb
          
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
          
          if(parallel) parallel::stopCluster(cl)
          if(show_progress) close(pb)
          
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
    
    activate_parallel(parallel = parallel)
    
    pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_data))
    opts <- pro_res$opts
    pb <- pro_res$pb
    
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
    
    if(parallel) parallel::stopCluster(cl)
    if(show_progress) close(pb)
    
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






