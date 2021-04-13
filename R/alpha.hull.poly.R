#' Internal function
#'
#' Alpha hull process
#'
#' @param XY data.frame coordinates
#' @param alpha integer if mode planar, in kilometers, if mode spheroid in decimal degrees
#' @param buff numeric
#' @param exclude.area logical
#' @param poly_exclude
#'
#' @details 
#' The functions ahull_to_SPLDF and alpha.hull.poly were originally posted in the website https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
#' in a now broken link. It is also used in functions written by David Bucklin, see https://github.com/dnbucklin/r_movement_homerange 
#'
#' @import sf
#' @importFrom rgeos gBuffer
#' @importFrom sp SpatialPolygons Polygons CRS proj4string
#' @importFrom stats median quantile dist
#' @importFrom methods slot
#' @importFrom utils installed.packages
#' 
alpha.hull.poly <-
  function(XY,
           alpha = 1,
           buff = 0.1,
           exclude.area = FALSE,
           poly_exclude = NULL,
           mode = "spheroid",
           proj_type = "cea")
  {
    
    
    if (mode == "planar") {
      
      if(class(proj_type) != "CRS") {
        
        projEAC <- proj_crs(proj_type = proj_type)
        
      } else {
        
        projEAC <- proj_type
        
      }
      
      XY <-
        sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(projEAC),
          pts = XY[, c(1, 2)]
        )
      
      if (alpha < median(dist(XY,  upper = F))/1000/10) {
        
        warning("alpha value is defined in kilometer given planar estimation and might be too small to get any polygon")
        
      }

      if (buff < quantile(dist(XY,  upper = F), probs = 0.01)/1000/5) {
        warning("buff value is defined in kilometer given planar estimation and might be too small to get any polygon")
      }
      
            
      buff <-
        buff*1000

      alpha <-
        alpha*1000

    }

    
    Used_data <- unique(XY)
    if (any(rownames(installed.packages()) == "alphahull")) {
      
      loadNamespace("alphahull")
      
      run_alpha <- TRUE
      while(run_alpha) {
        
        ahull.obj <-
          try(alphahull::ahull(Used_data[, c(1, 2)], alpha = alpha), silent = T)
        
        if(class(ahull.obj) != "try-error") {
          
          run_alpha <- FALSE
          
        } else {
          
          Used_data <- apply(Used_data, 2, jitter)
          
        }
      }
      
      y.as.spldf <- ahull_to_SPLDF(ahull.obj)
      y.as.spldf_buff <- rgeos::gBuffer(y.as.spldf, width = buff)
      
      NZp <- slot(y.as.spldf_buff, "polygons")
      holes <-
        lapply(NZp, function(x)
          sapply(slot(x, "Polygons"), slot,
                 "hole"))
      res <- lapply(1:length(NZp), function(i)
        slot(NZp[[i]],
                      "Polygons")[!holes[[i]]])
      IDs <- row.names(y.as.spldf_buff)
      NZfill <- sp::SpatialPolygons(
        lapply(1:length(res), function(i)
          sp::Polygons(res[[i]], ID = IDs[i])),
        proj4string = sp::CRS(sp::proj4string(y.as.spldf_buff), doCheckCRSArgs = TRUE)
      )
      
      # crs <- sp::CRS("+proj=longlat +datum=WGS84", doCheckCRSArgs = TRUE)
      # raster::crs(NZfill) <- crs
    
      NZfill <-
        suppressWarnings(rgeos::gBuffer(NZfill, byid = TRUE, width = 0))
      
      if (mode == "planar") {
        
        NZfill <- as(NZfill, "sf")
        
        sf::st_crs(NZfill) <-
          projEAC
        
        # sp::proj4string(NZfill) <- projEAC

        if (exclude.area) {

          poly_exclude_proj <-
            st_transform(st_make_valid(poly_exclude), crs = projEAC)


          ### work around to avoid bug appearing randomly
          # cont_try <- TRUE
          # while(cont_try) {
          #   NZfill <-
          #     try(suppressWarnings(suppressMessages(st_union(
          #       st_union(st_intersection(st_make_valid(NZfill), poly_exclude_proj))
          #     ))))
          #
          #   if (class(NZfill) != "try-error")
          #     cont_try <- FALSE
          # }


          # NZfill <-
          #   st_union(st_intersection(st_make_valid(NZfill), poly_exclude_proj))

          NZfill <-
            suppressWarnings(suppressMessages(st_union(
              st_intersection(st_make_valid(NZfill), poly_exclude_proj)
            )))
          
          
          sf::st_crs(NZfill) <-
            projEAC

          if(length(NZfill) == 0) {
            
            warning("After excluding areas, the alpha hull is empty. EOO is NA.")

            NZfill <- NA
          }

          # NZfill <- as(NZfill, "Spatial")

        }
      }

      if (mode == "spheroid") {
        
        NZfill <- as(NZfill, "sf")
        
        sf::st_crs(NZfill) <-
          4326
        
        # p1_lines <- suppressWarnings(sf::st_cast(NZfill, "MULTILINESTRING"))
        # p1_lines_seg <-
        #   sf::st_segmentize(p1_lines, units::set_units(20, km))
        # p1 <- sf::st_cast(p1_lines_seg, "POLYGON")
        
        # sp::proj4string(NZfill) <- sp::CRS(SRS_string = 'EPSG:4326')
        
        # NZfill_vertices <-
        #   try(suppressWarnings(geosphere::makePoly(NZfill)), silent = TRUE)
        # 
        # if(class(NZfill_vertices) == "try-error") {
        #   NZfill_vertices <- 
        #     NZfill
        #   
        #   warning("Failed to add vertices to polygon, be careful with EOO estimation if large EOO")
        # }
        
        if (exclude.area) {

          # NZfill_sf <- as(NZfill, "sf")

          # poly_exclude <-
          #   st_make_valid(poly_exclude)

          NZfill <-
            suppressWarnings(suppressMessages(st_union(
            st_intersection(st_make_valid(NZfill), poly_exclude)
          )))

          # NZfill <- NZfill_sf

          if(length(NZfill) == 0) {
            
            warning("After excluding areas, the alpha hull is empty. EOO is set to NA.")

            NZfill <- NA
            
          } 
          # else {
          # 
          #   sf::st_crs(NZfill) <-
          #     4326
          # 
          #   # NZfill <- as(NZfill, "Spatial")
          # 
          # }

        }
        # else {
        #   
        #   NZfill <- as(NZfill, "sf")
        #   
        # }
        
      }
      
      # raster::crs(NZfill) <- "+proj=longlat +datum=WGS84"
      
      NZfill <- st_sf(geom = NZfill)
      
    } else{
      
      stop("The package alphahull is required for this procedure, please install it")
      
    }
    
    return(NZfill)
}
