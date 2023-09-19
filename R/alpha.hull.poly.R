#' Internal function
#'
#' Alpha hull process
#'
#' @param XY data.frame coordinates
#' @param alpha integer if mode planar, in kilometers, if mode spheroid in decimal degrees
#' @param buff numeric
#' @param exclude.area logical
#' @param poly_exclude sf polygon
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param proj_type character string or numeric or object of crs class, by default is "cea"
#' 
#' @details 
#' The functions ahull_to_SPLDF and alpha.hull.poly were originally posted in the website https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
#' in a now broken link. It is also used in functions written by David Bucklin, see https://github.com/dnbucklin/r_movement_homerange 
#'
#' @import sf
#' @importFrom stats median quantile dist
#' @importFrom methods slot
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
    
    if (!try(requireNamespace("alphahull", quietly = F), silent = TRUE))
      stop("The package alphahull is required for this procedure, please install it")
    
    # if (!try(requireNamespace("rgeos", quietly = F), silent = TRUE))
    #   stop("The package rgeos is required for this procedure, please install it")
    
    requireNamespace("alphahull")
    # requireNamespace("rgeos")

    if (mode == "planar") {
      
      projEAC <- proj_crs(proj_type = proj_type)
      
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
    
    Used_data <- XY
    
    run_alpha <- TRUE
    while (run_alpha) {
      ahull.obj <-
        try(alphahull::ahull(Used_data[, c(1, 2)], alpha = alpha), silent = T)
      
      if (class(ahull.obj) != "try-error") {
        run_alpha <- FALSE
        
      } else {
        Used_data <- apply(Used_data, 2, jitter)
        
      }
    }
    
    y.as.spldf <- ahull_to_SPLDF(ahull.obj)
    y.as.spldf_sf <- st_as_sf(y.as.spldf)
    NZfill <- st_buffer(y.as.spldf_sf, buff)
    if (mode == "planar")
      st_crs(NZfill) <- projEAC
    if (mode == "spheroid")
      st_crs(NZfill) <- 4326
    
      
      # 
      # # y.as.spldf_buff <- rgeos::gBuffer(y.as.spldf, width = buff)
      # 
      # NZp <- slot(y.as.spldf_buff, "polygons")
      # holes <-
      #   lapply(NZp, function(x)
      #     sapply(slot(x, "Polygons"), slot,
      #            "hole"))
      # res <- lapply(1:length(NZp), function(i)
      #   slot(NZp[[i]],
      #                 "Polygons")[!holes[[i]]])
      # IDs <- row.names(y.as.spldf_buff)
      # NZfill <- sp::SpatialPolygons(
      #   lapply(1:length(res), function(i)
      #     sp::Polygons(res[[i]], ID = IDs[i])),
      #   proj4string = sp::CRS(sp::proj4string(y.as.spldf_buff), doCheckCRSArgs = TRUE)
      # )
      
      # crs <- sp::CRS("+proj=longlat +datum=WGS84", doCheckCRSArgs = TRUE)
      # raster::crs(NZfill) <- crs
    
      # NZfill <-
      #   suppressWarnings(rgeos::gBuffer(NZfill, byid = TRUE, width = 0))
      
      if (exclude.area) {
        
        if (mode == "planar") {
          poly_exclude_proj <-
            st_transform(st_make_valid(poly_exclude), crs = projEAC)
          
          NZfill <-
            suppressWarnings(suppressMessages(st_union(
              st_intersection(st_make_valid(NZfill), poly_exclude_proj)
            )))
          
          if(length(NZfill) == 0) {
            
            warning("After excluding areas, the alpha hull is empty. EOO is NA.")
            
            NZfill <- NA
          }
          
        }
        
        if (mode == "spheroid") {
          
          
          NZfill <-
            suppressWarnings(suppressMessages(st_union(
              st_intersection(st_make_valid(NZfill), poly_exclude)
            )))
          
          if(length(NZfill) == 0) {
            
            warning("After excluding areas, the alpha hull is empty. EOO is set to NA.")
            
            NZfill <- NA
            
          } 
          
        }
        
        
      }
      
      # NZfill <- st_sf(geom = NZfill)

    return(NZfill)
}
