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
#' @import sp
#' @importFrom rgeos gBuffer
#' @importFrom geosphere makePoly
#' @importFrom sf st_union st_intersection
#' 
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
        sf::sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(projEAC),
          pts = XY[, c(1, 2)]
        )
      
      if (alpha < median(dist(XY,  upper = F))/1000/10) {
        warning("alpha value is defined in kilometer given planar estimation and might be too small to get any polygon")
      }

      if (buff < quantile(dist(XY,  upper = F), probs = 0.01)/1000/5) {
        warning("alpha value is defined in kilometer given planar estimation and might be too small to get any polygon")
      }
      
            
      buff <-
        buff*1000

      alpha <-
        alpha*1000

    }

    
    Used_data = unique(XY)
    Used_data <- apply(Used_data, 2, jitter)
    
    if (any(rownames(utils::installed.packages()) == "alphahull")) {
      
      loadNamespace("alphahull")
      ahull.obj <-
        alphahull::ahull(Used_data[, c(1, 2)], alpha = alpha)
      y.as.spldf <- ahull_to_SPLDF(ahull.obj)
      y.as.spldf_buff <- rgeos::gBuffer(y.as.spldf, width = buff)
      
      NZp <- methods::slot(y.as.spldf_buff, "polygons")
      holes <-
        lapply(NZp, function(x)
          sapply(methods::slot(x, "Polygons"), slot,
                 "hole"))
      res <- lapply(1:length(NZp), function(i)
        methods::slot(NZp[[i]],
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
            sf::st_transform(sf::st_make_valid(poly_exclude), crs = projEAC)
          
          NZfill <-
            sf::st_union(sf::st_intersection(sf::st_make_valid(NZfill), poly_exclude_proj))
          
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
        
        sp::proj4string(NZfill) <- sp::CRS(SRS_string = 'EPSG:4326')
        
        NZfill <- suppressWarnings(geosphere::makePoly(NZfill))
        
        if (exclude.area) {
          NZfill_sf <- as(NZfill, "sf")
          
          poly_exclude <-
            sf::st_make_valid(poly_exclude)
          
          NZfill <-
            suppressWarnings(suppressMessages(sf::st_union(
              sf::st_intersection(sf::st_make_valid(NZfill_sf), poly_exclude)
            )))
          
          if(length(NZfill) == 0) {
            warning("After excluding areas, the alpha hull is empty. EOO is NA.")
            
            NZfill <- NA 
          } else {
            
            sf::st_crs(NZfill) <-
              "+proj=longlat +datum=WGS84"
            
            NZfill <- as(NZfill, "Spatial")
            
          }
          
        }
        
      }
      
      # raster::crs(NZfill) <- "+proj=longlat +datum=WGS84"
      
    } else{
      stop("The package alphahull is required for this procedure, please install it")
    }
    
    return(NZfill)
}
