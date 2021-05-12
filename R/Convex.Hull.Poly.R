
#' @title Internal function
#'
#' @description Build convex hull polygon
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#' 
#' @param XY data.frame
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param exclude.area logical
#' @param poly_exclude polygon
#' @param proj_type character string or numeric or object of CRS class, by default is "cea"
#' 
#' @import sf
#' @importFrom grDevices chull
#' @importFrom rgeos readWKT
#' 
Convex.Hull.Poly <-
  function(XY,
           mode = "spheroid",
           exclude.area = FALSE,
           poly_exclude = NULL,
           proj_type = "cea") {
    
    if (exclude.area & is.null(poly_exclude))
      stop("exclude.area is TRUE but no shape provided")
    
    if(any(grepl('Spatial', class(poly_exclude))))
      poly_exclude <- as(poly_exclude, "sf")
    
    if (mode == "spheroid") {
      
      hpts <- grDevices::chull(x =  XY[, 1], y = XY[, 2])
      hpts <- c(hpts, hpts[1])
      
      coord <- matrix(NA, length(hpts), 2)
      
      coord[,1] <- XY[hpts, 1]
      coord[,2] <- XY[hpts, 2]
      
      # POLY <- "POLYGON(("
      # for (i in 1:length(hpts)) {
      #   # POLY <- paste(POLY, XY[hpts[i], 1], " ", XY[hpts[i], 2], sep = "")
      #   # if (i != length(hpts))
      #   #   POLY <- paste(POLY, ", ", sep = "")
      #   # if (i == length(hpts))
      #   #   POLY <- paste(POLY, "))", sep = "")
      # 
      #   coord[i, 1] <- XY[hpts[i], 2]
      #   coord[i, 2] <- XY[hpts[i], 1]
      # 
      # }
      
      POLY <- st_polygon(x = list(coord))

      p1 <- sf::st_sf(a = 1, geometry = list(POLY))
      
      # p1 <- sf::st_as_sf(data.frame(a = 1, geometry = POLY), wkt = "geometry")
      sf::st_crs(p1) <-
        4326
      p1_lines <- suppressWarnings(sf::st_cast(p1, "LINESTRING"))
      p1_lines_seg <-
        sf::st_segmentize(p1_lines, units::set_units(20, km))
      p1 <- sf::st_cast(p1_lines_seg, "POLYGON")
      p1 <- sf::st_make_valid(p1)
      
      # p1_sp <- rgeos::readWKT(POLY) 
      # p1_sp <- suppressWarnings(geosphere::makePoly(p1_sp)) 
      
      ## Gilles, there is a warning here that may be a potential the problem... use sf to make poly?
      # sp::proj4string(p1) <-
      #   sp::CRS("+proj=longlat +datum=WGS84", doCheckCRSArgs = TRUE)
      
      # p1 <- suppressWarnings(geosphere::makePoly(p1))
      
      ##Gilles, why not define directly a polygon using 'sf'? See code examples commented below:
      #p1_sf <- sf::st_sfc(sf::st_polygon(list(coord[,2:1])))
      #sf::st_crs(p1_sf) <- sf::st_crs(poly_exclude)
      
      if (exclude.area) {
        # p1_sf <- as(p1, "sf")
        
        p1 <-
          suppressWarnings(suppressMessages(sf::st_union(
            sf::st_intersection(p1, poly_exclude)
          )))
        
        sf::st_crs(p1) <-
          4326
        
        if(length(p1) == 0) {
          warning("After excluding areas, the convex hull is empty. EOO is set as NA.")
         
          p1 <- NA 
        }
        
        
        
        p1 <- st_sf(geom = p1)
        
        # else {
        #   
        #   p1 <- 
        #     as(p1, "Spatial")
        #   
        # }
        
      }
      
    }
    
    if (mode == "planar") {
      
      if(class(proj_type) != "CRS") {
        
        projEAC <- proj_crs(proj_type = proj_type)
        
      } else {
        
        projEAC <- proj_type
        
      }
      
      XY_sf_proj <-
        sf::sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(projEAC),
          pts = XY[, c(1, 2)]
        )
      
      p1 <- sf::st_convex_hull(x = sf::st_multipoint(XY_sf_proj))
      # eoo <- st_area(p1)
      
      p1 <-
        sf::st_sfc(p1)
      
      sf::st_crs(p1) <- projEAC
      
      if (exclude.area) {

        poly_exclude_proj <-
          sf::st_transform(poly_exclude, crs = projEAC)
        
        p1 <-
          sf::st_union(sf::st_intersection(p1, poly_exclude_proj))
        
        if(length(p1) == 0) {
          warning("After excluding areas, the convex hull is empty. EOO is NA.")
          
          p1 <- NA 
        } 
        # else {
        #   
        #   p1 <- 
        #     as(p1, "Spatial")
        #   
        # }
        
      }
      
      p1 <- st_sf(geom = p1)
      
      # p1 <- suppressWarnings(geosphere::makePoly(as(p1, "Spatial")))
      
    }
    
    return(p1)
  }