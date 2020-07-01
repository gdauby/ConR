
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
#' @importFrom grDevices chull
#' @importFrom rgeos readWKT
#' @importFrom geosphere makePoly
#' 
Convex.Hull.Poly <-
  function(XY,
           mode = "spheroid",
           exclude.area = FALSE,
           poly_exclude = NULL,
           proj_type = "cea") {
    
    if (exclude.area & is.null(poly_exclude))
      stop("exclude.area is TRUE but no shape provided")
    
    
    if (mode == "spheroid") {
      hpts <- grDevices::chull(x =  XY[, 1], y = XY[, 2])
      hpts <- c(hpts, hpts[1])
      coord <- matrix(NA, length(hpts), 2)
      POLY <- "POLYGON(("
      for (i in 1:length(hpts)) {
        POLY <- paste(POLY, XY[hpts[i], 1], " ", XY[hpts[i], 2], sep = "")
        if (i != length(hpts))
          POLY <- paste(POLY, ", ", sep = "")
        if (i == length(hpts))
          POLY <- paste(POLY, "))", sep = "")
        
        coord[i, 1] <- XY[hpts[i], 2]
        coord[i, 2] <- XY[hpts[i], 1]
        
      }
      
      p1 <- rgeos::readWKT(POLY)
      
      sp::proj4string(p1) <-
        sp::CRS("+proj=longlat +datum=WGS84", doCheckCRSArgs = TRUE)
      
      p1 <- suppressWarnings(geosphere::makePoly(p1))
      
      
      if (exclude.area) {
        p1_sf <- as(p1, "sf")
        
        p1 <-
          suppressWarnings(suppressMessages(st_union(
            st_intersection(p1_sf, poly_exclude)
          )))
        
        
        st_crs(p1) <-
          "+proj=longlat +datum=WGS84"
        
        p1 <- as(p1, "Spatial")
        
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
      
      p1 <- st_convex_hull(x = st_multipoint(XY_sf_proj))
      # eoo <- st_area(p1)
      
      p1 <-
        st_sfc(p1)
      
      st_crs(p1) <- projEAC
      
      if (exclude.area) {
        poly_exclude_proj <-
          st_transform(poly_exclude, crs = projEAC)
        
        p1 <-
          st_union(st_intersection(p1, poly_exclude_proj))
        
      }
      
      
      # p1 <- suppressWarnings(geosphere::makePoly(as(p1, "Spatial")))
      
    }
    
    return(p1)
  }