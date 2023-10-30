#' @title Internal function
#'
#' @description Build convex hull polygon
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#' 
#' @param XY data.frame
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param proj_type crs
#' 
#' @import sf
#' @importFrom grDevices chull
#' @keywords internal
Convex.Hull.Poly <-
  function(XY,
           mode = "spheroid",
           proj_type = NULL) {
    

    # if(any(grepl('Spatial', class(poly_exclude))))
    #   poly_exclude <- as(poly_exclude, "sf")
    
    if (mode == "spheroid") {
      
      hpts <- grDevices::chull(x =  XY[, 1], y = XY[, 2])
      hpts <- c(hpts, hpts[1])
      
      coord <- matrix(NA, length(hpts), 2)
      
      coord[,1] <- XY[hpts, 1]
      coord[,2] <- XY[hpts, 2]
      
      
      POLY <- st_polygon(x = list(coord))

      p1 <- sf::st_sf(tax = "tax", geometry = list(POLY))
      
      
      # p1 <- sf::st_as_sf(data.frame(a = 1, geometry = POLY), wkt = "geometry")
      sf::st_crs(p1) <-
        4326
      p1_lines <- suppressWarnings(sf::st_cast(p1, "LINESTRING"))
      p1_lines_seg <-
        sf::st_segmentize(p1_lines, units::set_units(20, km))
      p1 <- sf::st_cast(p1_lines_seg, "POLYGON")
      p1 <- sf::st_make_valid(p1)
      

      
    }
    
    if (mode == "planar") {
      
      p1 <- sf::st_convex_hull(x = sf::st_multipoint(as.matrix(XY[, c(1, 2)])))
      # eoo <- st_area(p1)
      
      p1 <-
        sf::st_sfc(p1)
      
      sf::st_crs(p1) <- proj_type
      
      
      p1 <- st_sf(geom = p1)
      
      p1 <- p1[sf::st_is(p1, c("MULTIPOLYGON", "POLYGON")),]
      
      p1 <- sf::st_make_valid(p1)
      
      # p1 <- suppressWarnings(geosphere::makePoly(as(p1, "Spatial")))
      
    }
    
    return(p1)
  }