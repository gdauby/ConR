#' @title Compute Extent of Occurrence
#'
#' @description Compute EOO given a data frame of coordinate in decimals degrees
#'
#' @param XY data.frame
#' @param exclude.area logical. Default if FALSE
#' @param country_map SpatialPolygonDataframe. Default if NULL
#' @param Name_Sp string
#' @param method.range string, by default "convex.hull", can also take "alpha.hull"
#' @param alpha integer
#' @param buff.alpha numeric
#' @param method.less.than3 string
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param proj_type character string or numeric or object of CRS class, by default is "cea"
#'
#' @author Dauby, G. & Lima, R.A.F.
#' 
#' @return A list with one element being a numeric vector [[1]]EOO and the polygon as Simple feature collection [[2]]spatial.polygon
#'  
#' @details If `exclude.area` is TRUE and country_map is not provided, 
#' the world country polygons used comes from the package \href{https://www.rdocumentation.org/packages/rnaturalearth/versions/0.1.0/topics/ne_countries}{rnaturalearth}
#' 
#' By default (`mode = "spheroid"`),the area of the polygon is based 
#' on a polygon in longitude/latitude coordinates considering the earth as an ellipsoid.  
#' 
#' To make a polygon more accurate, the function use the function `st_segmentize` from the \href{https://CRAN.R-project.org/package=sf}{sf} package.
#' This adds vertices on the great circles (in order to make shortest distances between points, see example below) 
#' which can make difference for species with large distribution.  
#' 
#' An estimation of EOO based on projected data is also possible (`mode = "planar"`).
#' This allow the user to use its own projection.
#' By default, the projection (`proj_type = "cea"`) is a \href{https://epsg.io/6933}{global cylindrical equal-area projection}.
#' It can makes sense to use a planar projection to estimate the EOO. See example.
#' 
#' It is possible to use another projection by providing the EPSG code to `proj_type`. See example.
#' Check \href{https://www.rdocumentation.org/packages/rgdal/versions/0.4-10/topics/make_EPSG}{make_EPSG}
#' from the `rgdal` packag for getting the EPSG code.
#' 
#' For the very specific (and infrequent) case where all occurrences are
#' localized on a straight line (in which case EOO would be null), 
#' 'noises' are added to coordinates, using the `jitter` function, 
#' see \href{https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/jitter}. 
#' There is a warning when this happens. This means that EOO value will not be constant 
#' across multiple estimation (although the variation should be small)
#' 
#' @examples
#' 
#' country_map <-
#'   rnaturalearth::ne_countries(scale = 50, returnclass = "sp")
#'
#' ### example Large distribution
#' # Spheroid estimation
#' XY <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
#' p1 <- 
#'   EOO.comp(XY =  XY[,c(2, 1)])
#' plot(country_map)
#' plot(p1$spatial.polygon, lwd = 2, col = 'red', add = TRUE)
#' 
#' p2 <- 
#'   EOO.comp(XY = XY[,c(2, 1)], mode = "planar")
#' plot(st_as_sf(country_map), reset = FALSE, extent = as(p2$spatial.polygon, "sf"))
#' plot(p2$spatial.polygon, lwd = 2, col = 'red', add = TRUE)
#' 
#' # World Mercartor projection
#' p2 <- 
#'   EOO.comp(XY = XY[,c(2, 1)], mode = "planar", proj_type = 3395)
#' plot(st_as_sf(country_map), reset = FALSE, extent = as(p2$spatial.polygon, "sf"))
#' plot(p2$spatial.polygon, lwd = 2, col = 'red', add = TRUE)
#' 
#' ### example Antartic distribution
#' XY <- rbind(c(-150, -65), c(0, -62), c(120, -78), c(150, -65))
#' 
#' p1 <- 
#'   EOO.comp(XY = XY[,c(2, 1)], mode = "planar", proj_type = "Antarctic")
#' 
#' p1 <- st_transform(as(p1$spatial.polygon, "sf"), proj_crs(proj_type = "Antarctic"))
#' plot(st_transform(st_as_sf(country_map), proj_crs(proj_type = "Antarctic")), extent = p1)
#' plot(p1, lwd = 2, col = 'red', add = TRUE)
#' 
#' p1 <- 
#'   EOO.comp(XY = XY[,c(2, 1)])
#' 
#' p1 <- st_transform(as(p1$spatial.polygon, "sf"), proj_crs(proj_type = "Antarctic"))
#' plot(st_transform(st_as_sf(country_map), proj_crs(proj_type = "Antarctic")), extent = p1)
#' plot(p1, lwd = 2, col = 'red', add = TRUE)
#' 
#' 
#' @import sf
#' @importFrom sp proj4string CRS
#' @importFrom stats dist cor
#' @importFrom rgeos gBuffer
#' @importFrom rnaturalearth ne_countries
#' 
#' @export
EOO.comp <-  function(XY,
                      exclude.area = FALSE,
                      country_map = NULL,
                      Name_Sp = "tax",
                      method.range = "convex.hull",
                      alpha = 1,
                      buff.alpha = 0.1,
                      method.less.than3 = "not comp",
                      mode = "spheroid",
                      proj_type = "cea"
) {
  # , verbose=TRUE
  
  XY <- 
    coord.check(XY = XY, listing = FALSE)
  
  ### Getting by default land map if poly_borders is not provided
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  }else{
    
    if(any(grepl('sf', class(country_map))))
      country_map <- 
        suppressWarnings(as(country_map, "Spatial"))
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
    country_map <- 
      as(country_map, "sf")
  }
  
  ### Checking if the method of calculating EOO has been chosen
  # if (!convex.hull & !alpha.hull)
  #   stop("alpha.hull and convex.hull are both FALSE, one should TRUE")
  
  if (nrow(unique(XY)) > 1)
    if (max(dist(XY[, 2]), na.rm = T) >= 180)
      warning(
        paste(
          "Occurrences spans more than 180 degrees longitude for species",
          as.character(Name_Sp),
          ". EOO unlikely reliable, check the projection for a proper estimation and used a 'planar' (projected) mode"
        )
      )
  
  # if(alpha.hull) {
  #   convex.hull <- FALSE
  # }
  
  ## Check if there are less than 3 unique occurrences
  if (nrow(XY) < 3) {
    ## if there is only one occurrence, EOO is NA
    if (nrow(XY) < 2) {
      
      EOO <- NA
      message(
        paste(
          "\nEOO parameter cannot be estimated for",
          as.character(Name_Sp),
          "because there is only 1 unique occurrence"
        )
      )
      
    } else {
      if (method.less.than3 == "arbitrary") {
        
        projEAC <- proj_crs(proj_type = proj_type)
        
        # coordEAC <-
        #   as.data.frame(matrix(unlist(
        #     rgdal::project(
        #       as.matrix(unique(XY)[, 1:2]),
        #       proj = as.character(projEAC),
        #       inv = FALSE
        #     )
        #   ), ncol = 2))
        
        coordEAC <-
          sf_project(
          from = sf::st_crs(4326),
          to =
            sf::st_crs(projEAC),
          pts = XY[, c(1, 2)]
        )
        
        EOO <-
          as.numeric(dist(coordEAC / 1000) * 0.1 * dist(coordEAC / 1000))
      }
      
      if (method.less.than3 == "not comp") {
        ## if there are two unique occurences, EOO is not computed neither
        message(
          paste(
            "\nEOO parameter cannot be estimated for",
            as.character(Name_Sp),
            "because there is less than 3 unique occurrences"
          )
        )
        EOO <- NA
      }
    }
    
    OUTPUT <- list(EOO = EOO, spatial.polygon = NA)
    
  } else {
    
    ### Checking if all occurrences are on a straight line
    if (length(unique(XY[, 1])) == 1 ||
        length(unique(XY[, 2])) == 1 ||
        round(abs(cor(XY[, 1], XY[, 2])), 6) == 1) {
      
      message(
        paste(
          "\nOccurrences of",
          as.character(Name_Sp),
          "follow a straight line, 'noise' to coordinates is added"
        )
      )
      
      check_line <- TRUE
      while(check_line) {
        
        XY[,c(1, 2)] <- 
          apply(XY[,c(1, 2)], 2, function(x) jitter(x, factor = 0.1))
        
        if (round(abs(cor(XY[, 1], XY[, 2])), 6) != 1)
          check_line <- FALSE
        
      }
      
      # hpts <- unique(XY[, c(2, 1)])
      # POLY <- "LINESTRING("
      # for (Z in 1:dim(hpts)[1]) {
      #   POLY <- paste(POLY, hpts[Z, 1], " ", hpts[Z, 2], sep = "")
      #   if (Z != dim(hpts)[1])
      #     POLY <- paste(POLY, ", ", sep = "")
      #   if (Z == dim(hpts)[1])
      #     POLY <- paste(POLY, ")", sep = "")
      # }
      # p1 <- rgeos::readWKT(POLY)
      # 
      # p1 <- rgeos::readWKT(POLY)
      # sp::proj4string(p1) <- CRS(SRS_string='EPSG:4326')
      # 
      # # crs <- CRS("+proj=longlat +datum=WGS84")
      # # crs(p1) <- crs
      # 
      # p1 <-
      #   suppressWarnings(geosphere::makeLine(p1)) ### Add vertices to line
      # 
      # p1 <-
      #   suppressWarnings(rgeos::gBuffer(p1, width = buff_width)) ### Add buffer to line
      # 
      # ## If exclude.area is TRUE
      # if (exclude.area) {
      #   
      #   p1_sf <- as(p1, "sf")
      # 
      #   p1 <-
      #     suppressWarnings(suppressMessages(st_union(
      #       st_intersection(p1_sf, country_map)
      #     )))
      # 
      #   sf::st_crs(p1) <-
      #     "+proj=longlat +datum=WGS84"
      # 
      #   p1 <- as(p1, "Spatial")
      # 
      #   EOO <- 
      #     suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
      # 
      # } else {
      #   
      #   EOO <- suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
      #   
      # }
      
    }
    # else {
      
      if (method.range == "alpha.hull") {
        
        ### work around to avoid bug appearing randomly
        # cont_try <- TRUE
        # while(cont_try) {
          
          p1 <-
            alpha.hull.poly(
              XY = XY[, c(2, 1)],
              alpha = alpha,
              buff = buff.alpha,
              exclude.area = exclude.area,
              poly_exclude = country_map,
              mode = mode,
              proj_type = proj_type)

        #   if (!grepl("trye-error", class(p1)))
        #     cont_try <- FALSE
        # }
        
      }
        
      
      if (method.range == "convex.hull")
        p1 <-
          Convex.Hull.Poly(XY = XY[, c(2, 1)],
                            mode = mode,
                            proj_type = proj_type, 
                            exclude.area = exclude.area,
                            poly_exclude = country_map)
      
      # old <- function(XY) {
      #   p1 <-
      #     .Convex.Hull.Poly(XY[,c(2, 1)])
      #   eoo <- round(suppressWarnings(geosphere::areaPolygon(p1)),
      #         0)
      #   return(eoo)
      # }
      # 
      
      # 
      # library(sf)
      # 
      # microbenchmark::microbenchmark("old" = { eoo <- old(XY = XY) },
      #                "new" = {
      #                  eoo <- new(XY = XY)
      #                }, times = 400L)
      # 
      # microbenchmark::microbenchmark("old" = {
      #   p1 <-
      #     .Convex.Hull.Poly(XY[, c(2, 1)], mode = "spheroid")
      # },
      # "new" = {
      #   p1 <-
      #     .Convex.Hull.Poly(XY[, c(2, 1)], mode = "planar")
      # },
      # times = 400L)
      
      # if (exclude.area) {
      #   croped.EOO <- 
      #     .crop.poly(poly = p1, 
      #                crop = country_map)
      #   p1 <- croped.EOO[[2]]
      # }
      
      ## If exclude.area is TRUE
      # if (exclude.area) {
      #   EOO <-
      #     croped.EOO[[1]]
      # } else{
      #   EOO <-
      #     suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
      # }
      
      if(any(class(p1) == "SpatialPolygons") | any(class(p1) == "sfc") | any(class(p1) == "sf")) {
        # if (mode == "spheroid") {
        #   
        #   EOO <-
        #     as.numeric(st_area(p1)) / 1000000
        #   
        #   # EOO <-
        #   #   suppressWarnings(geosphere::areaPolygon(p1)) / 1000000
        # }
        # 
        # if (mode == "planar") {
        #   
        #   EOO <-
        #     as.numeric(st_area(p1)) / 1000000
        #   
        #   p1 <- 
        #     as(st_transform(p1, 4326), "Spatial")
        # }
        
        EOO <-
          as.numeric(st_area(p1)) / 1000000
        
        if (mode == "planar")
          p1 <-
            sf::st_transform(p1, 4326)
        
        p1$tax <- 
          Name_Sp
        
      } else  {
        
        EOO <- NA
        
      }
    # }
    
    OUTPUT <- list(EOO = EOO, spatial.polygon = p1)
  }
  
  digits <-
    c(6, 5, 4, 3, 2, 1, 0)[findInterval(OUTPUT$EOO, 
                                        c(0, 0.0001, 0.01, 0.1, 1, 10, 30000, Inf))]
  
  OUTPUT$EOO <- round(OUTPUT$EOO, digits)
  
  # if(verbose) cat(" ",paste(Name_Sp,"EOO comp."))
  
  return(OUTPUT)
}

