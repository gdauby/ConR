#' Internal function
#'
#' Get Occurrences Distance to a Polygon
#'
#' @param poly Spatial polygon
#' @param points XY data frame
#' @param proj_user projected coordinate system (in meters)
#' @param names_poly character string
#' @param names_taxa character string
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'#'
#' @param proj_type character string or numeric or object of CRS class, by default is "cea"
#' @param min.dist minimum tolerated distance between polygons and points.
#'   Default to 0.1 m.
#' @param value output value: proportional distance ("dist") or inside/outside
#'   the polygon ("flag")?

#' @details The spatial polygon must be a `SpatialPolygonsDataFrame` in which
#'   each polygon/feature is one taxon, an the data frame contains a column
#'   `tax` with the taxa name. The XY data frame has the same structure as other
#'   XY objects within `ConR` with the three first columns being `ddlat`,
#'   `ddlon` and `tax`, but also a column `classes` with the classes of
#'   conficence level, which is represented by numbers in increasing order of
#'   confidence. For instance, if confidence levels are "low" and "high", the
#'   corresponding values in columns `classes` should be 1 and 2.
#'   
#' @examples
#' 
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-28.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-17.2),
#'                    tax = rep("a", 7),
#'                    valid = c(c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE)),
#'                    stringsAsFactors = FALSE)
#' mydf$classes = as.double(mydf$valid)
#' shp <- EOO.computing(mydf[mydf$valid,], export_shp = TRUE)[[2]]
#' shp$a <- data.frame(tax = "a", stringsAsFactors = FALSE)
#' plot(mydf[,2:1])
#' plot(sf::st_geometry(shp), add=TRUE)
#' over.valid.poly(shp, mydf, names_poly = "a", names_taxa = "a")
#' over.valid.poly(shp, mydf, names_poly = "a", names_taxa = "a", value = "flag")  
#' over.valid.poly(shp, mydf, names_poly = "a", names_taxa = "b")
#' 
#' @import sf
#' @importFrom fields rdist
#' 
over.valid.poly <- function(poly,
                            points,
                            names_poly = NULL,
                            names_taxa = NULL,
                            mode = "spheroid",
                            proj_type = "cea",
                            min.dist = 0.1,
                            value = "dist") {
  
  poly.tax <- 
    poly[grep(names_taxa, names_poly),]
  
  if( nrow(poly.tax) > 0) {
    # if (is.null(proj_user)) {
    #   proj_user <- 3857
    #   warning("no projected coordinate reference system provided by the user: assuming WGS 84 Pseudo-Mercator (see https://epsg.io)"
    #   )
    # }
    
    # poly_sf <- sf::st_as_sf(poly)
    # if (is.na(sf::st_crs(poly_sf)[[1]]))
    #   sf::st_crs(poly_sf) <- proj_user
    
    points_sf <- sf::st_as_sf(points, coords = c("ddlon", "ddlat"))
    # if (is.na(sf::st_crs(points_sf)[[1]]))
    sf::st_crs(points_sf) <- sf::st_crs(4326)
    
    # points_sf <- sf::st_transform(points_sf, crs = proj_user)
    # poly_sf <- sf::st_transform(poly_sf, crs = proj_user)
    
    # Solution for multiple species at once# solution for single species at a time
    # if (nrow(poly.tax) == 1) {
    #   
    #   # dist_poly <- sf::st_distance(points_sf, poly_sf)
    #   # dists <- round(as.double(dist_poly[,1]), 2)
    #   
    # 
    # # Solution for multiple species EOO at once
    # } else {
    
    dist_poly <- sf::st_distance(points_sf, poly.tax)
    dist_poly <- matrix(dist_poly, nrow = dim(dist_poly)[1], ncol = dim(dist_poly)[2])
    colnames(dist_poly) <- names_poly[grep(names_taxa, names_poly)]
    dists <- dist_poly[,1]
    # row.names(dist_poly) <- as.character(points$tax)
    # # j <- match(rownames(dist_poly), colnames(dist_poly))
    # dists <- round(diag(dist_poly[,j]),2)
    # dists <- round(diag(dist_poly),2)
    # }
    
    flag <- 
      data.frame(points,
                 flag = dists < min.dist)
    # 
    # flag <- dists < min.dist
    
    
    if(value == "dist") {
      
      dist_ch_strict <- 
        as.vector(sf::st_distance(points_sf[points_sf$classes == 1,]))
      dist_ch_strict <- 
        dist_ch_strict[dist_ch_strict > 0]
      dist_ch_strict_quant <- 
        stats::quantile(dist_ch_strict, 0.95)
      
      rel_dist <- dists / dist_ch_strict_quant
      
      df_weights <- data.frame(points,
                               rel_dist = round(rel_dist, 3))
      
      # coords <- as.data.frame(sf::st_coordinates(points_sf))
      # true.ids <- points_sf$classes >= max(points_sf$classes, na.rm = TRUE)
      # coords.true <- coords[true.ids,]
      # coords.true$tax <- points_sf$tax[true.ids]
      # coords.true <- coords.true[!is.na(coords.true[, 1]), ]
      
      # med_dists <- tapply(1:nrow(coords.true), coords.true$tax, 
      #                     function(x) mean(dist(coords.true[x, 1:2]), na.rm = TRUE))
      # med_dists <- tapply(1:nrow(coords.true), coords.true$tax, 
      #                                 function(x) mean(as.matrix(distances::distances(coords.true[x, 1:2])), na.rm = TRUE))
      # med_dists <- tapply(1:nrow(coords.true), coords.true$tax, 
      #                                 function(x) mean(fields::rdist(coords.true[x, 1:2], compact = TRUE), na.rm = TRUE))
      # rob <- tapply(1:nrow(coords.true), coords.true$tax, 
      #               function(x) robustbase::covMcd(coords.true[x, 1:2], alpha = 1/2, tol=1e-20))
      # NOT WORKING...
      # maha <- tapply(1:nrow(coords), points_sf$tax, 
      #                function(x) sqrt(stats::mahalanobis(coords[x, 1:2], center = rob[x]$center, cov = rob[x]$cov, tol=1e-20)))
      #max_dists <- tapply(1:nrow(coords.true), coords.true$tax,
      #                                 function(x) quantile(distances::distance_matrix(distances::distances(coords.true[x, 1:2])), prob=0.95 , na.rm = TRUE))
      ## MAYBE USE FUNCTION rdist.eart INSTEAD
      # max_dists <- tapply(1:nrow(coords.true), coords.true$tax, 
      #                      function(x) quantile(fields::rdist(coords.true[x, 1:2], compact = TRUE), prob=0.95 , na.rm = TRUE))
      # tmp <- dplyr::left_join(data.frame(dist = dists, tax = points_sf$tax, stringsAsFactors = FALSE),
      #                         data.frame(
      #                           #med_dists = as.double(med_dists),
      #                           inter_dists = as.double(max_dists),
      #                           tax = names(max_dists), stringsAsFactors = FALSE), 
      #                         by = 'tax')
      # dist <- round(tmp$d / tmp$inter_dists, 5)
      return(df_weights)
    }
    
    
    if(value == "flag") {
      
      return(flag)
      
    }
    
  } else {
    #return(NA)
    return(rep(NA, nrow(points)))
  }
  
  
  
}  
