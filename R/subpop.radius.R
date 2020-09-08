#' @title Estimate Taxon Dispersal Radius
#'
#' @description Calculates the fraction of the maximum distance between
#'   species occurrences, which can be used as a proxy of species dispersal
#'   ability. This fraction can be used as the radius needed to apply the
#'   circular buffer method (Rivers et al 2010) to obtain subpopulations for
#'   conservation assessments of extinction risk.
#'
#' @param XY a data frame containing the geographical coordinates for each taxon.
#' @param factor.div. numeric. denominator value used to obtain the fraction of
#'   the maximum distance. Default to 10.
#' @param quant.max numeric. The upper-quantile of the inter-point distance to
#'   be considered as a threshold of maximum distance. Can very between 0 and 1.
#'   Default to 1.
#' @param proj_type   
#' 
#' @details 
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not matter}
#' 
#' \tabular{ccc}{
#'   [,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#'   [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#'   [,3] \tab tax \tab character or factor, taxa names\cr
#' }
#' 
#' @return The estimated radius in kilometers for each taxon.
#' 
#' @author Renato A. Ferreira de Lima & Gilles Dauby
#'   
#' @references Rivers MC, Bachman SP, Meagher TR, Lughadha EN, Brummitt NA
#'   (2010) Subpopulations, locations and fragmentation: applying IUCN red list
#'   criteria to herbarium specimen data. Biodiversity and Conservation 19:
#'   2071-2085. doi: 10.1007/s10531-010-9826-9   
#'   
#' @examples
#' 
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-28.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-17.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#'                    
#' subpop.radius(mydf)
#' subpop.radius(mydf, quant.max = 0.95)
#' subpop.radius(mydf, factor.div = 15)
#' 
#' @import sf
#' @import data.table
#' @importFrom stats quantile
#' 
#' @export subpop.radius
#' 
subpop.radius = function(XY,
                         factor.div = 10,
                         quant.max = 1) {
  
  XY <- 
    ConR:::coord.check(XY = XY, listing = FALSE)
  
  ## Getting the maximum inter-point distance
  #### GILLES: The function 'subpop.radius' is so small that I did not put 'f'
  # as a separate internal function
  f <- function(lat, lon) {
    
    x <- cbind.data.frame(lon, lat)
    
    if (dim(x)[1] > 3) {
      
      points_sf <- sf::st_as_sf(x, coords = c("lon", "lat"))
      sf::st_crs(points_sf) <- sf::st_crs(4326)
      
      dist <- sf::st_distance(points_sf)
      dist <- matrix(dist, nrow = dim(dist)[1], ncol = dim(dist)[2])
      dist <- dist[upper.tri(dist)]
      
      d.inter <- as.character(round(stats::quantile(
        as.double(dist / 1000),
        prob = quant.max,
        na.rm = TRUE
      ) / factor.div,
      3))
      
    } else {
      d.inter <- NA_character_
    }
    return(d.inter)
  }  

  ## Getting the maximum inter-point distance
  #### GILLES: I did not changed yet this part, which is done using
  #package data.table. Maybe use another approach to not import the
  #the package just for this small function? But I am not sure if
  #the foreach loop is really necessary.
  
  XY.dt <- data.table::data.table(XY)
  data.table::setkeyv(XY.dt, "tax") ## setting 'tax' as key/group to the data.table (makes computations faster)
  radius <- as.data.frame(XY.dt[ , f(ddlat, ddlon) , by = tax]) 
  names(radius) <- c("tax", "radius")
  
  return(radius)
}
