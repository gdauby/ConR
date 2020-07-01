
#' @title Number of subpopulations
#'
#' @description Estimate the number of locations following the method **circular buffer method**
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY string, indicating the method used for estimating the number of locations. Either "fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"
#' @param Resol_sub_pop numeric. Defines in kilometers the radius of the circles around each occurrence
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
#' @references Rivers MC, Bachman SP, Meagher TR, Lughadha EN, Brummitt NA (2010) Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation 19: 2071-2085. doi: 10.1007/s10531-010-9826-9
#'
#' @return A list with one list for each taxa containing [[1]]Number of subpopulation and [[2]]SpatialPolygons.
#' 
#' @examples 
#' data(dataset.ex)
#' \dontrun{
#'subpop <- subpop.comp(dataset.ex, Resol_sub_pop = 5)
#'}
#'
#' 
#' @export
subpop.comp <- function(XY, Resol_sub_pop = 5) {
  
  # if (is.null(Resol_sub_pop))
  #   stop("Resol_sub_pop is missing, please provide a value")
  # 
  # if (any(is.na(XY[, c(1, 2)]))) {
  #   length(which(rowMeans(is.na(XY[, 1:2])) > 0))
  #   unique(XY[which(rowMeans(is.na(XY[, 1:2])) > 0), 3])
  #   print(
  #     paste(
  #       "Skipping",
  #       length(which(rowMeans(is.na(
  #         XY[, 1:2]
  #       )) > 0)) ,
  #       "occurrences because of missing coordinates for",
  #       paste(as.character(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
  #                                            0), 3])), collapse = " AND ")
  #     )
  #   )
  #   XY <- XY[which(!is.na(XY[, 1])), ]
  #   XY <- XY[which(!is.na(XY[, 2])), ]
  # }
  # 
  # if (any(XY[, 1] > 180) ||
  #     any(XY[, 1] < -180) ||
  #     any(XY[, 2] < -180) ||
  #     any(XY[, 2] > 180))
  #   stop("coordinates are outside of expected range")
  # 
  # colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
  # XY$tax <- as.character(XY$tax)
  # list_data <- split(XY, f = XY$tax)
  
  list_data <- coord.check(XY = XY)
  
  OUTPUT <-
    lapply(list_data, function(x)
      subpop.comp(XY = x, Resol_sub_pop = Resol_sub_pop))
  
  if (length(OUTPUT) == 1)
    OUTPUT <- OUTPUT[[1]]
  
  return(OUTPUT)
}