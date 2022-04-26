
#' @title Area of occupancy
#'
#' @description Compute areas of occupancy (AOO) for multiple taxa in square kilometers
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY \code{"dataframe"} see Details
#' @param Cell_size_AOO numeric, value indicating the grid size in kilometers used for estimating Area of Occupancy.  By default, equal to 2 km (i.e. 4 km2 grid cells)
#' @param nbe.rep.rast.AOO numeric , indicate the number of raster with random starting position for estimating the AOO. By default, it is 0 but some minimal translation of the raster are still done
#' @param parallel logical, whether running in parallel. By default, it is FALSE
#' @param NbeCores string integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
#' @param export_shp logical, whether a shapefile of occupied cells should be exported. By default, it is FALSE
#' @param proj_type character string or numeric or object of CRS class, by default is "cea"
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
#' The argument \code{nbe.rep.rast.AOO} should ideally be higher than 20 for increasing 
#' the chance to get the minimal number of occupied cell. Increasing \code{nbe.rep.rast.AOO} however 
#' also increase the computing time. So this is a trade-off that depends on the importance to 
#' get the minimal AOO and the size of the dataset.
#' 
#' 
#' 
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9
#'
#' @return 
#' If \code{export_shp} if FALSE a vector of AOO estimates for each taxa
#' If \code{export_shp} if TRUE a list with two elements
#' \enumerate{
#'   \item a vector of AOO estimates for each taxa
#'   \item a list of SpatialPolygonsDataFrame for each taxa
#' }
#'   
#' @examples 
#' data(dataset.ex)
#' \dontrun{
#'AOO <- AOO.computing(dataset.ex)
#'}
#'
#'# This would estimate AOO for all taxa by overlaying randomly a 
#'# grid 100 times. For each taxa, the minimum value is kept
#' \dontrun{
#'AOO <- AOO.computing(dataset.ex, nbe.rep.rast.AO = 100)
#'}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export
AOO.computing <- function(XY,
                          Cell_size_AOO = 2,
                          nbe.rep.rast.AOO = 0,
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          export_shp = FALSE,
                          proj_type = "cea"
) {
  
  proj_type <- proj_crs(proj_type = proj_type)
  
  list_data <- 
    coord.check(XY = XY, proj_type = proj_type)
  
  if(parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    `%d%` <- foreach::`%dopar%`
  }else{
    `%d%` <- foreach::`%do%`
  }
  
  
  x <- NULL
  if(show_progress) {
    pb <-
      txtProgressBar(min = 0,
                            max = length(list_data),
                            style = 3)
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {opts <- NULL}
  
  
  # print(proj_type)
  
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c', .options.snow = opts
    ) %d% {
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      # source("./R/IUCNeval.functionv11.R")
      
      res <- AOO.estimation(
        coordEAC = list_data[[x]],
        cell_size = Cell_size_AOO,
        nbe_rep = nbe.rep.rast.AOO,
        export_shp = export_shp,
        proj_type = proj_type
      )
      
      if (export_shp) {
        names(res) <- c("aoo", "spatial")
        names(res)[1] <- list_data[[x]]$tax[1]
        res$spatial <- cbind(res$spatial, tax = list_data[[x]]$tax[1])
        res$spatial <- res$spatial[,-which(colnames(res$spatial) == "lyr.1")]
      } else {
        names(res)[1] <- list_data[[x]]$tax[1]
      }
      
      res
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  if(!export_shp) {
    
    res_df <-
      data.frame(aoo =  unlist(output[names(output) != "spatial"]))
    
  }
  
  if(export_shp) {
    
    res_df <-
      data.frame(aoo =  unlist(output[names(output) != "spatial"]))
    
    # res <- unlist(output[names(output) == "aoo"])
    # names(res) <- names(list_data)
    
    shapes <- output[names(output) == "spatial"]
    shapes <- do.call('rbind', shapes)
    row.names(shapes) <- 1:nrow(shapes)
    
    # names(shapes) <- names(list_data)
    
  }
  
  if(!export_shp) 
    return(res_df)
  
  if(export_shp) 
    return(list(AOO = res_df, 
                AOO_poly = shapes))
  
}