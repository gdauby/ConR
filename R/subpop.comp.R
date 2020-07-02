
#' @title Number of subpopulations
#'
#' @description Estimate the number of locations following the method **circular buffer method**
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY string, indicating the method used for estimating the number of locations. Either "fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"
#' @param Resol_sub_pop numeric. Defines in kilometers the radius of the circles around each occurrence
#' @param parallel logical, whether running in parallel. By default, it is FALSE
#' @param NbeCores string integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
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
subpop.comp <- function(XY, 
                        Resol_sub_pop = 5, 
                        proj_type = "cea",
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2) {
  
  proj_type <- proj_crs(proj_type = proj_type)
  
  list_data <- coord.check(XY = XY, listing = T, proj_type)
  
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else{
    `%d%` <- foreach::`%do%`
  }
  
  x <- NULL
  
  if(show_progress) {
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(list_data),
                            style = 3)
    
    progress <- function(n)
      utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{opts <- NULL}
  
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        utils::setTxtProgressBar(pb, x)
      
      res <- 
        subpop.estimation(
          XY = list_data[[x]], 
          Resol_sub_pop = Resol_sub_pop, 
          proj_type = proj_type
        )
      
      res
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  number_subpop <- 
    unlist(output[names(output) == "number_subpop"])
  poly <- 
    unlist(output[names(output) == "poly_subpop"])
  names(number_subpop) <-
    names(poly) <-
    gsub(pattern = " ",
         replacement = "_",
         names(list_data))
  
  # if (length(OUTPUT) == 1)
  #   OUTPUT <- OUTPUT[[1]]
  
  return(list(number_subpop = number_subpop, poly_subpop = poly))
}