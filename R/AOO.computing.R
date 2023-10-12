
#' @title Area of occupancy
#'
#' @description Compute areas of occupancy (AOO) for multiple taxa in square kilometers
#'
#' @author Gilles Dauby \email{gilles.dauby@@ird.fr}
#'
#' @param XY [data.frame][base::data.frame()] see Details.
#' @param Cell_size_AOO numeric, by default is 2. Value indicating the grid size
#'   in kilometres used for estimating Area of Occupancy.
#' @param nbe.rep.rast.AOO numeric, by default is 0. Indicate the number of
#'   raster with random starting position used for estimating the AOO. If 0 but
#'   some translation of the raster are still done.
#' @param parallel logical, by default is FALSE. Whether running in parallel.
#' @param NbeCores integer, by default is 2. Register the number of cores for
#'   parallel execution. Only used if parallel is TRUE.
#' @param show_progress logical, by default is TRUE. Whether a progress bar
#'   during computation is shown.
#' @param export_shp logical, by default is FALSE. Whether a shapefile of
#'   occupied cells should be exported.
#' @param proj_type character or numeric, by default is "cea", see Details.
#' 
#' @details 
#' # Input data
#' **XY** as a [data.frame][base::data.frame()] should have the following structure:
#' 
#' **It is mandatory to respect field positions, but column names do not matter**
#' 
#' \enumerate{
#'   \item The first column is contains numeric value i.e. latitude in decimal degrees
#'   \item The second column is contains numeric value i.e. longitude in decimal degrees
#'   \item The third column is contains character value i.e. the names of the species
#' }
#' 
#' See Examples.
#' 
#' # Iteration to get the minimal AOO
#' The argument `nbe.rep.rast.AOO` should ideally be higher than 10 for
#' increasing the chance to get the minimal number of occupied cell. However,
#' increasing `nbe.rep.rast.AOO` also increases the computing time. Note that if
#' `nbe.rep.rast.AOO = 0`, several translations of the grid overlaying
#' occurrences are still conducted
#' 
#' # proj_type
#' 
#' See [proj_type()]
#' 
#' 
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges,
#'   Journal of Applied Ecology, 49 1-9
#'
#' @return 
#' If `export_shp` if FALSE a vector of AOO estimates for each taxa
#' If `export_shp` if TRUE a list with two elements
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
#'AOO <- AOO.computing(dataset.ex, nbe.rep.rast.AO = 10)
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
    coord.check(XY = XY, proj_type = proj_type, cell_size = Cell_size_AOO, check_eoo = FALSE)
  
  issue_close_to_anti <- list_data$issue_close_to_anti
  list_data <- list_data$list_data
  
  res_df <-
    data.frame(aoo =  rep(NA, length(list_data)), 
               issue_aoo = rep(NA, length(list_data)))
  row.names(res_df) <- names(list_data)
  
  if (length(issue_close_to_anti) > 0) {
    
    list_data <- list_data[-issue_close_to_anti]
    
  }
  
  if (length(list_data) > 0) {
    
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
        
        # for (x in 1:length(list_data)) {
        # print(x)
        res <- AOO.estimation(
          coordEAC = list_data[[x]],
          cell_size = Cell_size_AOO,
          nbe_rep = nbe.rep.rast.AOO,
          export_shp = export_shp,
          proj_type = proj_type
        )
        
        if (export_shp) {
          names(res) <- c("aoo", "spatial")
          res$spatial <- cbind(res$spatial, tax = list_data[[x]]$tax[1])[, c("tax", "geometry")]
          # names(res)[1] <- list_data[[x]]$tax[1]
          # res$spatial <- cbind(res$spatial, tax = list_data[[x]]$tax[1])
          # res$spatial <- res$spatial[,-which(colnames(res$spatial) == "lyr.1")]
        } 
        # else {
        #   names(res)[1] <- list_data[[x]]$tax[1]
        # }
        names(res)[1] <- list_data[[x]]$tax[1]
        
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)

    res <- unlist(output[names(output) != "spatial"])
    
    res_df[which(row.names(res_df) %in% names(res)), 1] <-
      res
    
    if (length(issue_close_to_anti) > 0)
      res_df[issue_close_to_anti, 2] <-
      "AOO could not computed because grid cells would overlap with antimeridian"
    
    
    
    if(export_shp) {
      
      
      # res <- unlist(output[names(output) == "aoo"])
      # names(res) <- names(list_data)
      
      shapes <- output[names(output) == "spatial"]
      shapes <- do.call('rbind', shapes)
      row.names(shapes) <- 1:nrow(shapes)
      
      # names(shapes) <- names(list_data)
      
    }
    
  } else {
    
    if (export_shp)
      shapes <- NA
    
  }
  
  if(!export_shp) 
    return(res_df)
  
  if(export_shp) 
    return(list(AOO = res_df, 
                AOO_poly = shapes))
  
}
