
#' @title Area of occupancy
#'
#' @description Compute areas of occupancy (AOO) for multiple taxa in square
#'   kilometers
#'
#' @author Gilles Dauby \email{gilles.dauby@@ird.fr}
#'
#' @param XY [data.frame][base::data.frame()] see Details.
#' @param cell_size_AOO numeric, by default is 2. Value indicating the grid size
#'   in kilometres used for estimating Area of Occupancy.
#' @param nbe.rep.rast.AOO numeric, by default is 0. Indicate the number of
#'   raster with random starting position used for estimating the AOO. If 0 but
#'   some translation of the raster are still done.
#' @inheritParams activate_parallel
#' @param show_progress logical. Whether progress informations should displayed.
#'   TRUE by default
#' @param export_shp logical, by default is FALSE. Whether a shapefile of
#'   occupied cells should be exported.
#' @inheritParams proj_crs
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
#' 
#' 
#' See `proj_type`
#' 
#' 
#' @references Gaston & Fuller 2009 The sizes of species' geographic ranges,
#'   Journal of Applied Ecology, 49 1-9
#'
#' @return 
#' If `export_shp` if FALSE (the default) a numeric vector of AOO estimates for each taxa
#' If `export_shp` if TRUE a list with two elements
#' \enumerate{
#'   \item a vector of AOO estimates for each taxa
#'   \item a list of simple feature for each taxa
#' }
#'   
#' @examples 
#' data(dataset.ex)
#'
#' AOO <- AOO.computing(dataset.ex)
#'
#'
#'# This would estimate AOO for all taxa by overlaying randomly a 
#'# grid 3 times. For each taxa, the minimum value is kept
#'
#' AOO <- AOO.computing(dataset.ex, nbe.rep.rast.AO = 3)
#'
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel stopCluster
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export
AOO.computing <- function(XY,
                          cell_size_AOO = 2,
                          nbe.rep.rast.AOO = 0,
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          export_shp = FALSE,
                          proj_type = "cea"
) {
  
  proj_type <- proj_crs(proj_type = proj_type)
  
  list_data <- 
    coord.check(XY = XY, proj_type = proj_type, cell_size = cell_size_AOO, check_eoo = FALSE)
  
  issue_close_to_anti <- list_data$issue_close_to_anti
  list_data <- list_data$list_data
  
  res_df <-
    data.frame(tax = names(list_data),
               aoo =  rep(NA, length(list_data)), 
               issue_aoo = rep(NA, length(list_data)))
  
  # res_df <-
  #   data.frame(aoo =  rep(NA, length(list_data)), 
  #              issue_aoo = rep(NA, length(list_data)))
  # row.names(res_df) <- names(list_data)
  
  if (length(issue_close_to_anti) > 0) {
    
    list_data <- list_data[-issue_close_to_anti]
    
  }
  
  if (length(list_data) > 0) {
    
    cl <- activate_parallel(parallel = parallel, NbeCores = NbeCores)
    `%d%` <- c_par(parallel = parallel)
    pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_data))
    opts <- pro_res$opts
    pb <- pro_res$pb
    
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
          cell_size = cell_size_AOO,
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
    
    if(parallel) parallel::stopCluster(cl)
    if(show_progress) close(pb)
    
    res <- unlist(output[names(output) != "spatial"])
    
    res_df[res_df$tax %in% names(res), 2] <-
      res
    
    if (length(issue_close_to_anti) > 0)
      res_df[issue_close_to_anti, 3] <-
      "AOO could not computed because grid cells would overlap with antimeridian"
    
    if (export_shp) {
      
      shapes <- output[names(output) == "spatial"]
      shapes <- do.call('rbind', shapes)
      row.names(shapes) <- 1:nrow(shapes)
      
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


#' @importFrom foreach %dopar% %do% foreach
#' @keywords internal
#' 
AOO.estimation <- function(coordEAC,
                           cell_size = 2,
                           nbe_rep = 0,
                           export_shp = FALSE,
                           proj_type = proj_type
) {
  
  res <-
    cell.occupied(
      nbe_rep = nbe_rep,
      size = cell_size,
      coord = coordEAC[,c(2, 1)],
      export_shp = export_shp,
      proj_type = proj_type
    )
  
  
  
  AOO <- res[[2]] * cell_size * cell_size
  
  if (export_shp)
    return(list(AOO = AOO, poly_AOO = res[[1]]))
  
  if (!export_shp)
    return(AOO)
  
}
