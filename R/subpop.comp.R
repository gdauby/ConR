
#' @title Estimate the number of subpopulations for one or multiple species
#'
#' @description 
#' `r lifecycle::badge("stable")`
#'  Estimate the number of subpopulations following the method
#'   **circular buffer method** (overlapping buffered circles form a single 
#'   subpopulation)
#'
#' @author Gilles Dauby & Renato A. Ferreira de Lima
#'
#' @param XY a data frame containing the geographical coordinates for each taxon
#'   (see Details).
#' @param resol_sub_pop a value defining the radius of the circles around each
#'   occurrence (in kilometres) or data frame vector containing a column 'tax' 
#'   with the taxa names and a column 'radius' with the species-specific radius
#'   (in kilometre as well). Typically, this data frame is the output of
#'   ```ConR``` function ```subpop.radius```.
#' @param export_shp logical. Whether the resulting shapefiles should be
#'   exported. FALSE by default.
#' @inheritParams activate_parallel
#' @param show_progress logical. Whether progress informations should displayed. TRUE by default
#' @inheritParams proj_crs
#' 
#' @details 
#' `XY` as a [data.frame][base::data.frame()] should have the following structure:
#' 
#' **It is mandatory to respect field positions, but field names do not matter**
#' 
#' \enumerate{
#'   \item The first column is contains numeric value i.e. latitude in decimal degrees
#'   \item The second column is contains numeric value i.e. longitude in decimal degrees
#'   \item The third column is contains character value i.e. the names of the species
#' }
#' 
#' @references Rivers MC, Bachman SP, Meagher TR, Lughadha EN, Brummitt NA
#'   (2010) Subpopulations, locations and fragmentation: applying IUCN red list
#'   criteria to herbarium specimen data. Biodiversity and Conservation 19:
#'   2071-2085. doi: 10.1007/s10531-010-9826-9
#'
#' @return 
#' If `export_shp` is TRUE,
#' \enumerate{
#'   \item number_subpop a numeric vector of AOO estimates for each taxa
#'   \item poly_subpop a simple feature collection
#' }
#' 
#' a `Simple feature collection` with as many MULTIPOLYGON as taxa.
#' If `export_shp` is FALSE, a vector with estimated number of subpopulation per
#' taxa.
#' 
#' @examples 
#' data(dataset.ex)
#'
#' subpop.comp(dataset.ex, resol_sub_pop = 5)
#' rad.df <- data.frame(
#'     tax = unique(dataset.ex$tax),
#'     radius = seq(3,13, by=2),
#'     stringsAsFactors = FALSE
#'   )
#' subpop.comp(dataset.ex, resol_sub_pop = rad.df)
#' subpop.comp(dataset.ex, resol_sub_pop = rad.df, export_shp = TRUE)
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom parallel stopCluster
#' 
#' @export subpop.comp
#' 
subpop.comp <- function(XY, 
                        resol_sub_pop = NULL, 
                        proj_type = "cea",
                        export_shp = FALSE,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2) {
  
  if (is.null(resol_sub_pop)) 
    stop(" is missing, please provide a value for all species or a data frame with species-specific values")
  
  proj_type <- 
    proj_crs(proj_type = proj_type)
  
  if ("data.frame" %in% class(resol_sub_pop)) {
    XY <- merge(XY, resol_sub_pop, 
                by = "tax", all.X = TRUE, sort = FALSE)
    XY <- XY[,c("ddlat", "ddlon", "tax", "radius")]
  } else {
    XY$radius <- resol_sub_pop   
  }

  list_data <-
    coord.check(XY = XY, listing = TRUE, proj_type = proj_type)
  
  cl <- activate_parallel(parallel = parallel, NbeCores = NbeCores)
  `%d%` <- c_par(parallel = parallel)
  
  pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_data[[1]]))
  opts <- pro_res$opts
  pb <- pro_res$pb
  
  output <-
    foreach::foreach(
      x = 1:length(list_data[[1]]),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        utils::setTxtProgressBar(pb, x)
      
      res <- 
        subpop.estimation(
          XY = list_data[[1]][[x]], 
          resol_sub_pop = unique(list_data[[1]][[x]]$radius),
          proj_type = proj_type,
          export_shp = export_shp
        )
      
      if (export_shp) {
        names(res) <- c("subpop", "spatial")
        names(res)[1] <- list_data[[1]][[x]]$tax[1]
        res$spatial <- cbind(res$spatial, tax = list_data[[1]][[x]]$tax[1])
        
      } else {
        
        names(res)[1] <- list_data[[1]][[x]]$tax[1]
        
      }
      
      res
    }
  
  if(parallel) parallel::stopCluster(cl)
  if(show_progress) close(pb)
  
  if (export_shp) {
    
    number_subpop <-
      data.frame(tax =  names(output[names(output) != "spatial"]),
                 subpop =  as.numeric(unlist(output[names(output) != "spatial"])))

    shapes <- output[names(output) == "spatial"]
    shapes <- do.call('rbind', shapes)
    row.names(shapes) <- 1:nrow(shapes)
    
    shapes <- st_transform(shapes, 4326)
    
  } else {
    
    number_subpop <-
      data.frame(tax =  names(output),
                 subpop =  as.numeric(unlist(output)))
  }
  
  if (!export_shp) return(number_subpop)
  
  if (export_shp) return(list(number_subpop = number_subpop,
                              poly_subpop = shapes))
  
}
