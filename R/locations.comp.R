#' @title Number of locations
#'
#' @description Estimate the number of locations for multiple taxa
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY data.frame, see details
#' @param method string, indicating the method used for estimating the number of locations. Either "fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"
#' @param nbe_rep numeric , indicate the number of raster with random starting position for estimating the number of locations By default, it is 0 but some minimal translation of the raster are still done
#' @param protec.areas \code{SpatialPolygonsDataFrame}, shapefile with protected areas. If provided, this will be taken into account for calculating number of location (see Details and \code{method_protected_area}). By default, no shapefile is provided
#' @param Cell_size_locations numeric, value indicating the grid size in kilometers used for estimating the number of location. By default, equal to 10
#' @param method_protected_area string, by default is "no_more_than_one"", which means occurrences within protected areas (if provided) will not be taken into account for estimating the number of locations following the grid system, see Details. By default, it is "no_more_than_one"
#' @param ID_shape_PA string, indicating the field name of \code{protec.areas} with ID of the \code{SpatialPolygonsDataFrame} of protected areas
#' @param Rel_cell_size numeric, if \code{method_locations="sliding scale"}, \code{Cell_size_locations} is ignored and the resolution is given by the maximum distance separating two occurrences multiplied by \code{Rel_cell_size}. By default, it is 0.05
#' @param parallel logical, whether running in parallel. By default, it is FALSE
#' @param NbeCores string integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
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
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9
#'
#' @return A list with one list for each species containing [[1]]SpatialPolygonDataframe and [[2]]vector of the number of location.
#' 
#' @examples 
#' data(dataset.ex)
#' \dontrun{
#'locations <- locations.comp(dataset.ex)
#'}
#'
#'# This would estimate the number of locations for all taxa by overlaying 
#'# randomly a grid 100 times. For each taxa, the minimum value is kept
#' \dontrun{
#'locations <- locations.comp(dataset.ex, nbe_rep = 100)
#'}
#'
#' @importFrom rgdal project
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export
locations.comp <- function(XY,
                           method = "fixed_grid",
                           nbe_rep = 0,
                           protec.areas = NULL,
                           Cell_size_locations = 10,
                           method_protected_area = "no_more_than_one",
                           ID_shape_PA = "WDPA_PID",
                           Rel_cell_size = 0.05,
                           parallel = FALSE,
                           NbeCores = 2,
                           show_progress = TRUE) {
  
  if (!any(class(XY) == "data.frame"))
    XY <- as.data.frame(XY)
  if (any(XY[, 2] > 180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 1] > 180))
    stop("coordinates are outside of expected range")
  
  if(method != "fixed_grid" & method != "sliding scale")
    stop("method should be fixed_grid sliding scale")
  
  projEAC <- .proj_crs()
  
  coordEAC <-
    data.frame(matrix(unlist(
      rgdal::project(as.matrix(XY[, c(2, 1)]),
                     proj = as.character(projEAC), inv =
                       FALSE)
    ),
    ncol = 2),
    tax = XY[, 3])
  
  
  ## if any missing coordinates
  if (any(is.na(coordEAC[, c(1:2)]))) {
    print(
      paste(
        "Skipping",
        length(which(rowMeans(
          is.na(coordEAC[, 1:2])
        ) > 0)),
        "occurrences because of missing coordinates for",
        paste(as.character(unique(coordEAC[which(rowMeans(is.na(coordEAC[, 1:2])) >
                                                   0), 3])), collapse = " AND ")
      )
    )
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 1])),]
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 2])),]
  }
  
  coordEAC$tax <- as.character(coordEAC$tax)
  list_data <- split(coordEAC, f = coordEAC$tax)
  
  
  crs_proj <- projEAC
  
  ## geographical distances for all pairs of occurrences
  
  if (is.null(protec.areas)) {
    if (nrow(coordEAC) > 1)
      pairwise_dist <- stats::dist(coordEAC[, 1:2],  upper = F)
    
    ## resolution definition
    if (any(method == "fixed_grid"))
      Resolution <- Cell_size_locations
    if (any(method == "sliding scale")) {
      if (nrow(coordEAC) > 1) {
        Resolution <- max(pairwise_dist) * Rel_cell_size / 1000
      } else{
        Resolution <- 10
      }
    }
    
    if (parallel) {
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      
      # registerDoParallel(NbeCores)
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
        
        res <- .cell.occupied(
          size = Resolution,
          coord = list_data[[x]],
          nbe_rep = nbe_rep
        )
        
        names(res) <- c("spatial", "nbe_occ")
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    Locations <- unlist(output[names(output) == "nbe_occ"])
    r2 <- unlist(output[names(output) == "spatial"])
    names(Locations) <-
      names(r2) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
  }
  
  if (!is.null(protec.areas)) {
    ### Taking into account Protected Areas if provided
    DATA_SF <- as.data.frame(XY[, 1:2])
    colnames(DATA_SF) <- c("ddlat", "ddlon")
    sp::coordinates(DATA_SF) <-  ~ ddlon + ddlat
    
    
    sp::proj4string(DATA_SF) <- sp::CRS(SRS_string='EPSG:4326')
    
    # raster::crs(DATA_SF) <- raster::crs(protec.areas)
    
    Links_NatParks <- sp::over(DATA_SF, protec.areas)
    
    coordEAC_pa <- 
      coordEAC[!is.na(Links_NatParks[, 1]),]
    coordEAC_pa <-
      cbind(coordEAC_pa,
            id_pa =
              Links_NatParks[which(!is.na(Links_NatParks[, 1])), 
                             ID_shape_PA])
    
    LocNatParks <-
      vector(mode = "numeric", length = length(list_data))
    names(LocNatParks) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
    if (nrow(coordEAC_pa) > 0) {
      if (method_protected_area == "no_more_than_one") {
        ## if method is 'no_more_than_one' the number of location is the number of occupied protected areas
        
        loc_pa <-
          by(
            coordEAC_pa[, c("tax", "id_pa")],
            coordEAC_pa[, "tax"],
            FUN = function(x)
              length(unique(x$id_pa))
          )
        names(loc_pa) <-
          gsub(pattern = " ",
               replacement = "_",
               names(loc_pa))
        LocNatParks[names(LocNatParks) %in% names(loc_pa)] <-
          loc_pa
        r2_PA <- NA
        
      } else{
        coordEAC_pa$tax <- as.character(coordEAC_pa$tax)
        list_data_pa <- split(coordEAC_pa, f = coordEAC_pa$tax)
        
        ## geographical distances for all pairs of occurrences
        if (nrow(coordEAC_pa) > 1)
          pairwise_dist_pa <- stats::dist(coordEAC_pa[, 1:2],  upper = F)
        
        ## resolution definition
        if (any(method == "fixed_grid"))
          Resolution <- Cell_size_locations
        if (any(method == "sliding scale")) {
          if (nrow(coordEAC_pa) > 1) {
            Resolution <- max(pairwise_dist_pa) / 1000 * Rel_cell_size
          } else{
            Resolution <- 10
          }
        }
        
        if (parallel) {
          cl <- snow::makeSOCKcluster(NbeCores)
          doSNOW::registerDoSNOW(cl)
          
          # registerDoParallel(NbeCores)
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
            x = 1:length(list_data_pa),
            .combine = 'c',
            .options.snow = opts
          ) %d% {
            
            if (!parallel & show_progress)
              utils::setTxtProgressBar(pb, x)
            
            res <- .cell.occupied(
              size = Resolution,
              coord = list_data_pa[[x]],
              nbe_rep = nbe_rep
            )
            
            names(res) <- c("spatial", "nbe_occ")
            res
          }
        
        if(parallel) snow::stopCluster(cl)
        if(show_progress) close(pb)
        
        loc_pa <- unlist(output[names(output) == "nbe_occ"])
        r2_PA <- unlist(output[names(output) == "spatial"])
        names(loc_pa) <-
          names(r2_PA) <-
          gsub(pattern = " ",
               replacement = "_",
               names(list_data_pa))
        LocNatParks[names(LocNatParks) %in% names(loc_pa)] <-
          loc_pa
      }
    } else{
      r2_PA <- NA
    }
    
    coordEAC_not_pa <- coordEAC[is.na(Links_NatParks[, 1]),]
    LocOutNatParks <-
      vector(mode = "numeric", length = length(list_data))
    names(LocOutNatParks) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
    if (nrow(coordEAC_not_pa) > 0) {
      coordEAC_not_pa$tax <- as.character(coordEAC_not_pa$tax)
      list_data_not_pa <-
        split(coordEAC_not_pa, f = coordEAC_not_pa$tax)
      
      ## geographical distances for all pairs of occurrences
      if (nrow(coordEAC_pa) > 1)
        pairwise_dist_not_pa <- stats::dist(coordEAC_not_pa[, 1:2],  upper = F)
      
      ## resolution definition
      if (any(method == "fixed_grid"))
        Resolution <- Cell_size_locations
      if (any(method == "sliding scale")) {
        if (nrow(coordEAC_pa) > 1) {
          Resolution <- max(pairwise_dist_not_pa) * Rel_cell_size
        } else{
          Resolution <- 10
        }
      }
      
      if (parallel) {
        cl <- snow::makeSOCKcluster(NbeCores)
        doSNOW::registerDoSNOW(cl)
        
        # registerDoParallel(NbeCores)
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
        foreach::foreach(x = 1:length(list_data_not_pa),
                         .combine = 'c',
                         .options.snow = opts) %d% {
                           
                           if (!parallel & show_progress)
                             utils::setTxtProgressBar(pb, x)
                           
                           res <- .cell.occupied(
                             size = Resolution,
                             coord = list_data_not_pa[[x]],
                             nbe_rep = nbe_rep
                           )
                           
                           names(res) <- c("spatial", "nbe_occ")
                           res
                         }
      
      if(parallel) snow::stopCluster(cl)
      if(show_progress & show_progress) close(pb)
      
      loc_not_pa <- unlist(output[names(output) == "nbe_occ"])
      r2 <- unlist(output[names(output) == "spatial"])
      names(loc_not_pa) <-
        names(r2) <-
        gsub(pattern = " ",
             replacement = "_",
             names(list_data_not_pa))
      LocOutNatParks[names(LocOutNatParks) %in% names(loc_not_pa)] <-
        loc_not_pa
      
    } else{
      r2 <- NA
    }
    
  }
  
  if (!is.null(protec.areas))
    return(list(r2, r2_PA, LocNatParks, LocOutNatParks))
  if (is.null(protec.areas))
    return(list(r2, Locations))
  
}