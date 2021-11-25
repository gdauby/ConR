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
#' Locations are estimated by overlaying a grid of a given resolution (see \code{Cell_size_locations} for
#' specifying the resolution). The number of locations is simply the number of
#' occupied locations. Note that the grid position is overlaid in order to
#' minimize the number of locations (several translation of the grid are
#' performed and the one providing the minimum number of occupied cells is
#' provided).
#' 
#' If \code{method} is "fixed_grid" as it is by default, the resolution is fixed and determined 
#' by the argument \code{Cell_size_locations}.
#' If \code{method} is "sliding scale", the resolution is defind as 1/x*max.dist where max.dist is the maximum distance between any pairs of occurrences 
#' and x is a defined parameter. 1/x is defined by \code{Rel_cell_size} argument and is 0.05 by default. 
#' See Rivers M.C. et al. (2010) for more information on the methods.
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
                           show_progress = TRUE,
                           proj_type = "cea") {
  
  
  proj_type <- proj_crs(proj_type = proj_type)
  
  list_data <- 
    coord.check(XY = XY, proj_type = proj_type)
  
  ## geographical distances for all pairs of occurrences
  
  if (is.null(protec.areas)) {
    
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
    
    if (show_progress) {
      pb <-
        utils::txtProgressBar(min = 0,
                              max = length(list_data),
                              style = 3)
      
      progress <- function(n)
        utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } else{
      opts <- NULL
    }
    
    output <-
      foreach::foreach(
        x = 1:length(list_data),
        .combine = 'c',
        .options.snow = opts
      ) %d% {
        
        if (!parallel & show_progress)
          utils::setTxtProgressBar(pb, x)
        
        res <- 
          Locations.estimation(
          coordEAC = list_data[[x]], 
          cell_size = Cell_size_locations, 
          export_shp = TRUE, 
          proj_type = proj_type, 
          method = method,
          nbe_rep = nbe_rep,
          Rel_cell_size = Rel_cell_size
        )
        
        names(res) <- c("nbe_occ", "spatial")
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    Locations <- unlist(output[names(output) == "nbe_occ"])
    r2 <- output[names(output) == "spatial"]
    names(Locations) <-
      names(r2) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
  }
  
  if (!is.null(protec.areas)) {
    ### Taking into account Protected Areas if provided
    XY_ID <- 
      data.frame(XY, ID_prov_data = seq(1, nrow(XY), 1))
    
    DATA_SF <- st_as_sf(coord.check(XY = XY_ID, listing = F, proj_type = proj_type), 
                        coords = c(2, 1), crs = proj_type)
    
    if(!any(class(protec.areas) == "sf"))
      protec.areas <- 
        as(protec.areas, "sf")
    
    protec.areas_proj <- 
      st_transform(protec.areas, proj_type)
    
    Links_NatParks <- 
      suppressWarnings(st_intersection(DATA_SF, protec.areas_proj))
    
    Links_NatParks <- st_set_geometry(Links_NatParks, NULL)
    
    # LocNatParks <-
    #   vector(mode = "numeric", length = length(list_data))

    XY_all <-
      data.frame(
        st_coordinates(DATA_SF)[,c(2, 1)],
        tax = as.character(DATA_SF$tax),
        ID_prov_data = DATA_SF$ID_prov_data, 
        stringsAsFactors = F
      )
    
    XY_PA <- 
      XY_all[which(XY_all$ID_prov_data %in% Links_NatParks$ID_prov_data),]
    XY_NOT_PA <- 
      XY_all[which(!XY_all$ID_prov_data %in% Links_NatParks$ID_prov_data),]
    
    locations_pa <- vector(mode = "numeric", length = length(list_data))
    names(locations_pa) <- names(list_data)
    if (nrow(XY_PA) > 0) {
      if (method_protected_area == "no_more_than_one") {
        ## if method is 'no_more_than_one' the number of location is the number of occupied protected areas
        
        count_protec <- 
          table(Links_NatParks$tax, as.vector(Links_NatParks[, colnames(Links_NatParks) == ID_shape_PA]))
        
        loc_pa <- apply(count_protec, 1, function(x) sum(x > 0))
        
        locations_pa[which(names(locations_pa) %in% names(loc_pa))] <- loc_pa
        
        r2_PA <- NA
        
      } else{
        
        list_data_pa <- split(XY_PA, f = XY_PA$tax)
        
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
                                  max = length(list_data_pa),
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
            
            res <- 
              Locations.estimation(
              coordEAC = list_data_pa[[x]], 
              cell_size = Cell_size_locations, 
              export_shp = TRUE, 
              proj_type = proj_type, 
              method = method,
              nbe_rep = nbe_rep
            )
            
            names(res) <- c("nbe_occ", "spatial")
            res
          }
        
        if(parallel) snow::stopCluster(cl)
        if(show_progress) close(pb)
        
        loc_pa <- unlist(output[names(output) == "nbe_occ"])
        r2_PA <- unlist(output[names(output) == "spatial"])
        names(loc_pa) <-
          names(r2_PA) <-
          names(list_data_pa)
        
        locations_pa[which(names(locations_pa) %in% names(loc_pa))] <- loc_pa
        
        # loc_pa <- unlist(output[names(output) == "nbe_occ"])
        # r2_PA <- unlist(output[names(output) == "spatial"])
        # names(loc_pa) <-
        #   names(r2_PA) <-
        #   gsub(pattern = " ",
        #        replacement = "_",
        #        names(list_data_pa))
        # LocNatParks[names(LocNatParks) %in% names(loc_pa)] <-
        #   loc_pa
      }
    } else{
      r2_PA <- NA
    }
    
    names(locations_pa) <-
      gsub(pattern = " ",
           replacement = "_",
           names(locations_pa))
    
    locations_not_pa <- vector(mode = "numeric", length = length(list_data))
    names(locations_not_pa) <- names(list_data)
    if(nrow(XY_NOT_PA) > 0) {
      
      list_data_not_pa <- split(XY_NOT_PA, f = XY_NOT_PA$tax)
      
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
                                max = length(list_data_not_pa),
                                style = 3)
        
        progress <- function(n)
          utils::setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
      }else{opts <- NULL}
      
      output <-
        foreach::foreach(
          x = 1:length(list_data_not_pa),
          .combine = 'c',
          .options.snow = opts
        ) %d% {
          
          if (!parallel & show_progress)
            utils::setTxtProgressBar(pb, x)
          
          res <- 
            Locations.estimation(
              coordEAC = list_data_not_pa[[x]], 
              cell_size = Cell_size_locations, 
              export_shp = TRUE, 
              proj_type = proj_type, 
              method = method,
              nbe_rep = nbe_rep
            )
          
          names(res) <- c("nbe_occ", "spatial")
          res
        }
      
      if(parallel) snow::stopCluster(cl)
      if(show_progress) close(pb)
      
      loc_not_pa <- unlist(output[names(output) == "nbe_occ"])
      r2 <- unlist(output[names(output) == "spatial"])
      names(loc_not_pa) <-
        names(r2) <-
        names(list_data_not_pa)
      
      locations_not_pa[which(names(locations_not_pa) %in% names(loc_not_pa))] <- loc_not_pa
      
    } else{
      r2 <- NA
    }
    
    names(locations_not_pa) <-
      gsub(pattern = " ",
           replacement = "_",
           names(locations_not_pa))
  }
  
  if (!is.null(protec.areas))
    return(list(r2, r2_PA, locations_pa, locations_not_pa))
  if (is.null(protec.areas))
    return(list(r2, Locations))
  
}


