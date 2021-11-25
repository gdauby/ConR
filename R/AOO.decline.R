
#' @title Area of occupancy decline
#'
#' @description Estimate areas of occupancy (AOO) decline for multiple taxa in square kilometers and percentage
#' `r lifecycle::badge("experimental")
#' 
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY \code{"dataframe"} see Details
#' @param Cell_size_AOO numeric, value indicating the grid size in kilometers used for estimating Area of Occupancy.  By default, equal to 2 km (i.e. 4 km2 grid cells)
#' @param nbe.rep.rast.AOO numeric , indicate the number of raster with random starting position for estimating the AOO. By default, it is 0 but some minimal translation of the raster are still done
#' @param parallel logical, whether running in parallel. By default, it is FALSE
#' @param NbeCores string integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
#' @param proj_type character string or numeric or object of CRS class, by default is "cea"
#' @param hab.map raster, raster layer/stack or spatial polygons containing the
#'   habitat spatial information
#' @param hab.map.type logical, vector of same lenght of hab.map, if TRUE means hab.map is suitable for species, if FALSE unsuitable
#' @param hab.class classes of values in ```hab.map``` to be considered as suitable
#' @param all_individual_layers logical, by default FALSE, if TRUE compute AOO decline for each individual hab.map and threat.map provided
#'   
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
#' \enumerate{
#'   \item AOOs a data.frame of AOO estimates for each taxa and for each layer/spatial polygons if all_individual_layers is TRUE
#'   \item AOO_decline a data.frame of AOO.decline in percentages
#'   \item categories categories based on the sub-criteria of IUCN criterion A
#' }
#'   
#' @examples 
#' \dontrun{
#' data(land_cover_ex)
#' AOOs <- AOO.decline(XY = dataset.ex, hab.map = land_cover_ex, hab.class = c(40, 40))
#'}
#'
#' @importFrom terra project crop app vect extract
#' @import sf
#' 
#' @export AOO.decline
AOO.decline <- function(XY,
                        Cell_size_AOO = 2,
                        nbe.rep.rast.AOO = 0,
                        parallel = FALSE,
                        NbeCores = 2,
                        show_progress = TRUE,
                        proj_type = "cea",
                        hab.map = NULL,
                        hab.class = NULL,
                        hab.map.type = NULL,
                        all_individual_layers = FALSE
) {
  
  ### checking arguments validity
  if (is.null(hab.map)) 
    stop("Please provide hab.map or threat.map (a sf polygon object or a raster)")
  
  hab.map.checked <- check_hab_map(hab.map = hab.map, hab.map.type = hab.map.type)
  
  # if (!is.null(hab.map) & !is.null(hab.class)) {
  #   
  #   if (any(grepl("Raster", class(hab.map)))) {
  #     
  #     if (dim(hab.map)[3] != length(hab.class)) {
  #       
  #       stop ("The number of layers of hab.map is different to the number hab.class values. It should be identical")
  #       
  #     }
  #     
  #     if (any(hab.map@data@isfactor) & is.factor(hab.class)) {
  #       
  #       stop ("hab.map values are factors while hab.class are not. They must be of same types")
  #       
  #     }
  #   }
  # }
  
  ### converting XY inputs data
  proj_type_ <- proj_crs(proj_type = proj_type)
  
  XY_id <- data.frame(ddlat = XY[,1],
                      ddlon = XY[,2],
                      species = XY[,3],
                      id = 1:nrow(XY))
  
  XY_sf <-
    st_as_sf(XY_id, coords = c("ddlon", "ddlat"), 
             crs = 4326)
  
  XY_sf <- 
    st_transform(XY_sf, proj_type_)
  
  ids_to_keep <-
    vector('list', sum(hab.map.checked$classes.hab.map$nbe.layers))
  
  for (i in 1:nrow(hab.map.checked$classes.hab.map)) {
    
    if (hab.map.checked$classes.hab.map[i,]$sf) {
      
      hab.map_proj <- 
        st_transform(hab.map.checked$hab.map[[i]], proj_type_)
      
      intersect_treath_map <- suppressWarnings(st_intersection(XY_sf, hab.map_proj))
      
      if (hab.map.checked$hab.map.type[i]) {
        
        ids_to_keep[[i]] <- 
          XY_id$id[XY_id$id %in% st_set_geometry(intersect_treath_map, NULL)[,"id"]]
        
      } else {
        
        ids_to_keep[[i]] <- 
          XY_id$id[!XY_id$id %in% st_set_geometry(intersect_treath_map, NULL)[,"id"]]
        
      }
      
      names(ids_to_keep)[[i]] <- 
        paste0("hab.map.", i)
      
    }
    
    if (hab.map.checked$classes.hab.map[i,]$rast) {
      
      hab.map.selected <- hab.map.checked$hab.map[[i]]
      
      for (j in 1:hab.map.checked$classes.hab.map[i,"nbe.layers"]) {
        
        hab.map.selected_lay <- hab.map.selected[[j]]
        
        if (is.character(hab.class) | is.factor(hab.class)) {
          
          if (!is.factor(hab.map.selected_lay)) {
            
            val_rast <- terra::values(hab.map.selected_lay)[,1]
            val_rast <- unique(val_rast[which(!is.na(val_rast))])
            
            levels(hab.map.selected_lay) <- as.character(val_rast)
            
          }
            
          
          if (hab.map.checked$hab.map.type[i]) {
            ## if suitable
            hab.map.selected_lay <- terra::app(hab.map.selected_lay,
                                           fun = function(x) {
                                             x[which(!x %in% hab.class)] <- NA
                                             return(x)
                                           })
            
          } else {
            ## if unsuitable
            hab.map.selected_lay <- terra::app(hab.map.selected_lay,
                                           fun = function(x) {
                                             x[which(x %in% hab.class)] <- NA
                                             return(x)
                                           })
            
          }
          
        } else {
          
          stop("implement for hab.class numeric")
          
        }
        
        hab.map.selected_lay <-
          terra::crop(hab.map.selected_lay, terra::vect(XY_sf))
        
        if (as.character(crs(hab.map.selected_lay)) != proj_type_@projargs) {
          
          hab.map.selected_lay <- terra::project(hab.map.selected_lay, as.character(proj_type_))
          
        }
        
        intersect_rast_points <-
          terra::extract(hab.map.selected_lay, terra::vect(XY_sf))
        
        index_ <- which(unlist(lapply(ids_to_keep, is.null)))[1]
        
        ### ids stored are for which points do not fall within defined layer (i.e. under threat)
        ids_to_keep[[as.numeric(index_)]] <-
          which(!is.na(intersect_rast_points[, 2]))
        
        names(ids_to_keep)[[as.numeric(index_)]] <-
          paste0("hab.map.", i + j -1)
        
      }
      
    }
    
  }
  
  # if (!is.null(threat.map)) {
  #   
  #   for (i in 1:length(threat.map)) {
  #     
  #     threat.map_proj <- 
  #       st_transform(threat.map[[i]], proj_type_)
  #     
  #     intersect_treath_map <- suppressWarnings(st_intersection(XY_sf, threat.map_proj))
  #     
  #     ids_to_keep[[i]] <- 
  #       XY_id$id[!XY_id$id %in% st_set_geometry(intersect_treath_map, NULL)[,"id"]]
  #     
  #     names(ids_to_keep)[[i]] <- 
  #       paste0("threat.map.", i)
  #     
  #   }
  #   
  # }
  
  # if (!is.null(hab.map)) {
  #   
  #   if (any(grepl("Raster", class(hab.map)))) {
  #     
  #     if (st_crs(crs(hab.map))$input != proj_type_@projargs) {
  #       
  #       hab.map <- raster::projectRaster(hab.map, crs = proj_type_)
  #       
  #     }
  #     
  #     ## masking input raster given hab.class
  #     if (!is.null(hab.class)) {
  #       
  #       hab.map_masked <- vector('list', dim(hab.map)[3])
  #       for (i in 1:dim(hab.map)[3]) {
  #         
  #         hab.map_masked[[i]] <- 
  #           raster::calc(hab.map[[i]], function(x) ifelse(x < hab.class[i], NA, 1))
  #         
  # 
  #         
  #       }
  #       
  #       hab.map_masked <- raster::stack(hab.map_masked)
  #       
  #     } else {
  #       
  #       hab.map_masked <- 
  #         hab.map
  #       
  #     }
  #     
  #     ext <- raster::extent(sf::st_bbox(XY_sf)[c(1, 3, 2, 4)])
  #     hab.map_masked_crop <-
  #     raster::crop(hab.map_masked, ext) # cropping raster to the extent of the polygon
  #     
  #     # }
  #       
  #     intersect_rast_points <- raster::extract(hab.map_masked_crop, XY_sf)
  #       
  #     for (i in 1:ncol(intersect_rast_points)) {
  #       
  #       ids_to_keep[[length(threat.map) + i]] <- 
  #         which(!is.na(intersect_rast_points[,i]))
  #       
  #       names(ids_to_keep)[[length(threat.map) + i]] <- 
  #         paste0("hab.map.", length(threat.map) + i)
  #       
  #     }
  #     
  #   }
  #   
  # }
  
  message("Total Area of occupancy")
  AOO <-
    AOO.computing(
      XY = XY,
      Cell_size_AOO = Cell_size_AOO,
      nbe.rep.rast.AOO = nbe.rep.rast.AOO,
      export_shp = FALSE,
      parallel = parallel,
      show_progress = show_progress,
      NbeCores = NbeCores,
      proj_type = proj_type
    )
  
  all_data <- 
    data.frame(AOO = AOO, species = names(AOO))
  
  ### whether the AOO should be estimated for all layers individually
  if (all_individual_layers & length(ids_to_keep) > 1) {
    
    message("Area of occupancy minus threatened areas for each individual layers")
    
    AOO_res_list <- vector('list', length(ids_to_keep))
    for (i in 1:length(ids_to_keep)) {
      
      print(names(ids_to_keep)[i])
      
      XY_ok <- 
        XY_id[which(XY_id$id %in% ids_to_keep[[i]]),]
      
      AOO_treatened <-
        AOO.computing(
          XY = XY_ok,
          Cell_size_AOO = Cell_size_AOO,
          nbe.rep.rast.AOO = nbe.rep.rast.AOO,
          export_shp = FALSE,
          parallel = parallel,
          show_progress = show_progress,
          NbeCores = NbeCores,
          proj_type = proj_type
        )
      
      AOO_res_list[[i]] <- 
        AOO_treatened
      
    }
    
    for (i in 1:length(AOO_res_list)) {
      AOO_treath <- data.frame(AOO_treatened = AOO_res_list[[i]], 
                               species = names(AOO_res_list[[i]]))
      colnames(AOO_treath)[1] <- names(ids_to_keep)[i]
      all_data <- merge(all_data,
                        AOO_treath,
                        by = "species", all.x = T)
      
      all_data[which(is.na(all_data[,i + 2])), i + 2] <- 0
    }
  }
  
  # XY_ok <-
  #   XY_id[which(!XY_id$id %in% as.numeric(names(which(table(unlist(ids_to_keep)) == length(ids_to_keep))))),]
  
  ## any point that is under at least one threat
  XY_ok <-
    XY_id[which(XY_id$id %in% as.numeric(names(table(unlist(ids_to_keep))))),]
  
  if (nrow(XY_ok) > 0) {
    
    AOO_treatened <-
      AOO.computing(
        XY = XY_ok,
        Cell_size_AOO = Cell_size_AOO,
        nbe.rep.rast.AOO = nbe.rep.rast.AOO,
        export_shp = FALSE,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type
      )
    
    AOO_treath <- data.frame(AOO_treath = AOO_treatened, 
                             species = names(AOO_treatened))
    # colnames(AOO_treath)[1] <- "AOO_treath"
    
  } else {
    
    AOO_treath <- data.frame(AOO_treath = NA, 
                             species = names(AOO))
    
  }
  
  all_data <- merge(all_data,
                    AOO_treath,
                    by = "species", 
                    all.x = T)
  
  all_data[which(is.na(all_data[,ncol(all_data)])), ncol(all_data)] <- 0
  
  all_data_decline <- data.frame(species = all_data[,1], 
                                 aoo_decline = 100 - all_data[,3:ncol(all_data)]/all_data$AOO*100)
  
  if (ncol(all_data_decline) > 3) {
    categories <- apply(all_data_decline[,2:ncol(all_data_decline)], MARGIN = 2, FUN = function(x) cat_criterion_a(A2_val = x))
    categories <- data.frame(species = all_data[,1], as.data.frame(lapply(categories, function(x) as.data.frame(x))))
  } else {
    categories <- data.frame(species = all_data[,1], as.data.frame(cat_criterion_a(A2_val = all_data_decline[,2:ncol(all_data_decline)])))
  }
  
  return(list(
    AOOs = all_data,
    AOO_decline = all_data_decline,
    categories = categories
  ))
  
}



