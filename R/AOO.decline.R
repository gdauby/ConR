
#' @title Area of occupancy decline
#'
#' @description 
#'  `r lifecycle::badge("experimental")`
#'  Estimate areas of occupancy (AOO) decline for multiple taxa in square kilometres and percentage
#'  
#' 
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#'
#' @inheritParams AOO.computing
#' @param hab.map SpatRaster or sf polygons containing the
#'   habitat spatial information
#' @param hab.map.type logical, vector of same length of `hab.map`,
#'  * `TRUE` means the habitat of `hab.map` is suitable
#'  * `FALSE` means the habitat of `hab.map` is unsuitable
#' @param hab.class classes of values in ```hab.map``` to be considered as suitable
#' @param all_individual_layers logical
#'  * `TRUE` compute AOO decline for each individual `hab.map`
#'  * `FALSE` (the default) compute AOO decline for all layers together
#'   
#' @details 
#' **Input** as a [data.frame][base::data.frame()] should have the following structure:
#' 
#' **It is mandatory to respect field positions, but field names do not matter**
#' 
#' | latitude | longitude | species |
#' | --- | --- | --- |
#' | numeric | numeric | character |
#' 
#' The argument `nbe.rep.rast.AOO` should ideally be higher than 20 for increasing 
#' the chance to get the minimal number of occupied cell. Increasing `nbe.rep.rast.AOO` however 
#' also increase the computing time. So this is a trade-off that depends on the importance to 
#' get the minimal AOO and the size of the dataset.
#' 
#' 
#' 
#' 
#' @return 
#' \enumerate{
#'   \item AOOs a `dataframe` of AOO estimates for each taxa and for each layer/spatial polygons if `all_individual_layers` is TRUE
#'   \item AOO_decline a `dataframe` of AOO.decline in percentages
#'   \item categories based on the sub-criteria of IUCN criterion A
#' }
#' 
#' 
#' @examples
#' 
#' m <- matrix(1:25, nrow=5, ncol=5)
#' rm <- terra::rast(m)
#' terra::values(rm) <- sample(c("forest", "cities", "roads"), 25, replace = TRUE)
#' cls <- data.frame(id=1:3, cover=c("forest", "cities", "roads"))
#' levels(rm) <- cls
#' terra::crs(rm) <- "epsg:4326"
#' 
#' 
#' test_data <- dummy_dist(n = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5)
#' 
#' 
#' res <- AOO.decline(
#' XY = test_data,
#' hab.map = rm,
#' hab.class = c("forest"),
#' all_individual_layers = TRUE
#' )
#' 
#' 
#' res <- AOO.decline(
#' XY = test_data,
#' hab.map = rm,
#' hab.class = c("cities", "roads"),
#' all_individual_layers = TRUE, 
#' hab.map.type = FALSE ### this means the provided hab.map is unsuitable
#' )
#' 
#' 
#' @importFrom terra project crop app vect extract
#' @importFrom lifecycle badge
#' @import sf
#' @rawNamespace import(terra, except = c(points, median, na.omit, quantile, head, tail, predict))
#' 
#' @export AOO.decline
AOO.decline <- function(XY,
                        hab.map,
                        cell_size_AOO = 2,
                        nbe.rep.rast.AOO = 0,
                        parallel = FALSE,
                        NbeCores = 2,
                        show_progress = TRUE,
                        proj_type = "cea",
                        hab.class = NULL,
                        hab.map.type = NULL,
                        all_individual_layers = FALSE
) {
  
  if (!exists("hab.map")) 
    stop("Please provide hab.map (a sf polygon object or a SpatRaster)")
  
  hab.map.checked <- check_hab_map(hab.map = hab.map, 
                                   hab.map.type = hab.map.type)
  
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
    
    if (hab.map.checked$classes.hab.map[i,]$SpatRaster) {
      
      hab.map.selected <- hab.map.checked$hab.map[[i]]
      
      for (j in 1:hab.map.checked$classes.hab.map[i,"nbe.layers"]) {
        
        hab.map.selected_lay <- hab.map.selected[[j]]
        
        if (is.character(hab.class) | is.factor(hab.class)) {
          
          if (!is.factor(hab.map.selected_lay)) {
            
            if (hab.map.checked$hab.map.type[i]) {
              
              hab.map.selected_lay <- subst(x = hab.map.selected_lay, from = as.numeric(hab.class), to = 1, others= NA)
              
            } else {
              
              hab.map.selected_lay <- subst(x = hab.map.selected_lay, from = as.numeric(hab.class), to = NA, others= 1)
              
            }
            
            # val_rast <- terra::values(hab.map.selected_lay)[,1]
            # length(val_rast[which(!is.na(val_rast))])
            # ct <- table(val_rast[which(!is.na(val_rast))])
            # ct[names(ct) %in% hab.class]
            # val_rast <- unique(val_rast[which(!is.na(val_rast))])
            # 
            # levels(hab.map.selected_lay) <-
            #   data.frame(id = 1:length(val_rast), cover = as.character(val_rast))

          } else {
            
            levels_rast <- levels(hab.map.selected_lay)[[1]]
            fac_levels <- levels_rast[which(levels_rast[,2] %in% hab.class),][,1]
            
            if (hab.map.checked$hab.map.type[i]) {
              ## if suitable
              
              hab.map.selected_lay <- subst(x = hab.map.selected_lay, from = fac_levels, to = 1, others= NA, raw = TRUE)
              
              # hab.map.selected_lay <- terra::app(hab.map.selected_lay,
              #                                    fun = function(x) {
              #                                      x[which(!x %in% hab.class)] <- NA
              #                                      return(x)
              #                                    })
              
              
            } else {
              ## if unsuitable
              # hab.map.selected_lay <- terra::app(hab.map.selected_lay,
              #                                    fun = function(x) {
              #                                      x[which(x %in% fac_levels)] <- NA
              #                                      return(x)
              #                                    })
              hab.map.selected_lay <- subst(x = hab.map.selected_lay, from = fac_levels, to = NA, others= 1, raw = TRUE)
              
            }
            
          }
        } else {
          
          stop("implement for hab.class numeric")
          
          val_rast <- terra::values(hab.map.selected_lay)[,1]
          val_rast <- unique(val_rast[which(!is.na(val_rast))])
          
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
          
          
        }
        
        if (paste(terra::crs(hab.map.selected_lay, describe = T)[,c("authority", "code")], collapse = ":") !=
             proj_type_$input) {
          
          hab.map.selected_lay <- 
            terra::project(hab.map.selected_lay, proj_type_$wkt)
          
        }
        
        XY_te <- terra::vect(XY_sf)
        
        index_ <- which(unlist(lapply(ids_to_keep, is.null)))[1]
        
        is_intersect <- is.related(hab.map.selected_lay, XY_te, "intersects")
        
        if (is_intersect) {
          
          hab.map.selected_lay <-
            terra::crop(hab.map.selected_lay, XY_te)
          
          intersect_rast_points <-
            terra::extract(hab.map.selected_lay, XY_te)
          
          ### ids stored are for which points do not fall within defined layer (i.e. under threat)
          ids_to_keep[[as.numeric(index_)]] <-
            which(!is.na(intersect_rast_points[, 2]))
          
        }
        
        names(ids_to_keep)[[as.numeric(index_)]] <-
          paste0("hab.map.", i + j -1)
        
      }
    }
  }
  
  message("Total Area of occupancy")
  AOO <-
    AOO.computing(
      XY = XY,
      cell_size_AOO = cell_size_AOO,
      nbe.rep.rast.AOO = nbe.rep.rast.AOO,
      export_shp = FALSE,
      parallel = parallel,
      show_progress = show_progress,
      NbeCores = NbeCores,
      proj_type = proj_type
    )
  
  all_data <-
    data.frame(
      AOO = AOO$aoo,
      AOO_issue = AOO$issue_aoo,
      tax = AOO$tax
    )
  
  ### whether the AOO should be estimated for all layers individually
  if (all_individual_layers & any(unlist(lapply(ids_to_keep, length)) > 1)) {
    
    message("Area of occupancy minus threatened areas for each layers individually")
    
    AOO_res_list <- vector('list', length(ids_to_keep))
    for (i in 1:length(ids_to_keep)) {
      
      print(names(ids_to_keep)[i])
      
      XY_ok <- 
        XY_id[which(XY_id$id %in% ids_to_keep[[i]]),]
      
      AOO_treatened <-
        AOO.computing(
          XY = XY_ok,
          cell_size_AOO = cell_size_AOO,
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
      AOO_treath <- data.frame(AOO_treatened = AOO_res_list[[i]]$aoo, 
                               tax = AOO_res_list[[i]]$tax)
      colnames(AOO_treath)[1] <- names(ids_to_keep)[i]
      all_data <- merge(all_data,
                        AOO_treath,
                        by = "tax", all.x = T)
      
      all_data[which(is.na(all_data[,i + 3])), i + 3] <- 0
    }
  }
  
  XY_ok <-
    XY_id[which(XY_id$id %in% as.numeric(names(table(unlist(ids_to_keep))))),]
  
  if (nrow(XY_ok) > 0) {
    
    AOO_treatened <-
      AOO.computing(
        XY = XY_ok,
        cell_size_AOO = cell_size_AOO,
        nbe.rep.rast.AOO = nbe.rep.rast.AOO,
        export_shp = FALSE,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type
      )
    
    AOO_treath <- data.frame(AOO_treath = AOO_treatened$aoo, 
                             AOO_issue_threat = AOO_treatened$issue_aoo,
                             tax = AOO_treatened$tax)
    # colnames(AOO_treath)[1] <- "AOO_treath"
    
  } else {
    
    AOO_treath <- data.frame(AOO_treath = NA, 
                             tax = AOO$tax)
    
  }
  
  all_data <- merge(all_data,
                    AOO_treath,
                    by = "tax", 
                    all.x = T)
  
  all_data[which(is.na(all_data[,which(grepl("AOO_treath", colnames(all_data)))])), 
           which(grepl("AOO_treath", colnames(all_data)))] <- 0
  
  # print(all_data$tax)
  # print(all_data$AOO)
  # print(all_data)
  # print(all_data[, which(grepl("AOO_threat", colnames(all_data)))])
  
  all_data_decline <-
    data.frame(tax = all_data$tax,
               aoo_decline = round(100 - all_data[, which(grepl("AOO_treath", colnames(all_data)))] /
                                     all_data$AOO * 100, 2))
  
  if (ncol(all_data_decline) > 3) {
    categories <- apply(all_data_decline[,2:ncol(all_data_decline)], 
                        MARGIN = 2, FUN = function(x) cat_criterion_a(A2_val = x))
    categories <- data.frame(tax = all_data$tax, 
                             as.data.frame(lapply(categories, function(x) as.data.frame(x))))
  } else {
    categories <- data.frame(tax = all_data$tax, 
                             as.data.frame(cat_criterion_a(A2_val = all_data_decline[,2:ncol(all_data_decline)])))
  }
  
  return(list(
    AOOs = all_data,
    AOO_decline = all_data_decline,
    categories = categories
  ))
  
}



