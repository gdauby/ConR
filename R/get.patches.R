


#' @title Assessment of fragmentation subpopulations
#'
#' @description 
#'  `r lifecycle::badge("experimental")`
#'  Assess whether species populations are severaly fragmented
#'
#' @param XY data frame
#' @param cell_size numeric
#' @param nbe_rep numeric
#' @param AOO numeric. AOO of the species in squared kilometres.
#' @param Resol_sub_pop numeric. Defines in kilometres the radius of the circles
#'   around each occurrence
#' @param subpop_poly Simple feature collection, output of
#'   ```subpop.estimation``` function
#' @param dist_isolated numeric vector providing a value for each species. Distance in kilometres for identifying
#'   subpopulations considered to be isolated
#' @param proj_type projected coordinate system (in meters)
#' @param export_shp logical
#' 
#' @details This function evaluates the proportion of the total area of
#' occupancy in habitat patches separated from others by a large distance. Based
#' on IUCN guidelines page(48): "A taxon can be considered to be severely
#' fragmented if most (>50%) of its total area of occupancy is in habitat
#' patches that are (...) (2) separated from other habitat patches by a large
#' distance."
#'
#' This function interpret subpopulations as obtained using
#' ```subpop.estimation``` function as habitat patches. First, subpopulations
#' that are isolated i.e. distant to all others subpopulations by at least
#' ```dist_isolated``` kilometres, are identified. Then, the percentage of the
#' AOO concerned by those "isolated" subpopulations is calculated.
#'   
#' @examples
#' 
#' \donttest{
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-48.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-37.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#' 
#' get.patches(XY = mydf, dist_isolated = 200)
#' }
#' 
#' 
#' 
#' @importFrom sf st_as_sf st_transform st_coordinates st_distance st_geometry
#' @export
frag_assess <- function(XY, 
                        cell_size = 2,
                        nbe_rep = 0,
                        AOO = NULL,
                        Resol_sub_pop = NULL,
                        subpop_poly = NULL,
                        dist_isolated = NULL,
                        proj_type = "cea",
                        export_shp = FALSE,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2) {
  
  
  if (!is.null(AOO) & is.null(cell_size))
    stop("Please provide the size (in km) of the grid cells used for AOO")
  
  
  if (length(dist_isolated) != length(unique(XY[,3])))
    stop("dist_isolated should be a vector with as many as values as there are species in the provided dataset")
  
  # proj_type <- 
  #   proj_crs(proj_type = proj_type)
  # 
  # list_data <- coord.check(XY = XY, proj_type = proj_type)
  
  ## if no AOO is provided, getting the AOO raster for the points
  if (is.null(AOO)) {
    
    if (show_progress) message("Area of occupancy computation")
    res_aoo <-
      AOO.computing(
        XY = XY,
        cell_size_AOO = cell_size,
        nbe.rep.rast.AOO = nbe_rep,
        export_shp = TRUE,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type
      )
    
    # res_aoo <-
    #   #For Casearia sylvestris (AOO = 3104): 24.5 secs; For Annona neosalicifoli (AOO = 604): 1.11 secs
    #   AOO.estimation(
    #     coordEAC = XY_proj_coord[, c(2, 1)],
    #     cell_size = cell_size,
    #     nbe_rep = nbe_rep,
    #     export_shp = TRUE,
    #     proj_type = proj_type
    #   )
    
    res_aoo_poly <-
      res_aoo$AOO_poly
    
    cells <- table(st_drop_geometry(res_aoo_poly)) ### RENATO, it does not make sense to use the number of row in the polygon
    
    # cells <- nrow(res_aoo_poly)
    
  }
  # else {
  #   cells <- round(AOO / (cell_size * cell_size), 0)  ### RENATO, It do not undderstand this neither
  # }  
  ## if no subpopulation sf provided, getting it
  if (is.null(subpop_poly)) {
    
    if (show_progress) message("Subpopulations computation")
    res_subpop <-
      subpop.comp(
        XY = XY, 
        Resol_sub_pop = Resol_sub_pop,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type, 
        export_shp = TRUE
      )
    
    # res_subpop <-          #For Casearia sylvestris: 0.30 secs; For Annona neosalicifolia: 0.07 secs
    #   subpop.estimation(
    #     XY = XY_proj_coord,
    #     Resol_sub_pop = Resol_sub_pop,
    #     proj_type = proj_type,
    #     export_shp = TRUE
    #   )
    
    res_subpop_poly <-
      res_subpop$poly_subpop
    
  } else {
    
    res_subpop_poly <- 
      subpop_poly
    
  }
  
  ### This provide the number of isolated patches for each species
  if (show_progress) message("Isolated patches computation")
  res <- get_isolated_subpop <- get.patches(subpop_poly = res_subpop_poly, dist_isolated = dist_isolated)
  
  return(res)
}






#' @title Internal function
#'
#' @description 
#'  `r lifecycle::badge("experimental")`
#'  Get Occupied Cells and Patches
#'
#' @param subpop_poly Simple feature collection, output of
#'   ```subpop.estimation``` function
#' @param dist_isolated numeric. Distance in kilometres for identifying
#'   subpopulations considered to be isolated
#' @param proj_type projected coordinate system (in meters)
#' @param export_shp logical
#' 
#' @details This function evaluates the proportion of the total area of
#' occupancy in habitat patches separated from others by a large distance. Based
#' on IUCN guidelines page(48): "A taxon can be considered to be severely
#' fragmented if most (>50%) of its total area of occupancy is in habitat
#' patches that are (...) (2) separated from other habitat patches by a large
#' distance."
#'
#' This function interpret subpopulations as obtained using
#' ```subpop.estimation``` function as habitat patches. First, subpopulations
#' that are isolated i.e. distant to all others subpopulations by at least
#' ```dist_isolated``` kilometres, are identified. Then, the percentage of the
#' AOO concerned by those "isolated" subpopulations is calculated.
#'   
#' @examples
#' 
#' \donttest{
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-48.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-37.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#' 
#' get.patches(XY = mydf, dist_isolated = 200)
#' }
#' 
#' 
#' @keywords internal
#' 
#' @importFrom sf st_as_sf st_transform st_coordinates st_distance st_geometry
#' 
get.patches <- function(
    # XY, 
                        # aoo_poly = NULL,
                        subpop_poly = NULL,
                        dist_isolated = NULL,
                        proj_type = "cea",
                        export_shp = FALSE,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2) {
  
  # proj_type <- 
  #   proj_crs(proj_type = proj_type)
  # 
  # ## Creating and projecting the points  
  # XY_proj <-
  #   sf::st_as_sf(XY, coords = c("ddlon", "ddlat"))
  # sf::st_crs(XY_proj) <- 4326   
  # XY_proj <- sf::st_transform(XY_proj, crs = proj_type)
  # XY_proj_coord <- sf::st_coordinates(XY_proj)
  # XY_proj_coord <- as.data.frame(XY_proj_coord)
  
  res_subpop_poly_proj <- st_transform(res_subpop_poly, crs = proj_crs(proj_type = proj_type))
  
  list_dataset <- split(res_subpop_poly_proj, f = res_subpop_poly_proj$tax)
  
  cl <- activate_parallel(parallel = parallel, NbeCores = NbeCores)
  `%d%` <- c_par(parallel = parallel)
  pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_dataset))
  opts <- pro_res$opts
  pb <- pro_res$pb
  
  output <-
    foreach::foreach(
      x = 1:length(list_dataset),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      
      subpop <- list_dataset[[x]]
      
      dist_btw_poly <- sf::st_distance(subpop) #For Casearia sylvestris (11 subpops): 0.85 secs; For Annona neosalicifolia (12 subpops): 0.45 secs 
      
      dist_btw_poly <- matrix(dist_btw_poly,
                              nrow = nrow(dist_btw_poly), ncol = nrow(dist_btw_poly))
      dist_btw_poly[upper.tri(dist_btw_poly, diag = TRUE)] <- NA
      
      mat_above_thres <- dist_btw_poly > dist_isolated[x] * 1000
      
      isolated_subpop1 <- 
        apply(mat_above_thres, 1, FUN = function(x) all(x, na.rm = T))
      isolated_subpop2 <- 
        apply(mat_above_thres, 2, FUN = function(x) all(x, na.rm = T))
      
      isolated_subpop_poly <-
        res_subpop_poly[isolated_subpop1 & isolated_subpop2, ]
      
      if (nrow(isolated_subpop_poly) > 0) {
        df <- data.frame(tax = subpop$tax[1],
                         frag = "isolated", stringsAsFactors = FALSE)
        isolated_subpop_poly <-
          sf::st_as_sf(data.frame(sf::st_geometry(isolated_subpop_poly), df))
      }
      
      connected_subpop_poly <- res_subpop_poly[!(isolated_subpop1 & isolated_subpop2),]
      
      if (nrow(connected_subpop_poly) > 0) {
        df <- data.frame(tax = subpop$tax[1],
                         frag = "connected")
        connected_subpop_poly <- 
          sf::st_as_sf(data.frame(sf::st_geometry(connected_subpop_poly), df))
      }
      
      # fraction_aoo <- 100 * nrow(isolated_subpop_poly)/cells
      
      list(poly = rbind(isolated_subpop_poly, connected_subpop_poly),
           nbe_isolated_poly = data.frame(tax = subpop$tax[1],
                                          isolated_subpop_poly = nrow(isolated_subpop_poly)))
      
    }
  
  isolated_subpop_poly <- output[names(output) == "nbe_isolated_poly"]
  isolated_subpop_poly <- do.call('rbind', isolated_subpop_poly)
  row.names(isolated_subpop_poly) <- 1:nrow(isolated_subpop_poly)
  
  if (export_shp) {
    
    polys <- output[names(output) == "poly"]
    polys <- do.call('rbind', polys)
    row.names(polys) <- 1:nrow(polys)
    
    return(list(isolated_subpop_poly = isolated_subpop_poly, subpop_poly = polys))
    
  } else {
    # names(fraction_aoo) <- as.character(unique(XY$tax))
    return(isolated_subpop_poly)
  }
  
  

}