#' @title Assessment of fragmentation subpopulations
#'
#' @description 
#'  `r lifecycle::badge("experimental")`
#'  Assess whether species populations are severely fragmented
#'
#' @param XY data frame
#' @param cell_size numeric
#' @param nbe_rep numeric, by default is 0. Indicate the number of
#'   raster with random starting position used for estimating the AOO. If 0 but
#'   some translation of the raster are still done.
#' @param AOO_poly numeric. AOO of the species in squared kilometres.
#' @param resol_sub_pop numeric. Defines in kilometres the radius of the circles
#'   around each occurrence
#' @param habitat_poly Simple feature collection documenting habitat of species, see Details
#' @param dist_isolated numeric vector providing a value for each species. Distance in kilometres for identifying
#'   subpopulations considered to be isolated. See details
#' @inheritParams proj_crs
#' @param export_shp logical
#' @inheritParams activate_parallel
#' @param show_progress logical. Whether progress informations should displayed. TRUE by default
#' @param threshold_severe numeric of one value indicates the threshold to identify severe fragmentation. By default is 50
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
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-48.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-37.2),
#'                    tax = rep("a", 7),
#'                    stringsAsFactors = FALSE)
#' 
#' severe_frag(XY = mydf, dist_isolated = 200, resol_sub_pop = 10)
#' 
#' 
#' 
#' 
#' @importFrom sf st_intersects
#' @export
severe_frag <- function(XY = NULL, 
                        cell_size = 2,
                        nbe_rep = 0,
                        AOO_poly = NULL,
                        resol_sub_pop = NULL,
                        habitat_poly = NULL,
                        dist_isolated = NULL,
                        proj_type = "cea",
                        export_shp = FALSE,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2,
                        threshold_severe = 50) {
  
  if (is.null(AOO_poly) & is.null(cell_size))
    stop("Please provide the 'cell_size' (in km) of the grid cells used for AOO")
  
  if (is.null(habitat_poly) & is.null(resol_sub_pop))
    stop("'habitat_poly' is not provided so you must provide 'resol_sub_pop'")
  
  if (!is.null(XY)) {
    if (length(dist_isolated) != length(unique(XY[,3]))) {
      if (length(dist_isolated) == 1) {
        dist_isolated <- rep(dist_isolated, length(unique(XY[,3])))
      } else {
        stop("'dist_isolated' should be a vector with as many as values as there are species in the provided dataset or with one value")
      }
    }
  }
  
  if (is.null(AOO_poly)) {
    
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
    
    AOO_poly <-
      res_aoo$AOO_poly
    
  }

  if (is.null(habitat_poly)) {
    
    if (show_progress) message("Subpopulations computation")
    res_subpop <-
      subpop.comp(
        XY = XY, 
        resol_sub_pop = resol_sub_pop,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type, 
        export_shp = TRUE
      )
    
    habitat_poly <-
      res_subpop$poly_subpop
    
  }
  
  ### This provide the number of isolated patches for each species
  if (show_progress) message("Isolated patches computation")
  
  all_tax <- unique(habitat_poly$tax)
  if (identical(class(dist_isolated), "data.frame")) {
    all_tax_ok <- dist_isolated[which(!is.na(dist_isolated[,2])),1]
    dist_isolated <- dist_isolated[which(!is.na(dist_isolated[,2])),2]
    habitat_poly <- habitat_poly[which(habitat_poly$tax %in% all_tax_ok),]
    AOO_poly <- AOO_poly[which(AOO_poly$tax %in% all_tax_ok),]
  }
  
  res_isol <- 
    get.patches(subpop_poly = habitat_poly, 
                dist_isolated = dist_isolated, export_shp = TRUE, show_progress = show_progress)
  
  list_frag <- split(res_isol$subpop_poly, f = res_isol$subpop_poly$tax)
  list_aoo <- split(AOO_poly, f = AOO_poly$tax)
  
  cl <- activate_parallel(parallel = parallel, NbeCores = NbeCores)
  `%d%` <- c_par(parallel = parallel)
  pro_res <- display_progress_bar(show_progress = show_progress, max_pb = length(list_frag))
  opts <- pro_res$opts
  pb <- pro_res$pb
  
  output <-
    foreach::foreach(
      x = 1:length(list_aoo),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      
      aoo <- list_aoo[[x]]
      isol <- list_frag[[x]]
      isol <- isol[which(isol$frag == "isolated"),]
      
      if (nrow(isol) > 0) {
        
        intersects_aoo <- sf::st_intersects(isol, aoo, sparse = TRUE)
        frac <- nrow(aoo[unique(unlist(intersects_aoo)),])/nrow(aoo)*100
        
      } else {
        
        frac <- 0
        
      }
      
      frac
      
    }
  
  res <- data.frame(tax = names(list_aoo), frac = output, severe_frag = output > threshold_severe)
  
  severe_freg_results <- merge(data.frame(tax = all_tax), res, by = c("tax"), all.x = T)
  
  if (export_shp) {
    return(list(res = severe_freg_results, isolated_patches = res_isol$subpop_poly))
  } else {
    return(severe_freg_results)
  }
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
#' @param export_shp logical
#' @inheritParams activate_parallel
#' @inheritParams proj_crs
#' 
#' 
#' @keywords internal
#' 
#' @importFrom sf st_as_sf st_transform st_coordinates st_distance st_geometry
#' 
get.patches <- function(
                        subpop_poly = NULL,
                        dist_isolated = NULL,
                        proj_type = "cea",
                        export_shp = FALSE,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2) {
  
  
  res_subpop_poly_proj <- st_transform(subpop_poly, crs = proj_crs(proj_type = proj_type))
  
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
      
      dist_btw_poly <- sf::st_distance(subpop) 
      
      dist_btw_poly <- matrix(dist_btw_poly,
                              nrow = nrow(dist_btw_poly), ncol = nrow(dist_btw_poly))
      dist_btw_poly[upper.tri(dist_btw_poly, diag = TRUE)] <- NA
      
      mat_above_thres <- dist_btw_poly > dist_isolated[x] * 1000
      
      isolated_subpop1 <- 
        apply(mat_above_thres, 1, FUN = function(x) all(x, na.rm = T))
      isolated_subpop2 <- 
        apply(mat_above_thres, 2, FUN = function(x) all(x, na.rm = T))
      
      isolated_subpop_poly <-
        subpop[isolated_subpop1 & isolated_subpop2, ]
      
      if (nrow(isolated_subpop_poly) > 0) {
        df <- data.frame(tax = subpop$tax[1],
                         frag = "isolated", stringsAsFactors = FALSE)
        isolated_subpop_poly <-
          sf::st_as_sf(data.frame(sf::st_geometry(isolated_subpop_poly), df))
      }
      
      connected_subpop_poly <- subpop[!(isolated_subpop1 & isolated_subpop2),]
      
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
    
    polys <- st_transform(polys, crs = 4326)
    
    return(list(isolated_subpop_poly = isolated_subpop_poly, subpop_poly = polys))
    
  } else {
    return(isolated_subpop_poly)
  }
  
}