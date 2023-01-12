#' @title Estimate the number of locations.
#'
#' @description 
#' Estimate the number of locations (sensu IUCN) for multiple taxa, 
#' taking into account spatial threats if provided.
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY [data.frame][base::data.frame()], see details
#' @param method string, indicating the method used for estimating the number of locations. Either "fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"
#' @param nbe_rep numeric , indicate the number of raster with random starting position for estimating the number of locations By default, it is 0 but some minimal translation of the raster are still done
#' @param threat_list list or sfc objects POLYGON or MULTIPOLYGON documenting threats. If provided, this will be taken into account for calculating number of location (see Details and `method_polygons`). By default, no shapefile is provided
#' @param names_threat character vector, indicating names of threats, optional
#' @param threat_weight numeric vector, indicating weight given to each threat
#' @param Cell_size_locations numeric, value indicating the grid size in kilometres used for estimating the number of location. By default, equal to 10
#' @param method_polygons string. Used if `threat_list` is provided. See Details.
#'   * `"no_more_than_one"` (the default): each single POLYGON will be considered as a single location
#'   * `"grid"`: a grid of `Cell_size_locations` size will be used to estimate the number of location within polygons
#' @param id_shape string
#' @param Rel_cell_size numeric, if `method_locations="sliding scale"`, `Cell_size_locations` is ignored 
#' and the resolution is given by the maximum distance separating two occurrences multiplied by `Rel_cell_size`. By default, it is 0.05
#' @param parallel logical, whether running in parallel. By default, it is FALSE
#' @param NbeCores integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
#' @param proj_type character string or numeric or object of CRS class, by default is `"cea"`
#' 
#' @details 
#' ## 
#' `XY` as a `dataframe` should have the following structure:
#' 
#' **It is mandatory to respect field positions, but field names do not matter**
#' 
#' \tabular{ccc}{
#'   [,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#'   [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#'   [,3] \tab tax \tab character or factor, taxa names\cr
#' }
#' 
#' 
#' 
#' Locations are estimated by overlaying a grid of a given resolution (see `Cell_size_locations` for
#' specifying the resolution). The number of locations is  the number of
#' occupied locations. Note that the grid position is overlaid in order to
#' minimize the number of locations (several translation of the grid are
#' performed and the one providing the minimum number of occupied cells is
#' provided).
#' 
#' If `threat_list` is provided, 
#' 
#' which means occurrences within polygon documenting threats (if provided) will not be taken into account for estimating the number of locations following the grid system,
#' 
#' If `method` is "fixed_grid" as it is by default, the resolution is fixed and determined 
#' by the argument `Cell_size_locations`.
#' If `method` is "sliding scale", the resolution is defind as 1/x*max.dist where max.dist is the maximum distance between any pairs of occurrences 
#' and x is a defined parameter. 1/x is defined by `Rel_cell_size` argument and is 0.05 by default. 
#' See Rivers M.C. et al. (2010) for more information on the methods.
#' 
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9
#'
#' @return A list with one list for each species containing [[1]]data.frame with the number of locations and 
#' potential issue for each species and [[2]]sf with representing the squared polygons. If threat_list is not null, then [[2]] is a list for each threat
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
                           threat_list = NULL,
                           names_threat = NULL,
                           threat_weight = NULL,
                           Cell_size_locations = 10,
                           method_polygons = c("no_more_than_one"),
                           id_shape = "id_orig",
                           Rel_cell_size = 0.05,
                           parallel = FALSE,
                           NbeCores = 2,
                           show_progress = TRUE,
                           proj_type = "cea") {
  
  proj_type <- proj_crs(proj_type = proj_type, wkt = T)
  
  if (length(method_polygons) > 1)
    stop('Choose only one method_polygons, either "no_more_than_one" or "grid"')
  
  match.arg(method_polygons, c("no_more_than_one", "grid"))
  
  list_data <- 
    coord.check(XY = XY, proj_type = proj_type, cell_size = Cell_size_locations, check_eoo = FALSE)
  
  issue_close_to_anti <- list_data$issue_close_to_anti
  list_data <- list_data$list_data
  
  if (!is.null(threat_list)) {
    
    if (is.null(threat_weight)) stop("Provide a threat_weight vector of same length of threat_list with score comprised of integers 1, 2 or 3")
    
    if (!is.null(threat_weight)) {
      if (length(threat_weight) != length(threat_list))
        stop("threat_weight should be of same length of threat_list")

      if (any(threat_weight < 1))
        stop("threat_weight must be integers > 1")
    }
    
    if (!is.list(threat_list))
      threat_list <- list(threat_list)
    
    which_rast <- unlist(lapply(threat_list, function(x) inherits(x, "SpatRaster")))
    which_sf <- unlist(lapply(threat_list, function(x) inherits(x, "sf")))
    
    # if (!any(which_rast) & !any(which_sf))
    #   stop("threat_list should be either SpatRaster of sf class")
    
    if (all(!which_rast) & all(!which_sf))
      stop("threat_list should be either SpatRaster of sf class")
    
    # if (any(unlist(lapply(threat_list, function(x) class(x$geometry)[2])) != "sfc"))
    #   stop("threat_list should be sfc class objects")
    
    if (any(which_sf))
      threat_list[which(which_sf)] <-
      lapply(threat_list[which(which_sf)],
             function(x)
               st_transform(x, crs = proj_type))
    
    if (any(which_rast))
      threat_list[which(which_rast)] <- lapply(threat_list[which(which_rast)],
             function(x)
               terra::project(x, proj_type))
    
    if (!is.null(names_threat)) {
      
      if (length(names_threat) != length(threat_list))
        stop("names_threat and threat_list must of the same length")
      
      names(threat_list) <- names(threat_weight) <- names_threat
      
    } else {
      
      if (is.null(names(threat_list)))
        names(threat_list) <- paste0("threat_layer.", seq(1, length(threat_list), 1))
      
      names(threat_weight) <- names(threat_list)
    }
    
  }
  
  ## geographical distances for all pairs of occurrences
  
  if (is.null(threat_list)) {
    
    res_df <-
      data.frame(locations =  rep(NA, length(list_data)), 
                 issue_locations = rep(NA, length(list_data)))
    row.names(res_df) <- names(list_data)
    
    if (length(issue_close_to_anti) > 0) {
      
      list_data <- list_data[-issue_close_to_anti]
      
    }
    
    res_list <- .generate_loc(dataset = list_data,
                              method = method,
                              nbe_rep = nbe_rep,
                              Cell_size_locations = Cell_size_locations,
                              Rel_cell_size = Rel_cell_size,
                              parallel = parallel,
                              NbeCores = NbeCores,
                              show_progress = show_progress,
                              proj_type = proj_type)
    
    res_df[which(row.names(res_df) %in% names(res_list$res_df)), 1] <-
      res_list$res_df
    
    shapes_loc <- res_list$shapes
    
  }
  
  if (!is.null(threat_list)) {
    
    ### Taking into account polygon threats if provided
    XY_ID <- 
      data.frame(XY, ID_prov_data = seq(1, nrow(XY), 1))
    
    DATA_SF <- st_as_sf(coord.check(XY = XY_ID, 
                                    listing = F, 
                                    proj_type = proj_type, 
                                    check_eoo = FALSE)$list_data, 
                        coords = c(2, 1), crs = proj_type)
    
    if (any(which_rast)) {
      threat_list[which_rast] <-
        lapply(threat_list[which_rast], function(x)
          st_make_valid(st_as_sf(stars::st_as_stars(x), merge = TRUE)))
      
      threat_list[which_rast] <- 
        lapply(threat_list[which_rast], function(x) 
          cbind(x, id_orig = 1:nrow(x)))
      
      which_sf[which(which_rast)] <- TRUE
      
    }
    
    if (any(which_sf)) {
      intersects_poly <- lapply(lapply(threat_list[which_sf],
                                       function(x)
                                         st_intersects(x = x, y = DATA_SF)),
                                function(y)
                                  length(unlist(y)))
      
      intersects_poly <-
        unlist(lapply(intersects_poly, function(x)
          x > 0))
      
      XY_all <-
        data.frame(
          st_coordinates(DATA_SF)[, c(2, 1)],
          tax = as.character(DATA_SF$tax),
          ID_prov_data = DATA_SF$ID_prov_data,
          stringsAsFactors = F
        )
      
      ### find for each threats spatial data which id intersect with occurrences
      if (any(intersects_poly)) {
        
        threat_list_sf <- threat_list[which_sf][intersects_poly]
        
        threat_list_inter <-
          lapply(threat_list[which_sf], function(x)
            suppressWarnings(st_set_geometry(st_intersection(DATA_SF, x), NULL)))
        
        threat_list_inter <- 
          lapply(threat_list_inter, function(x) cbind(x, combined = paste0(x$tax, x$ID_prov_data)))
        
        nbe_occ_tax <- data.frame(table(DATA_SF$tax))
        colnames(nbe_occ_tax)[2] <- "tot"
        
        .freq_threats <- function(x, y) {
          
          fre_df <- data.frame(table(unique(x[,c("tax", "ID_prov_data")])[,c("tax")]))
          res <- merge(fre_df, y)
          freq <- res$Freq/res$tot
          names(freq) <- res$Var1
          
          freq_score <- freq
          
          freq_score[freq < 0.5] <- 1
          freq_score[freq >= 0.5] <- 2
          freq_score[freq >= 0.9] <- 3
          
          return(freq_score)
        }
        
        scores_threats <- lapply(threat_list_inter, 
                                 function(x) .freq_threats(x = x, y = nbe_occ_tax))
        
        for (i in 1:length(scores_threats)) scores_threats[[i]] <- 
          scores_threats[[i]] + threat_weight[which(names(threat_weight) == names(scores_threats)[i])]
        
        tax_list_scor <- vector('list', length(unique(DATA_SF$tax)))
        for (i in 1:length(unique(DATA_SF$tax)))
          tax_list_scor[[i]] <- data.frame(tax = rep(unique(DATA_SF$tax)[i], length(threat_list_inter)),
                                           threat = names(threat_list_inter))
        names(tax_list_scor) <- unique(DATA_SF$tax)
        
        
        for (i in 1:length(tax_list_scor)) {
          
          test <- unlist(scores_threats)[unlist(lapply(scores_threats, function(x) names(x) == names(tax_list_scor)[i]))]
          
          if (length(test) > 0) {
            
            names(test) <- gsub("\\..*", "", names(test))
            
            test <- sort(test, decreasing = T)
            test[1:length(test)] <- 1:length(test)
            
            tax_list_scor[[i]] <-
              merge(
                tax_list_scor[[i]],
                data.frame(threat = names(test), rank = test),
                by = "threat",
                all.x = T
              )
            
          } else {
            
            tax_list_scor[[i]] <-
              merge(tax_list_scor[[i]],
                    data.frame(threat = NA, rank = NA),
                    by = "threat",
                    all.x = T)
            
          }
        }
        
        main_threat <- unlist(lapply(tax_list_scor, function(x) {ee = x$threat[x$rank == 1]; ee[!is.na(ee)]}))
        
        tax_df_scor <- do.call('rbind', tax_list_scor)
        
      }
      

      if (any(intersects_poly)) {
        
        rank_locations <-
          data.frame(rank = 0, ID_prov_data = DATA_SF$ID_prov_data, tax = DATA_SF$tax)
        
        rank_threats <- data.frame(id = 1:length(threat_list_inter), 
                                   threat = names(threat_list_inter))
        
          for (j in 1:length(unique(tax_df_scor$rank[!is.na(tax_df_scor$rank)]))) {
          
            rank_selected <- unique(tax_df_scor$rank[!is.na(tax_df_scor$rank)])[j]
            cor_id_threat <- merge(rank_threats, tax_df_scor[which(tax_df_scor$rank == rank_selected),], 
                  by = "threat")
            
            for (i in 1:length(unique(cor_id_threat$id))) {
              
              threat_id <- unique(cor_id_threat$id)[i]
              
              rank_locations[which(
                rank_locations$rank == 0 &
                  rank_locations$tax %in% cor_id_threat$tax &
                  rank_locations$ID_prov_data %in% threat_list_inter[[threat_id]]$ID_prov_data
              ), "rank"] <- threat_id
              
            }
            # which_tax <- tax_df_scor[which(tax_df_scor$rank == j),]
            # which_tax <- which_tax[which(which_tax$threat == names(threat_list_inter)[i]),]
            # which_tax
            
            # if (length(which_tax) > 0)
            #   rank_locations[which(
            #     rank_locations$rank == 0 &
            #       rank_locations$tax %in% which_tax$tax &
            #       rank_locations$ID_prov_data %in% threat_list_inter[[i]]$ID_prov_data
            #   ), "rank"] <- i
            
            # rank_locations %>% filter(tax == "tax11", rank > 0)
            
            # which_tax <- tax_df_scor[which(tax_df_scor$threat == names(threat_list_inter)[i]),]
            # which_tax <- which_tax[which(!is.na(which_tax$rank)),]
            # which_tax <- which_tax[which(!is.na(which_tax$rank)),]
           
          }
        
        unique_ranks <- sort(unique(rank_locations$rank))
        names(unique_ranks)[unique_ranks == 0] <- "not_threatened"
        names(unique_ranks)[unique_ranks != 0] <-
          names(threat_list_inter)[unique_ranks[unique_ranks != 0]]
        
      } else {
        
        unique_ranks <- setNames(0, "not_threatened")
        
      }
    }
    
    res_df <-
      data.frame(matrix(NA,  
                        nrow = length(list_data),
                        ncol = 3 + length(unique_ranks)))
    colnames(res_df) <- c("tax",
      "locations", 
                          "issue_locations", 
                          names(unique_ranks))
    res_df[,c(1, 4:(ncol(res_df)))] <- rep(0, nrow(res_df))
    res_df$tax <- names(list_data)
    
    if (any(intersects_poly))
      res_df <-
      merge(res_df,
            data.frame(main_threat, tax = names(main_threat)),
            by = "tax",
            all.x = T)
    
    if (length(issue_close_to_anti) > 0) {
      
      list_data <- list_data[-issue_close_to_anti]
      
    }
    
    locations_shp <- vector('list', length(unique_ranks))
    shapes_loc <- vector('list', length(unique_ranks))
    # row.names(locations_pa) <- names(list_data)
    for (i in 1:length(unique_ranks)) {
      
      # XY_shp_s <- XY_shp[[i]]
      
      if (method_polygons == "no_more_than_one" & unique_ranks[i] > 0) {
        
        selected_ranked_occ <- rank_locations[rank_locations$rank == unique_ranks[i],]
        
        threat_list_inter_selected <- threat_list_inter[[unique_ranks[i]]]
        threat_list_inter_selected <- threat_list_inter_selected[which(threat_list_inter_selected$ID_prov_data %in% selected_ranked_occ$ID_prov_data),]
        
        ## get rid of duplicated occurrences when polygons overlap
        threat_list_inter_selected_ <- threat_list_inter_selected[!duplicated(threat_list_inter_selected$combined),]

        count_protec <- 
          table(threat_list_inter_selected_$tax, 
                as.vector(threat_list_inter_selected_[, colnames(threat_list_inter_selected_) == id_shape]))
        
        # count_protec <- 
        #   table(threat_list_inter[[unique_ranks[i]]]$tax, 
        #         as.vector(threat_list_inter[[unique_ranks[i]]][, colnames(threat_list_inter[[unique_ranks[i]]]) == id_shape]))
        
        locations_shp[[i]] <- 
          apply(count_protec, 1, function(x) sum(x > 0))
        
        # locations_pa[which(row.names(locations_pa) %in% names(loc_pa)),1] <- loc_pa
        
        ### selecting polygons intersecting with occurrences
        occ_poly <-
          threat_list[[unique_ranks[i]]][which(threat_list[[unique_ranks[i]]]$id_orig %in% unique(threat_list_inter_selected$id_orig)), ][, c("id_orig")]
        
        all_un <- unique(threat_list_inter_selected[,c("tax", "id_orig")])
        all_un <- data.frame(all_un, threat = names(threat_list_inter)[unique_ranks[i]])
        
        occ_poly <- st_sf(merge(x = all_un, occ_poly, by = "id_orig"))
        
        # occ_poly <-
        #   cbind(occ_poly,
        #         threat = names(threat_list_inter)[unique_ranks[i]],
        #         tax = threat_list_inter[[unique_ranks[i]]]$tax)
        
        occ_poly <- st_transform(occ_poly, crs = 4326)
        
        shapes_loc[[i]] <- occ_poly
          
        
      } else {
        
        ### selecting occurrences not threatened
        XY_shp <- 
          XY_all[which(XY_all$ID_prov_data %in% 
                         rank_locations[which(rank_locations$rank == unique_ranks[i]), "ID_prov_data"]),]
        
        # XY_shp <- lapply(threat_list_inter,
        #                  function(x) XY_all[which(XY_all$ID_prov_data %in% x$ID_prov_data),])
        
        list_data_pa <- split(XY_shp, f = XY_shp$tax)
        
        res_list <- .generate_loc(dataset = list_data_pa,
                                  method = method,
                                  nbe_rep = nbe_rep,
                                  Cell_size_locations = Cell_size_locations,
                                  Rel_cell_size = Rel_cell_size,
                                  parallel = parallel,
                                  NbeCores = NbeCores,
                                  show_progress = show_progress,
                                  proj_type = proj_type)
        
        locations_shp[[i]] <- res_list$res_df
        
        occ_poly <-
          cbind(res_list$shapes,
                threat = names(unique_ranks[i]),
                id_orig = NA)
        
        shapes_loc[[i]]  <- occ_poly
        
      }
    }
    
    test <- do.call('rbind', lapply(locations_shp, function(x) data.frame(nbe = x, taxa = names(x))))
    
    summed_locations <- unlist(by(test[,"nbe"], test[,"taxa"], sum, simplify = F))
    
    res_df[which(res_df$tax %in% names(summed_locations)), "locations"] <-
      summed_locations
    
    for (i in 1:length(locations_shp))
      res_df[which(res_df$tax %in% names(locations_shp[[i]])), i + 3] <-
      locations_shp[[i]]
    
    # names(shapes_loc) <- names(unique_ranks)
    
    shapes_loc <- do.call('rbind', shapes_loc)
    
    if (method_polygons != "no_more_than_one")
      shapes_loc <- shapes_loc[ ,-which(colnames(shapes_loc) == "id_orig")]
    
  }
  
  if (length(issue_close_to_anti) > 0)
    res_df[issue_close_to_anti, "issue_locations"] <-
    "AOO could not computed because grid cells would overlap with antimeridian"
  
  # if (!is.null(protec.areas))
  #   return(list(locations_pa = locations_pa,
  #               locations_not_pa = locations_not_pa,
  #               locations_poly_pa = r2_PA,
  #               locations_poly_not_pa = r2))
  
  # if (is.null(protec.areas))
  return(list(locations = res_df,
              locations_poly = shapes_loc))
  
}




#' @title Internal function.
#' 
#' @param dataset list
#' 
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
.generate_loc <- function(dataset,
                          method = "fixed_grid",
                          nbe_rep = 0,
                          Cell_size_locations = 10,
                          Rel_cell_size = 0.05,
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          proj_type = "cea") {
  
  if (length(method) > 1)
    stop('Choose only one method_polygons, either "fixed_grid" or "sliding scale"')
  
  match.arg(method, c("fixed_grid", "sliding scale"))
  
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
                            max = length(dataset),
                            style = 3)
    
    progress <- function(n)
      utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }
  
  output <-
    foreach::foreach(
      x = 1:length(dataset),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        utils::setTxtProgressBar(pb, x)
      
      res <- 
        Locations.estimation(
          coordEAC = dataset[[x]], 
          cell_size = Cell_size_locations, 
          export_shp = TRUE, 
          proj_type = proj_type, 
          method = method,
          nbe_rep = nbe_rep,
          Rel_cell_size = Rel_cell_size
        )
      
      names(res) <- c("nbe_occ", "spatial")
      names(res)[1] <- dataset[[x]]$tax[1]
      res$spatial <- cbind(res$spatial[, "geometry"], tax = dataset[[x]]$tax[1])
      
      res
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  # Locations <- unlist(output[names(output) != "spatial"])
  
  res_df <-
    unlist(output[names(output) != "spatial"])
  
  shapes <- output[names(output) == "spatial"]
  shapes <- do.call('rbind', shapes)
  row.names(shapes) <- 1:nrow(shapes)
  
  # r2 <- output[names(output) == "spatial"]
  # names(Locations) <-
  #   names(r2) <-
  #   gsub(pattern = " ",
  #        replacement = "_",
  #        names(list_data))
  
  return(list(res_df = res_df,
               shapes = shapes))
}

