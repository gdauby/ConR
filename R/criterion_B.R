#' @title Assess IUCN Criterion B
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion B, which is based on species geographic distribution (i.e. extent
#'  of occurrence - EOO, and area of occupancy, AOO)
#'
#' @param x a `dataframe` or an object of class `spgeoIN` see
#' <https://github.com/azizka/speciesgeocodeR>. See Details
#' @param EOO.threshold numeric vector indicating the thresholds used to categorize EOO in IUCN categories
#' @param AOO.threshold numeric vector indicating the thresholds used to categorize AOO in IUCN categories
#' @param Loc.threshold numeric vector indicating the thresholds used to categorize the number of locations in IUCN categories
#' @param SubPop logical if the number of sub-populations should be computed. By default is TRUE
#' @inheritParams locations.comp
#' @inheritParams EOO.computing
#' @inheritParams AOO.computing
#' 
#' 
#' @return A data frame containing, for each of taxon, (EOO, AOO, n.locs, n.subpops?),
#'   the IUCN categories associated with the sub-criteria and the consensus category
#'   for criterion B.
#' 
#' @details The function ... 
#' 
#' @author Dauby, G. & Lima, R.A.F.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#'
#' @examples To be included
#' 
#' @importFrom tibble is_tibble
#' @importFrom rnaturalearth ne_countries
#' 
#' @export criterion_B
criterion_B <- function(x, 
                       #add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                       EOO.threshold = c(20000, 5000, 100), 
                       AOO.threshold = c(2000, 500, 10), 
                       Loc.threshold = c(10, 5, 1),
                       SubPop = TRUE,
                       Resol_sub_pop = 5,
                       Cell_size_locations = 10,
                       method_locations = "fixed_grid",
                       method_polygons = "no_more_than_one",
                       Rel_cell_size = 0.05,
                       threat_list = NULL,
                       names_threat = NULL,
                       threat_weight = NULL,
                       id_shape = "id_orig",
                       country_map = NULL,
                       method.range = "convex.hull",
                       alpha = 1,
                       buff.alpha = 0.1,
                       exclude.area = FALSE,
                       Cell_size_AOO = 2,
                       nbe.rep.rast.AOO = 0,
                       parallel = FALSE,
                       show_progress = TRUE,
                       NbeCores = 2,
                       proj_type = "cea",
                       mode = "spheroid",
                       DrawMap = FALSE) {
  
  if(class(x)[1] == "spgeoIN") {
    x <- cbind(x$species_coordinates, x$identifier)
    x <- x[,c(2,1,3)]
  }
  colnames(x)[1:3] <- c("ddlat","ddlon","tax")
  
  if(tibble::is_tibble(x)) x <- as.data.frame(x)
  
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  } else {
    
    # if(any(grepl('sf', class(country_map))))
    #   country_map <- 
    #     as(country_map, "Spatial")
    
    if (any(!st_is_valid(country_map)))
      country_map <- sf::st_make_valid(country_map)
    
    # country_map <-
    #   suppressWarnings(sf::st_buffer(country_map, dist = 0))
    
    # country_map <- 
    #   as(country_map, "sf")
  }
  
  # if (!is.null(threat_list)) {
  #   # if (!sp::identicalCRS(protec.areas, country_map)) {
  #   #   
  #   #   sp::proj4string(protec.areas) <- 
  #   #     sp::CRS(SRS_string = 'EPSG:4326')
  #   #   sp::proj4string(country_map) <-
  #   #     sp::CRS(SRS_string = 'EPSG:4326')
  #   #   
  #   # }
  #   
  #   if (st_crs(protec.areas) != st_crs(country_map)) {
  #     
  #     st_crs(protec.areas) <- 4326
  #     st_crs(country_map) <- 4326
  #     
  #   }
  #   
  #   
  # }
  
  
  ##########################################################################################
  ##############################  Sub-populations estimations ############################## 
  
  if(SubPop) {
    
    message("Subpopulations computation")
    
    subpop_stats <-
      subpop.comp(
        XY = x,
        Resol_sub_pop = 10,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type, 
        export_shp = ifelse(DrawMap, TRUE, FALSE)
      )
    
    if (DrawMap) {
      
      SubPopPoly <- subpop_stats$poly_subpop
      NbeSubPop <- subpop_stats$number_subpop
      
    } else {
      
      NbeSubPop <- subpop_stats
      
    }
    
  }
  
  ##########################################################################################
  ##############################  Estimations of number of Locations ####################### 
  
  message("Locations computation")
  
  locations_res <-
    locations.comp(
      XY = x,
      method = method_locations, 
      threat_list = threat_list, 
      threat_weight = threat_weight,
      names_threat = names_threat,
      Cell_size_locations = Cell_size_locations, 
      method_polygons = method_polygons, 
      id_shape = id_shape,
      show_progress = show_progress,
      Rel_cell_size = Rel_cell_size,
      parallel = parallel,
      NbeCores = NbeCores, 
      proj_type = proj_type
    )
  
  
  ##########################################################################################
  ##############################  EOO ####################### 
  message("Extent of occurrences computation")
  EOO <-
    EOO.computing(
      XY = x,
      exclude.area = exclude.area,
      country_map = country_map,
      export_shp = ifelse(DrawMap, TRUE, FALSE),
      alpha = alpha,
      buff.alpha = buff.alpha,
      method.range = method.range,
      write_results = FALSE, 
      parallel = parallel,
      show_progress = show_progress,
      NbeCores = NbeCores, 
      proj_type = proj_type, 
      mode = mode
    ) # , verbose=FALSE
  
  
  ################### AOO estimation #######################################################
  message("Area of occupancy computation")
  AOO <-
    AOO.computing(
      XY = x,
      Cell_size_AOO = Cell_size_AOO,
      nbe.rep.rast.AOO = nbe.rep.rast.AOO,
      export_shp = ifelse(DrawMap, TRUE, FALSE),
      parallel = parallel,
      show_progress = show_progress,
      NbeCores = NbeCores,
      proj_type = proj_type
    )
  
  
  if (any(EOO$eoo - AOO$aoo < 0) || any(is.na(EOO$eoo))) {
    message("Some EOO values are lower than AOO OR null and are thus set equal to AOO")
    
    EOO$issue_eoo[which(EOO$eoo - AOO$aoo < 0)] <-
      paste(EOO$issue_eoo[which(EOO$eoo - AOO$aoo < 0)], "EOO is set equal to AOO because it is lower than AOO", sep = "|")
    
    EOO$eoo[which(EOO$eoo - AOO$aoo < 0)] <- 
      AOO$aoo[which(EOO$eoo - AOO$aoo < 0)]
    
    EOO$issue_eoo[which(is.na(EOO$eoo))] <-
      paste(EOO$issue_eoo[which(is.na(EOO$eoo))], "EOO is set equal to AOO because it is null", sep = "|")
    
    EOO$eoo[which(is.na(EOO$eoo))] <- 
      AOO$aoo[which(is.na(EOO$eoo))]
    
    EOO$issue_eoo <- gsub("NA|", "", EOO$issue_eoo)
    
  }
  
  
  categories <- 
    cat_criterion_b(EOO = EOO$eoo, 
                    AOO = AOO$aoo, 
                    locations = locations_res$locations$locations)
  
  results_full <-
    data.frame(
      taxa = row.names(AOO),
      EOO = EOO$eoo,
      AOO = AOO$aoo,
      locations = locations_res$locations$locations,
      category = categories$ranks_B,
      cat_codes = categories$cats_code,
      subpop = if (SubPop) NbeSubPop$subpop else NA,
      issue_aoo = AOO$issue_aoo,
      issue_eoo = EOO$issue_eoo,
      issue_locations = locations_res$locations$issue_locations,
      main_threat = if (any(colnames(locations_res$locations) == "main_threat")) locations_res$locations$main_threat else NA,
      locations_res$locations[colnames(locations_res$locations) %in% names(threat_list)]
    )
  
  return(results_full)
}






