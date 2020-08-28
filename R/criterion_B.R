#' @title Assess IUCN Criterion B
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion B, which is based on species geographic distribution (i.e. extent
#'  of occurrence - EOO, and area of occupancy, AOO)
#'
#' @param x a \code{dataframe} or an object of class \code{spgeoIN} see
#' \url{https://github.com/azizka/speciesgeocodeR}. See Details
#' @param protec.areas a \code{SpatialPolygonsDataFrame}, shapefile with
#' protected areas.  If provided, this will be taken into account for
#' calculating number of location (see Details and
#' \code{method_protected_area}).  By default, no shapefile is provided
#' @param EOO.threshold numeric vector indicating the thresholds used to categorize EOO in IUCN categories
#' @param AOO.threshold numeric vector indicating the thresholds used to categorize AOO in IUCN categories
#' @param Loc.threshold numeric vector indicating the thresholds used to categorize the number of locations in IUCN categories
#' @param SubPop logical if the number of sub-populations should be computed. By default is TRUE
#' @param Resol_sub_pop a numeric, value indicating the radius size in
#' kilometers used for estimating the number of sub-population. By default,
#' equal to 5
#' @param Cell_size_locations a numeric, value indicating the grid size in
#' kilometers used for estimating the number of location. By default, equal to
#' 10
#' @param method_locations a character string, indicating the method used for
#' estimating the number of locations.  "fixed_grid" or "sliding scale". See
#' details. By default, it is "fixed_grid"
#' @param Rel_cell_size a numeric, if \code{method_locations="sliding scale"},
#' \code{Cell_size_locations} is ignored and the resolution is given by the
#' maximum distance separating two occurrences multiplied by
#' \code{Rel_cell_size}. By default, it is 0.05
#' @param method_protected_area a character string. By default is
#' "no_more_than_one"", which means occurrences within protected areas (if
#' provided) will not be taken into account for estimating the number of
#' locations following the grid system, see Details. By default, it is
#' "no_more_than_one"
#' @param ID_shape_PA a character string, indicating the field name of
#' \code{protec.areas} with ID of the \code{SpatialPolygonsDataFrame} of
#' protected areas
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
#' @export criterion_B
#'
#' @examples
#' 
#' @importFrom tibble is_tibble
#' @importFrom rnaturalearth ne_countries
#' @importFrom rgeos gBuffer
#' @importFrom sp identicalCRS proj4string CRS
#' 
#'
criterion_B <- function(x, 
                       protec.areas = NULL,
                       #add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                       EOO.threshold = c(20000, 5000, 100), 
                       AOO.threshold = c(2000, 500, 10), 
                       Loc.threshold = c(10, 5, 1),
                       SubPop = TRUE,
                       Resol_sub_pop = 5,
                       Cell_size_locations = 10,
                       method_locations = "fixed_grid",
                       Rel_cell_size = 0.05,
                       method_protected_area = "no_more_than_one",
                       ID_shape_PA = "WDPA_PID",
                       country_map = NULL,
                       method.range = "convex.hull",
                       alpha = 1,
                       buff.alpha = 0.1,
                       exclude.area = FALSE,
                       buff_width = 0.1,
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
    
  }else{
    
    if(any(grepl('sf', class(country_map))))
      country_map <- 
        as(country_map, "Spatial")
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
    country_map <- 
      as(country_map, "sf")
  }
  
  if (!is.null(protec.areas)) {
    if (!sp::identicalCRS(protec.areas, country_map)) {
      
      sp::proj4string(protec.areas) <- 
        sp::CRS(SRS_string = 'EPSG:4326')
      sp::proj4string(country_map) <-
        sp::CRS(SRS_string = 'EPSG:4326')
      
    }
  }
  
  
  ##########################################################################################
  ##############################  Sub-populations estimations ############################## 
  
  if(SubPop) {
    
    message("Subpopulations computation")
    
    subpop_stats <-
      subpop.comp(
        x,
        Resol_sub_pop = 10,
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores, 
        proj_type = proj_type
      )
    
    SubPopPoly <- subpop_stats$number_subpop
    NbeSubPop <- subpop_stats$poly_subpop
    
  }else{
    
    SubPopPoly <- vector(mode = "numeric", length = length(unique(x$tax)))
    
  }
  
  ##########################################################################################
  ##############################  Estimations of number of Locations ####################### 
  
  message("Locations computation")
  
  locations_res <-
    locations.comp(
      XY = x,
      method = method_locations,
      protec.areas = protec.areas,
      method_protected_area = method_protected_area,
      Cell_size_locations = Cell_size_locations,
      ID_shape_PA = ID_shape_PA,
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
      buff_width = buff_width,
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
  
  
  EOO <- EOO$EOO
  
  categories <- 
    cat_criterion_b(EOO = EOO, AOO = AOO, locations = locations_res[[2]])
  
  results_full <-
    data.frame(
      taxa = names(AOO),
      EOO = EOO,
      AOO = AOO,
      locations = locations_res[[2]],
      category = categories$ranks_B12a,
      cat_codes = categories$cat_codes,
      subpop = SubPopPoly
    )
  
  return(results_full)
}






