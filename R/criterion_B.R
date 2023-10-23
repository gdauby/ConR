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
#' @inheritParams subpop.comp
#' 
#' 
#' @return A data frame containing, for each of taxon, (EOO, AOO, n.locs, n.subpops?),
#'   the IUCN categories associated with the sub-criteria and the consensus category
#'   for criterion B.
#' 
#' @details The function ... 
#' 
#' @author Gilles Dauby & Renato A. Ferreira de Lima
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#'
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
                       cell_size_locations = 10,
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
                       cell_size_AOO = 2,
                       nbe.rep.rast.AOO = 0,
                       parallel = FALSE,
                       show_progress = TRUE,
                       NbeCores = 2,
                       proj_type = "cea",
                       mode = "spheroid",
                       DrawMap = FALSE) {
  
  if(identical(class(x)[1], "spgeoIN")) {
    x <- cbind(x$species_coordinates, x$identifier)
    x <- x[,c(2,1,3)]
  }
  colnames(x)[1:3] <- c("ddlat","ddlon","tax")
  
  if (!requireNamespace("lwgeom", quietly = TRUE))
    stop(
      "The 'lwgeom' package is required to run this function.",
      "Please install it first."
    )
  
  if(tibble::is_tibble(x)) x <- as.data.frame(x)
  
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  } else {
    
    if (any(!st_is_valid(country_map)))
      country_map <- sf::st_make_valid(country_map)
  }
  
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
      cell_size_locations = cell_size_locations, 
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
      parallel = parallel,
      show_progress = show_progress,
      NbeCores = NbeCores, 
      proj_type = proj_type, 
      mode = mode
    ) # , verbose=FALSE
  
  if (DrawMap) {
    
    eoo_poly <- EOO$spatial
    EOO <- EOO$results
    
  }
  
  
  ################### AOO estimation #######################################################
  message("Area of occupancy computation")
  AOO <-
    AOO.computing(
      XY = x,
      cell_size_AOO = cell_size_AOO,
      nbe.rep.rast.AOO = nbe.rep.rast.AOO,
      export_shp = ifelse(DrawMap, TRUE, FALSE),
      parallel = parallel,
      show_progress = show_progress,
      NbeCores = NbeCores,
      proj_type = proj_type
    )
  
  if (DrawMap) {
    
    aoo_poly <- AOO$AOO_poly
    
    
    AOO <- AOO$AOO
    
  }
  
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
      tax = AOO$tax,
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
  
  list_data <- coord.check(XY = x)
  
  
  for (i in 1:length(list_data$list_data)) {
    
    name_sp = results_full$tax[i]
    
    draw_map_cb(XY = list_data$list_data[[i]], 
                name_sp = name_sp, 
                eoo_poly = eoo_poly[which(eoo_poly$tax == name_sp),], 
                aoo_poly = aoo_poly[which(aoo_poly$tax == name_sp),], 
                locations_poly = locations_res$locations_poly[which(locations_res$locations_poly$tax == name_sp),], 
                subpop = SubPopPoly[which(SubPopPoly$tax == name_sp),])
    
  }
  
  
  name_sp <- results_full$tax[1]
  eoo_poly <- eoo_poly[which(eoo_poly$tax == name_sp),]
  aoo_poly <- aoo_poly[which(aoo_poly$tax == name_sp),]
  locations_poly <- locations_res$locations_poly[which(locations_res$locations_poly$tax == name_sp),]
  XY <- x[which(x[,3] == name_sp),]
  SubPopPoly <- SubPopPoly[which(SubPopPoly$tax == name_sp),]
    
  draw_map_cb(XY = x[which(x[,3] == name_sp),], 
              name_sp = results_full$tax[1], 
              eoo_poly = eoo_poly[which(eoo_poly$tax == name_sp),], 
              aoo_poly = aoo_poly[which(aoo_poly$tax == name_sp),], 
              locations_poly = locations_res$locations_poly[which(locations_res$locations_poly$tax == name_sp),])
  
  return(results_full)
}




#' Internal function
#'
#' Coordinates check
#'
#' @param XY data.frame, of at least two columns (coordinates), third is taxa
#' @param name_sp character vector
#' @param eoo_poly sf
#' @param aoo_poly sf
#' @param locations_poly sf
#' @param subpop sf
#' 
#' @return a plot
#' 
#' @keywords internal
#'
#'
draw_map_cb <- function(XY, name_sp, eoo_poly, aoo_poly, locations_poly, subpop) {
  
  name_file <-
    paste("IUCN_", gsub(pattern = " ", replacement = "_", name_sp), sep =
            "")
  
  directory_name <- 
    "IUCN_"
  
  dir.create(file.path(paste(
    getwd(), paste("/", directory_name, "criterion_b_map", sep = ""), sep = ""
  )), showWarnings = FALSE)
  
  grDevices::png(
    paste(file.path(paste(
      getwd(), paste("/", 
                     directory_name, 
                     "criterion_b_map", sep = ""), sep = ""
    )), "/", name_file, ".png", sep = "")
    ,
    width = 2000,
    height = 2000
  )
  
  graphics::par(mar = c(10, 12, 10, 2),
                xpd = FALSE,
                las = 1)
  
  XY_sf <- st_as_sf(XY, coords = c("ddlon", "ddlat"))
  st_crs(XY_sf) <- 4326
  
  plot(
    st_geometry(XY_sf),
    pch = 19,
    cex = 2,
    col = rgb(
      red = 0,
      green = 0,
      blue = 0,
      alpha = 0.2
    )
  )
  
  if (!is.null(eoo_poly)) {
    plot(
      st_geometry(eoo_poly),
      col = rgb(
        red = 0.2,
        green = 0.2,
        blue = 0.2,
        alpha = 0.1
      ),
      extent = XY_sf,
      add= TRUE
    )
  }
  
  plot(st_geometry(subpop), 
       add =T, 
       border="black", 
       lwd=1, lty=1, 
       extent = XY_sf)
  
  plot(locations_poly,
       add = T,
       col = rgb(
         red = 1,
         green = 0,
         blue = 0,
         alpha = 0.2
       ))
  
  
  graphics::axis(1, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
  graphics::axis(1, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
  graphics::axis(2, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
  graphics::axis(2, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
  graphics::box()
  
  
  xlim <- graphics::par("xaxp")[1:2]
  xlim <- abs(xlim[1]-xlim[2])
  border_to_center <- as.data.frame(matrix(NA, 2, 2))
  border_to_center[,1] <- c(xlim/10, 0)
  border_to_center[,2] <- c(0,0)
  scaleBAR <- round(matrix(unlist(sf::sf_project(
    from = sf::st_crs(4326),
    to =
      proj_crs(proj_type),
    pts = border_to_center
  )), ncol = 2)/ 1000, 0)[1, 1]
  
  plot(
    terra::rast(matrix(NA, nrow=2, ncol=2)),
    col = NA,
    add = T
  )
  
  terra::sbar(d = scaleBAR, 
              type = "bar", 
              below = "kilometers", 
              cex=2.5, xy = "bottomright")
  mtext(name_sp, side=3, cex=3, line=3)
  
  grDevices::dev.off()
}

