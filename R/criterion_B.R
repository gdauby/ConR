#' @title Assess IUCN Criterion B
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion B, which is based on species geographic distribution (i.e. extent
#'  of occurrence - EOO, and area of occupancy, AOO)
#'
#' @param x a `dataframe` or an object of class `spgeoIN` see
#' <https://github.com/azizka/speciesgeocodeR>. See Details
#' @param AOO a vector of species AOO, if available
#' @param EOO a vector of species EOO, if available
#' @param locations a vector of species number of locations, if available
#' @param severe.frag a vector indicating if species is severely fragmented, if
#'   available
#' @param subpops a vector of the number of subpopulations (sensu IUCN), if
#'   available
#' @param decline a vector providing the status of the species continuing
#'   decline in EOO, AOO, habitat, locations or subpopulations or population
#'   size (i.e. condition 'b') to be passed on to function `cat_criterion_b()`.
#'   If different of 'Decreasing', the condition 'b' of criterion B will not be
#'   met.
#' @param EOO.threshold numeric vector indicating the thresholds used to
#'   categorize EOO in IUCN categories
#' @param AOO.threshold numeric vector indicating the thresholds used to
#'   categorize AOO in IUCN categories
#' @param Loc.threshold numeric vector indicating the thresholds used to
#'   categorize the number of locations in IUCN categories
#' @param comp_subpop logical if the number of sub-populations should be
#'   computed. By default is TRUE
#' @param comp_severe.frag logical if severe fragmentation should be
#'   computed. By default is FALSE
#' @param method_locations string, indicating the method used for estimating the number of locations. See Details
#'  * `"fixed_grid"` (by default)
#'  * `"sliding_scale"`
#' @inheritParams locations.comp
#' @inheritParams EOO.computing
#' @inheritParams AOO.computing
#' @inheritParams subpop.comp
#' @param threshold_severe numeric of one value indicates the threshold 
#'    to identify severe fragmentation. By default is 50
#' @param DrawMap logical, by default is FALSE, if TRUE, a png map is created for each species 
#'    in a directory of the working environment
#' @param add.legend logical, whether legend should be added to map
#' @inheritParams activate_parallel
#' 
#' @return A data frame containing, for each of taxon, (EOO, AOO, n.locs, n.subpops?),
#'   the IUCN categories associated with the sub-criteria and the consensus category
#'   for criterion B.
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
criterion_B <- function(x = NULL,
                        AOO = NULL,
                        EOO = NULL,
                        locations = NULL,
                        severe.frag = NULL,
                        subpops = NULL,
                        decline = NULL,
                        #add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                        EOO.threshold = c(20000, 5000, 100),
                        AOO.threshold = c(2000, 500, 10),
                        Loc.threshold = c(10, 5, 1),
                        comp_subpop = TRUE,
                        comp_severe.frag = FALSE,
                        resol_sub_pop = 5,
                        cell_size_locations = 10,
                        method_locations = "fixed_grid",
                        method_polygons = "no_more_than_one",
                        rel_cell_size = 0.05,
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
                        threshold_severe = 50,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2,
                        proj_type = "cea",
                        mode = "spheroid",
                        DrawMap = FALSE,
                        add.legend = TRUE) {
  
  if (!requireNamespace("lwgeom", quietly = TRUE))
    stop("The 'lwgeom' package is required to run this function.",
         "Please install it first.")
  
  if (is.null(country_map)) {
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    country_map <- sf::st_make_valid(country_map)
  } else {
    if (any(!st_is_valid(country_map)))
      country_map <- sf::st_make_valid(country_map)
  }
  
  if (is.null(AOO) | is.null(EOO) | is.null(locations) | is.null(severe.frag) | ifelse(comp_subpop, is.null(subpops), FALSE))
    if (is.null(x)) stop("The argument 'x' is needed because at least one of AOO, EOO, locations, severe.frag and subpops is not provided")
  
  if (!is.null(x)) {
    if (identical(class(x)[1], "spgeoIN")) {
      x <- cbind(x$species_coordinates, x$identifier)
      x <- x[, c(2, 1, 3)]
    }
    colnames(x)[1:3] <- c("ddlat", "ddlon", "tax")    
  }
  
  ##########################################################################################
  ##############################  Sub-populations estimations ############################## 
  
  if (comp_subpop) {
  
    if(!is.null(x) & is.null(subpops)) {
      
      if (show_progress) message("Subpopulations computation")
      
      subpop_stats <-
        subpop.comp(
          XY = x,
          resol_sub_pop = resol_sub_pop,
          parallel = parallel,
          show_progress = show_progress,
          NbeCores = NbeCores,
          proj_type = proj_type, 
          export_shp = any(c(DrawMap, comp_severe.frag))
        )
      
      if (DrawMap | comp_severe.frag) {
        
        SubPopPoly <- subpop_stats$poly_subpop
        NbeSubPop <- subpop_stats$number_subpop
        
      } else {
        
        NbeSubPop <- subpop_stats
        
      }
    } else {
      
      NbeSubPop <- subpops
      
      if (!identical(class(NbeSubPop), "data.frame"))
        stop("subpops should be a data.frame similar to the results of the function 'subpop.comp'")
      
      if (!all(colnames(NbeSubPop[,c(1:2)]) == c("tax", "subpop")))
        stop("subpops data.frame first two columns should be named as tax, subpop")
      
      if (DrawMap) {
        DrawMap <- FALSE
        message("Cannot draw map if subpops is provided")
      }
      
    }
      
  }
  
  
  
  ##########################################################################################
  ##############################  Estimations of number of Locations ####################### 
  
  if (!is.null(x) & is.null(locations)) {
    
    if (show_progress) message("Locations computation")
    
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
        rel_cell_size = rel_cell_size,
        parallel = parallel,
        NbeCores = NbeCores, 
        proj_type = proj_type
      )
    
    if (DrawMap) locations_poly <- locations_res$locations_poly
    locations_res <- locations_res$locations
    
  } else {
    
    locations_res <- locations
    
    if (!identical(class(locations_res), "data.frame"))
      stop("locations should be a data.frame similar to the results of the function 'locations.comp'")
    
    if (!all(colnames(locations_res[,c(1:3)]) == c("tax", "locations", "issue_locations")))
      stop("locations data.frame first three columns should be named as tax, locations, issue_locations")
    
    if (DrawMap) {
      DrawMap <- FALSE
      message("Cannot draw map if locations is provided")
    }
    
  }
  
  
  
  ##########################################################################################
  ##############################  EOO ####################### 
  
  if (!is.null(x) & is.null(EOO)) {
    if (show_progress) message("Extent of occurrences computation")
    EOO <-
      EOO.computing(
        XY = x,
        exclude.area = exclude.area,
        country_map = country_map,
        export_shp = DrawMap,
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
      EOO <- EOO$results
      eoo_poly <- EOO$spatial
    }
    
  } else {
    
    if (!identical(class(EOO), "data.frame"))
      stop("EOO should be a data.frame similar to the results of the function 'EOO.computing'")
    
    if (!all(colnames(EOO[,c(1:3)]) == c("tax", "eoo", "issue_eoo")))
      stop("EOO data.frame first three columns should be named as tax, eoo, issue_eoo")
    
    if (DrawMap) {
      DrawMap <- FALSE
      message("Cannot draw map if EOO is provided")
    }
    
  }
  
  ################### AOO estimation #######################################################
  if (!is.null(x) & is.null(AOO)) {
    if (show_progress) message("Area of occupancy computation")
    AOO <-
      AOO.computing(
        XY = x,
        cell_size_AOO = cell_size_AOO,
        nbe.rep.rast.AOO = nbe.rep.rast.AOO,
        export_shp = any(c(DrawMap, comp_severe.frag)),
        parallel = parallel,
        show_progress = show_progress,
        NbeCores = NbeCores,
        proj_type = proj_type
      )
    
    if (DrawMap | comp_severe.frag) {
      
      aoo_poly <- AOO$AOO_poly
      AOO <- AOO$AOO
      
    }
  } else {
    
    if (!identical(class(AOO), "data.frame"))
      stop("AOO should be a data.frame similar to the results of the function 'AOO.computing'")
    
    if (!all(colnames(AOO[,c(1:3)]) == c("tax", "aoo", "issue_aoo")))
      stop("AOO data.frame first three columns should be named as tax, aoo, issue_aoo")
    
    if (DrawMap) {
      DrawMap <- FALSE
      message("Cannot draw map if AOO is provided")
    }
  }
  
  if (comp_severe.frag & !is.null(x)) {
    
    if (is.null(severe.frag)) {
      if (show_progress) message("Assessment of severe fragmentation")
      
      radius <- subpop.radius(XY = x[, c(1:3)], 
                              quant.max = 0.9)
      
      severe.frag <- 
        severe_frag(AOO_poly = aoo_poly, 
                    habitat_poly = SubPopPoly, 
                    dist_isolated = radius, 
                    show_progress = show_progress, 
                    threshold_severe = threshold_severe)
      
    } else {
      
      if (!identical(class(severe.frag), "data.frame"))
        stop("severe.frag should be a data.frame similar to the results of the function 'severe_frag'")
      
      if (!all(colnames(AOO[,c(1:3)]) == c("tax", "frac", "severe_frag")))
        stop("AOO data.frame first three columns should be named as tax, aoo, issue_aoo")
    }
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
                    locations = locations_res$locations, 
                    sever.frag = if (any(c(comp_severe.frag, !is.null(severe.frag)))) severe.frag$severe_frag else NA,
                    decline = if (!is.null(decline)) decline else NA)
  
  if (!is.null(x)) list_data <- coord.check(XY = x)
  
  results_full <-
    data.frame(
      tax = AOO$tax,
      EOO = EOO$eoo,
      AOO = AOO$aoo,
      locations = locations_res$locations,
      nbe_unique_occs = if (!is.null(x)) list_data$unique_occs else NA_integer_,
      category_B = categories$ranks_B,
      category_B_code = categories$cats_code,
      subpop = if (comp_subpop) NbeSubPop$subpop else NA_real_,
      severe_frag = if (any(c(comp_severe.frag, !is.null(severe.frag)))) severe.frag$severe_frag else NA,
      issue_aoo = AOO$issue_aoo,
      issue_eoo = EOO$issue_eoo,
      issue_locations = locations_res$issue_locations,
      main_threat = if (any(colnames(locations_res) == "main_threat")) locations_res$main_threat else NA_character_,
      locations_res[colnames(locations_res) %in% names(threat_list)]
    )
  
  row.names(results_full) <- NULL
  
  if (!is.null(x)) {
    
    if (DrawMap) {
      
      poly_borders <- 
        suppressWarnings(sf::st_crop(country_map, c(xmin = min(x[,2]), ymin = min(x[,1]), xmax = max(x[,2]), ymax = max(x[,1]))))
      poly_borders <- st_make_valid(poly_borders)
      
      for (i in 1:length(list_data$list_data)) {
        
        name_sp = results_full$tax[i]
        
        draw_map_cb(XY = list_data$list_data[[i]], 
                    name_sp = name_sp, 
                    eoo_poly = eoo_poly[which(eoo_poly$tax == name_sp),], 
                    aoo_poly = aoo_poly[which(aoo_poly$tax == name_sp),], 
                    locations_poly = locations_poly[which(locations_poly$tax == name_sp),], 
                    subpop = SubPopPoly[which(SubPopPoly$tax == name_sp),], 
                    proj_type = proj_crs(proj_type = proj_type), 
                    results = results_full[which(results_full$tax == name_sp),], 
                    add.legend = add.legend,
                    poly_borders = poly_borders)
        
      }
      
      # name_sp <- results_full$tax[1]
      # eoo_poly <- eoo_poly[which(eoo_poly$tax == name_sp),]
      # aoo_poly <- aoo_poly[which(aoo_poly$tax == name_sp),]
      # locations_poly <- locations_res$locations_poly[which(locations_res$locations_poly$tax == name_sp),]
      # XY <- x[which(x[,3] == name_sp),]
      # SubPopPoly <- SubPopPoly[which(SubPopPoly$tax == name_sp),]
      # 
      # draw_map_cb(XY = x[which(x[,3] == name_sp),], 
      #             name_sp = results_full$tax[1], 
      #             eoo_poly = eoo_poly[which(eoo_poly$tax == name_sp),], 
      #             aoo_poly = aoo_poly[which(aoo_poly$tax == name_sp),], 
      #             locations_poly = locations_res$locations_poly[which(locations_res$locations_poly$tax == name_sp),])
      
    }
  }
  
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
#' @param proj_type crs
#' @param results data.frame
#' @param add.legend logical
#' 
#' @importFrom grDevices rgb
#' @importFrom graphics mtext layout legend box
#' @import sf
#' 
#' @return a plot
#' 
#' @keywords internal
#'
#'
draw_map_cb <- function(XY, 
                        name_sp, 
                        eoo_poly, 
                        aoo_poly, 
                        locations_poly, 
                        subpop, 
                        proj_type, 
                        results, 
                        add.legend =TRUE,
                        poly_borders) {
  
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
  
  if (add.legend) nf <-
    graphics::layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), c(4, 1.5), c(4, 1.5))
  
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
  
  plot(st_geometry(locations_poly),
       add = T,
       col = rgb(
         red = 1,
         green = 0,
         blue = 0,
         alpha = 0.2
       ))
  
  if (nrow(poly_borders) > 0)
    plot(
      st_geometry(poly_borders),
      col = NA,
      axes = FALSE,
      add = T
    )
  
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
  # scaleBAR <- round(matrix(unlist(sf::sf_project(
  #   from = sf::st_crs(4326),
  #   to =
  #     proj_type,
  #   pts = border_to_center
  # )), ncol = 2)/ 1000, 0)[1, 1]
  
  # plot(
  #   terra::rast(matrix(NA, nrow=2, ncol=2)),
  #   col = NA,
  #   add = T
  # )
  # 
  # terra::sbar(d = scaleBAR, 
  #             type = "bar", 
  #             below = "kilometers", 
  #             cex=2.5, xy = "bottomright")
  graphics::mtext(name_sp, side=3, cex=3, line=3)
  
  if (add.legend) {
    graphics::par(mar=c(1,1,1,1), xpd=T)
    plot(1:10, 1:10, type="n", bty='n', xaxt='n', yaxt='n')

    graphics::legend(1,10,  c(paste("EOO=", ifelse(!is.na(results$EOO), format(round(as.numeric(results$EOO),1), scientific = 5), NA), "km2"),
                      paste("AOO=", format(results$AOO, scientific = 5),"km2"),
                      # paste("Number of unique occurrences=", results$),
                      paste("Number of sub-populations=", results$subpop),
                      paste("Number of locations=", results$locations),
                      paste("IUCN category according to criterion B:", results$category_B)), cex=3.5, bg = grDevices::grey(0.9))
    
    # graphics::par(mar=c(4,1,1,1))
    # plot(full_poly_borders, lty=1, lwd=1,axes=FALSE)
    # graphics::points(XY[,1],XY[,2], pch=8, cex=2, col="red")
  }
  
  
  
  grDevices::dev.off()
}

