#' @title Internal function
#'
#' @description Compute IUCN eval
#' 
#' @param DATA data.frame
#' @param poly_borders SpatialPolygonDataFrame
#' @param Cell_size_AOO integer
#' @param Cell_size_locations integer
#' @param Resol_sub_pop integer
#' @param method_locations integer
#' @param Rel_cell_size integer
#' @param protec.areas SpatialPolygonDataFrame
#' @param exclude.area logical
#' @param method_protected_area string
#' @param ID_shape_PA string
#' @param buff_width numeric
#' @param NamesSp string
#' @param write_shp logical
#' @param file_name string
#' @param add.legend logical
#' @param DrawMap logical
#' @param map_pdf logical
#' @param draw.poly.EOO logical
#' @param SubPop logical
#' @param MinMax numeric vector
#' @param alpha integer
#' @param buff.alpha numeric
#' @param method.range string
#' @param nbe.rep.rast.AOO integer
#' @param verbose logical
#' @param showWarnings logical
#' 
#' @importFrom rgdal project writeOGR
#' @importFrom rnaturalearth ne_countries
#' 
IUCN.comp <- function(DATA,
                       poly_borders = NULL,
                       Cell_size_AOO = 2,
                       Cell_size_locations = 10,
                       Resol_sub_pop = 5,
                       method_locations = c("fixed_grid"),
                       Rel_cell_size = 0.05,
                       protec.areas = NULL,
                       exclude.area = FALSE,
                       method_protected_area = "no_more_than_one",
                       ID_shape_PA = "WDPA_PID",
                       buff_width = 0.1,
                       NamesSp = "species1",
                       write_shp = FALSE,
                       file_name = NULL,
                       add.legend = TRUE,
                       DrawMap = TRUE,
                       map_pdf = FALSE,
                       draw.poly.EOO = TRUE,
                       SubPop = TRUE,
                       MinMax,
                       alpha = 1,
                       buff.alpha = 0.1,
                       method.range = "convex.hull",
                       nbe.rep.rast.AOO = 0,
                       verbose = TRUE,
                       showWarnings = TRUE)
{
  
  
  ### Getting by default land map if poly_borders is not provided
  if (is.null(poly_borders)) {
    
    land <- 
      rnaturalearth::ne_countries(scale = 50, returnclass = "sp")
    
    # data('land', package = 'ConR', envir = environment())
    # land <- get("land", envir = environment())
    # data(land, envir = environment())
    poly_borders <- land
  }
  
  ### cropping poly_borders according to range of occurrences shapefile for producing lighter maps
  if (DrawMap) {
    full_poly_borders <- poly_borders
    if (!is.null(poly_borders))
      poly_borders <- suppressWarnings(raster::crop(poly_borders, raster::extent(MinMax) + 30))
  }
  
  projEAC <- .proj_crs()
  
  ## Initialization of data frame for stocking results
  if(is.null(protec.areas)){
    Results <- as.data.frame(matrix(NA,9,1))
    colnames(Results) <- NamesSp
    rownames(Results) <- c("EOO","AOO","Nbe_unique_occ.", "Nbe_subPop", "Nbe_loc", "Category_CriteriaB", "Category_code",
                           "Category_AOO","Category_EOO")
  }else{
    Results <- as.data.frame(matrix(NA,11,1))
    colnames(Results) <- NamesSp
    rownames(Results) <- c("EOO","AOO","Nbe_unique_occ.", "Nbe_subPop", "Nbe_loc", "Nbe_loc_PA", 
                           "Category_CriteriaB", "Category_code",  "Ratio_occ_within_PA",
                           "Category_AOO","Category_EOO")
  }
  
  # if(verbose) cat(" ",paste("Evaluation of", as.character(NamesSp)))
  
  XY <- DATA[,c(2:1)]
  ### Projection into equal area cylindrical
  coordEAC <- as.data.frame(matrix(unlist(rgdal::project(as.matrix(XY),proj=as.character(projEAC),inv=FALSE)), ncol=2))
  rownames(coordEAC) <- seq(1,nrow(coordEAC),1)
  
  ##########################################################################################
  ##############################  Sub-populations estimations ############################## 
  
  if(SubPop) {
    subpop_stats <- subpop.comp(DATA, Resol_sub_pop=Resol_sub_pop)
    SubPopPoly <- subpop_stats[[2]]
    NbeSubPop <- subpop_stats[[1]]
  }
  
  ##########################################################################################
  ##############################  Estimations of number of Locations ####################### 
  
  # ## range of lat and long
  # Corners <- rbind(c(min(XY[,1]), max(XY[,1])), c(min(XY[,2]), max(XY[,2])))
  
  locations_res <-
    locations.comp(
      XY = DATA,
      method = method_locations,
      protec.areas = protec.areas,
      method_protected_area = method_protected_area,
      Cell_size_locations = Cell_size_locations,
      ID_shape_PA = ID_shape_PA,
      show_progress = FALSE,
      Rel_cell_size = Rel_cell_size
    )
  
  if (is.null(protec.areas)) {
    
    r2 <- locations_res[[1]][[1]]
    Locations <- locations_res[[2]]
    
  } else{
    
    r2 <- locations_res[[1]][[1]]
    r2_PA <- locations_res[[2]]
    LocNatParks <- locations_res[[3]]
    LocOutNatParks <- locations_res[[4]]
    
  }
  
  if(!is.null(protec.areas)) {
    DATA_SF <- as.data.frame(unique(XY))
    colnames(DATA_SF) <- c("ddlon","ddlat")
    sp::coordinates(DATA_SF) <-  ~ddlon+ddlat
    
    sp::proj4string(DATA_SF) <- sp::CRS(SRS_string='EPSG:4326')
    
    
    # raster::crs(DATA_SF) <- raster::crs(protec.areas)
    Links_NatParks <- sp::over(DATA_SF, protec.areas)    
  }
  
  if(any(method_locations=="fixed_grid")) Resolution <- Cell_size_locations*1000
  if(any(method_locations=="sliding scale")){
    pairwise_dist <- stats::dist(coordEAC,  upper = F)
    if(nrow(coordEAC) > 1) {Resolution <- max(pairwise_dist)*Rel_cell_size
    }else{
      Resolution <- 10000
    }
  }
  
  if(nrow(unique(XY)) > 2) { ### if more than 2 uniques occurrences
    
    ##############################  EOO estimation ##############################      
    
    EOO_ <-
      EOO.computing(
        DATA[, 1:2],
        exclude.area = exclude.area,
        country_map = poly_borders,
        Name_Sp = NamesSp,
        buff_width = buff_width,
        export_shp = TRUE,
        alpha = alpha,
        buff.alpha = buff.alpha,
        method.range = method.range,
        write_results = FALSE, 
        show_progress = FALSE
      ) # , verbose=FALSE
    
    p1 <- EOO_[[2]]
    EOO <- EOO_[[1]]
    
    ################### AOO estimation #######################################################
    
    AOO <- 
      .AOO.estimation(coordEAC, 
                      cell_size = Cell_size_AOO, 
                      nbe_rep = nbe.rep.rast.AOO)
    # , 
    # poly_borders = poly_borders
    
    if(EOO<AOO) EOO <- AOO ### If EOO is < AOO, EOO is put equal to AOO
    
    ## recording results
    Results["EOO",1] <- as.numeric(EOO)
    Results["AOO",1] <- as.numeric(AOO)
    if(SubPop) Results["Nbe_subPop",1] <- as.numeric(NbeSubPop)
    Results["Nbe_unique_occ.",1] <- nrow(unique(XY))
    
    if(!is.null(protec.areas)) Results["Nbe_loc",1] <- as.numeric(LocNatParks + LocOutNatParks)
    if(is.null(protec.areas)) Results["Nbe_loc",1] <- as.numeric(Locations)
    
    if(!is.null(protec.areas)) {Results["Nbe_loc_PA",1] <- LocNatParks}
    
    if(!is.null(protec.areas)) Results["Ratio_occ_within_PA",1] <- round(length(which(!is.na(Links_NatParks[,1])))/nrow(Links_NatParks)*100,1)
    Nbe_Loc <- as.numeric(Results["Nbe_loc",1])
    
    ### Criteria B assessment following IUCN thresholds   ##################
    if (EOO < 20000) {
      Rank_EOO <- 3
      if (EOO < 5000) {
        Rank_EOO <- 2
        if (EOO < 100) {
          Rank_EOO <- 1
        }
      }
    } else
      (Rank_EOO <- 4)
    
    if(AOO<2000){
      Rank_AOO <- 3
      if(AOO<500){
        Rank_AOO <- 2
        if(AOO>10){
          Rank_AOO <-1
        }}}else{Rank_AOO <- 4}
    
    if (Nbe_Loc <= 10) {
      Rank_Loc <- 3
      if (Nbe_Loc <= 5) {
        Rank_Loc <- 2
        if (Nbe_Loc == 1) {
          Rank_Loc <- 1
        }
      }
    } else{
      Rank_Loc <- 4
    }
    
    Rank_B1a <- max(Rank_EOO, Rank_Loc)
    Rank_B2a <- max(Rank_AOO, Rank_Loc)
    Rank_CriteriaB <- min(Rank_B1a, Rank_B2a)
    
    if (Rank_CriteriaB == 1)
      Cat <- "CR"
    if (Rank_CriteriaB == 2)
      Cat <- "EN"
    if (Rank_CriteriaB == 3 && Nbe_Loc > 0 &&
        Nbe_Loc < 11)
      Cat <- "VU"
    
    if (Rank_CriteriaB > 3 && Nbe_Loc >= 0)
      Cat <- "LC or NT" ###
    
    if(Rank_B1a>Rank_B2a) Cat_Code <- paste(Cat,"B2a")
    if(Rank_B1a<Rank_B2a) Cat_Code <- paste(Cat,"B1a")
    if(Rank_B1a==Rank_B2a) Cat_Code <- paste(Cat,"B1a+B2a")
    
    if(!is.null(protec.areas)) {
      if(as.numeric(Results["Ratio_occ_within_PA",1])==100){
        Results["Category_CriteriaB",1] <- "LC or NT"
        Results["Category_code",1] <- Cat_Code 
      }else{
        Results["Category_CriteriaB",1] <- Cat
        Results["Category_code",1] <- Cat_Code
      }
    }else{
      Results["Category_CriteriaB",1] <- Cat
      Results["Category_code",1] <- Cat_Code
    }
    
    if(Rank_B2a==1) Results["Category_AOO",1] <- "CR"
    if(Rank_B2a==2) Results["Category_AOO",1] <- "EN"
    if(Rank_B2a==3) Results["Category_AOO",1] <- "VU"
    if(Rank_B2a>3) Results["Category_AOO",1] <- "LC or NT"
    
    if(Rank_B1a==1) Results["Category_EOO",1] <- "CR"
    if(Rank_B1a==2) Results["Category_EOO",1] <- "EN"
    if(Rank_B1a==3) Results["Category_EOO",1] <- "VU"
    if(Rank_B1a>3) Results["Category_EOO",1] <- "LC or NT"
    
    
    p1 <- as(p1, "SpatialPolygonsDataFrame")
    
  }else{ ### if less than 3 uniques occurrences
    
    p1 <- NULL ## EOO shapefile is NULL
    
    if (nrow(coordEAC) == 2) {
      ## if two uniques occurrences
      pairwise_dist <- stats::dist(coordEAC,  upper = F)
      
      if (pairwise_dist <= Resolution) {
        AOO <-
          Cell_size_AOO * Cell_size_AOO ## 1 occupied cell if distance <= resolution
      } else{
        AOO <-
          2 * Cell_size_AOO * Cell_size_AOO ## 2 occupied cells if distance > resolution
      }
    } else{
      AOO <-
        Cell_size_AOO * Cell_size_AOO ## 1 occupied cell if one unique occurrence
    }
    
    Results["AOO",1] <- AOO
    if(SubPop) Results["Nbe_subPop",1] <- NbeSubPop
    Results["Nbe_unique_occ.",1] <- nrow(unique(XY))
    
    if(!is.null(protec.areas)) Results["Ratio_occ_within_PA",1] <- round(length(which(!is.na(Links_NatParks[,1])))/nrow(Links_NatParks)*100,2)
    
    if(is.null(protec.areas)) Results["Nbe_loc",1] <- as.numeric(Locations)
    if(!is.null(protec.areas)) Results["Nbe_loc",1] <- LocNatParks + LocOutNatParks
    if(!is.null(protec.areas)) Results["Nbe_loc_PA",1] <- LocNatParks
    
    
    if(!is.null(protec.areas)){
      if(as.numeric(Results["Ratio_occ_within_PA",1]) == 100){
        ### If all occurences are found within protected areas, the species is considered as not threatened
        Results["Category_CriteriaB",1] <- "LC or NT"
      }else{
        if(as.numeric(Results["AOO",1]) < 10 & as.numeric(Results["Nbe_loc",1]) == 1) {
          Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "CR"
          
        }else{
          if(as.numeric(Results["AOO",1]) < 500 & as.numeric(Results["Nbe_loc",1]) < 6) {
            Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "EN"
            
          }else{
            if(as.numeric(Results["AOO",1]) < 2000 & as.numeric(Results["Nbe_loc",1]) < 11) {
              Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "VU"
              
            }else{
              Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "LC or NT"
            }
          }
        }
      }
    }else{
      
      if(as.numeric(Results["AOO",1]) < 10 & as.numeric(Results["Nbe_loc",1])==1) {
        
        Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "CR"
        
      }else{
        
        if(as.numeric(Results["AOO",1]) < 500 & as.numeric(Results["Nbe_loc",1])<6) {
          
          Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "EN"
          
        }else{
          
          if(as.numeric(Results["AOO",1]) < 2000 & as.numeric(Results["Nbe_loc",1])<11) {
            
            Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "VU"
            
          }else{
            
            Results["Category_AOO",1] <- Results["Category_CriteriaB",1] <- "LC or NT"
            
          }
        }
      }
    }
    
    if(is.na(Results["Category_AOO",1])) {
      if(as.numeric(Results["AOO",1]) < 10 & as.numeric(Results["Nbe_loc",1]) == 1) {
        
        Results["Category_AOO",1] <- "CR"
        
      }else{
        if(as.numeric(Results["AOO",1]) < 500 & as.numeric(Results["Nbe_loc",1]) < 6) {
          
          Results["Category_AOO",1] <- "EN"
          
        }else{
          if(as.numeric(Results["AOO",1]) < 2000 & as.numeric(Results["Nbe_loc",1]) < 11) {
            
            Results["Category_AOO",1] <- "VU"
            
          }else{
            
            Results["Category_AOO",1] <- "LC or NT"
            
          }
        }
      }
    }
    
    Results["Category_code", 1] <-
      paste(Results["Category_CriteriaB", 1], "B2a")
    if (showWarnings)
      message(paste(
        "\nEOO parameter cannot be estimated for",
        NamesSp,
        "because there is less than 3 records"
      ))
  } ## End less than 3 records
  
  ############ Map ###########
  if(DrawMap) {
    
    ## pdf or png format initialization
    if (!map_pdf) {
      if (!is.null(file_name)) {
        NAME_FILE <-
          paste(file_name, gsub(" ", replacement = "_", as.character(NamesSp)) , sep =
                  "")
      } else{
        NAME_FILE <-
          paste("IUCN_", gsub(" ", replacement = "_", as.character(NamesSp)), sep =
                  "")
      }
      
      NAME_FILE <- gsub("[[:punct:]]", " ", NAME_FILE)
      
      DIRECTORY_NAME <- 
        ifelse(!is.null(file_name), file_name, "IUCN_")
      
      dir.create(file.path(paste(
        getwd(), paste("/", DIRECTORY_NAME, "results_map", sep = ""), sep = ""
      )), showWarnings = FALSE)
      
    }
    
    if (!map_pdf)
      grDevices::png(
        paste(file.path(paste(
          getwd(), paste("/", 
                         DIRECTORY_NAME, 
                         "results_map", sep = ""), sep = ""
        )), "/", NAME_FILE, ".png", sep = ""),
        width = 2000,
        height = 2000
      )
    
    ### Layout of the map
    graphics::par(mar = c(10, 12, 10, 2),
                  xpd = FALSE,
                  las = 1)
    if (add.legend &
        !any(colnames(DATA) == "coly"))
      nf <-
      layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), c(4, 1.5), c(4, 1.5))
    if (any(colnames(DATA) == "coly") &
        add.legend)
      nf <-
      layout(matrix(c(1, 1, 1, 1, 1, 1, 2, 3, 4), 3, 3, byrow = TRUE), c(2, 1.5, 1.5), c(4, 1.5, 1.5))
    
    ### Mapping 
    if(!is.null(protec.areas)) {
      
      if (LocOutNatParks == 0) {
        sp::plot(
          poly_borders,
          xlim = c(range(XY[, 1])[1] - 1, range(XY[, 1])[2] + 1),
          ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1),
          axes = FALSE,
          xlab = "",
          ylab = ""
        )
        
      } else{
        
        # r2_pol <- r2
        if(LocOutNatParks == 1){
          
          sp::plot(
            r2,
            col = rgb(
              red = 1,
              green = 0,
              blue = 0,
              alpha = 0.2
            ),
            xlim = c(range(XY[, 1])[1] - 1, range(XY[, 1])[2] + 1),
            ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1)
          )
        }else{
          
          sp::plot(
            r2,
            col = rgb(
              red = 1,
              green = 0,
              blue = 0,
              alpha = 0.2
            ),
            xlim = c(range(XY[, 1])[1] - 1, range(XY[, 1])[2] + 1),
            ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1)
          )
          
        }
      }
    }else{
      # r2_pol <- rasterToPolygons(r2, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=FALSE)
      sp::plot(
        r2,
        col = rgb(
          red = 1,
          green = 0,
          blue = 0,
          alpha = 0.2
        ),
        xlim = c(range(XY[, 1])[1] - 1, range(XY[, 1])[2] + 1),
        ylim = c(range(XY[, 2])[1] - 1, range(XY[, 2])[2] + 1)
      )
    }
    
    if(SubPop) sp::plot(SubPopPoly, add=T, border="black", lwd=2, lty=1)
    
    if (!is.null(protec.areas)) {
      if (LocNatParks > 0) {
        if (method_protected_area != "no_more_than_one") {
          # r2_PA_pol <- rasterToPolygons(r2_PA, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=FALSE)
          sp::plot(r2_PA[[1]],
                   add = T,
                   col = rgb(
                     red = 0,
                     green = 0,
                     blue = 1,
                     alpha = 0.2
                   ))
        }
      }
    }
    
    if (!is.null(p1) &
        draw.poly.EOO)
      sp::plot(p1,
               add = T,
               col = rgb(
                 red = 0.2,
                 green = 0.2,
                 blue = 0.2,
                 alpha = 0.1
               ))
    
    sp::plot(
      poly_borders,
      axes = FALSE,
      lty = 1,
      add = T,
      lwd = 1
    )
    
    if (!is.null(protec.areas))
      sp::plot(
        protec.areas,
        add = T,
        col = rgb(
          red = 0.2,
          green = 0.2,
          blue = 0.2,
          alpha = 0.05
        ),
        lty = 2
      )
    
    if(!is.null(protec.areas)){
      
      colnames(XY) <- c("ddlon", "ddlat")
      XY_sp <- XY[which(is.na(Links_NatParks[, 1])), ]
      if (nrow(XY_sp) > 0) {
        sp::coordinates(XY_sp) <-  ~ ddlon + ddlat
        sp::plot(
          XY_sp,
          pch = 19,
          cex = 2,
          col = "black",
          add = T
        )
      }
      XY_sp <- XY[which(!is.na(Links_NatParks[, 1])), ]
      if (nrow(XY_sp) > 0) {
        sp::coordinates(XY_sp) <-  ~ ddlon + ddlat
        sp::plot(
          XY_sp,
          pch = 19,
          cex = 2,
          col = "blue",
          add = T
        )
      }
    }else{
      colnames(XY) <- c("ddlon", "ddlat")
      XY_sp <- XY
      sp::coordinates(XY_sp) <-  ~ ddlon + ddlat
      sp::plot(
        XY_sp,
        pch = 19,
        cex = 2,
        col = "black",
        add = T
      )
    }
    
    graphics::axis(1, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
    graphics::axis(1, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
    graphics::axis(2, outer=FALSE, cex.axis=3, tick = FALSE, line=1.5)  #pos=min(range(XY[,2]))-2)
    graphics::axis(2, outer=FALSE,labels=FALSE, cex.axis=3, tick = TRUE, line=0)  #pos=min(range(XY[,2]))-2)
    graphics::box()
    
    if(Results["Nbe_loc",1]>1) {
      xlim <- graphics::par("xaxp")[1:2]
      xlim <- abs(xlim[1]-xlim[2])
      border_to_center <- as.data.frame(matrix(NA, 2, 2))
      border_to_center[,1] <- c(xlim/10, 0)
      border_to_center[,2] <- c(0,0)
      scaleBAR <-
        round(matrix(unlist(
          rgdal::project(
            as.matrix(border_to_center),
            proj = as.character(projEAC),
            inv = F
          )
        ), ncol = 2) / 1000, 0)[1, 1]
    }else{
      scaleBAR <- Resolution/1000
    }
    raster::scalebar(scaleBAR, type="bar", below="kilometers", cex=2.5)
    
    mtext(NamesSp, side=3, cex=3, line=3)
    if(any(colnames(DATA)=="higher.tax.rank")) mtext(DATA[which(DATA[,3]==NamesSp),"higher.tax.rank"][1], side=3, cex=3, line=0.4)
    
    if(add.legend) {
      graphics::par(mar=c(1,1,1,1), xpd=T)
      plot(1:10, 1:10, type="n", bty='n', xaxt='n', yaxt='n')
      if(is.null(protec.areas)){
        legend(1,10,  c(paste("EOO=", ifelse(!is.na(Results["EOO",1]), round(as.numeric(Results["EOO",1]),1), NA), "km2"),
                        paste("AOO (grid res.",Cell_size_AOO,"km)=", format(Results["AOO",1], scientific = 5),"km2"),
                        paste("Number of unique occurrences=", Results["Nbe_unique_occ.",1]),
                        paste("Number of sub-populations (radius",Resol_sub_pop,"km)=", Results["Nbe_subPop",1]),
                        paste("Number of locations (grid res.:",round(Resolution/1000,1)," km)","=", Results["Nbe_loc",1]),
                        paste("IUCN category according to criterion B:", Results["Category_CriteriaB",1])), cex=3.5, bg = grDevices::grey(0.9))
      }
      if(!is.null(protec.areas)){
        legend(1,10,  c(paste("EOO=", ifelse(!is.na(Results["EOO",1]), round(as.numeric(Results["EOO",1]),1), NA),"km2"),
                        paste("AOO (grid res.",Cell_size_AOO,"km)=", format(Results["AOO",1], scientific = 5),"km2"),
                        paste("Number of unique occurrences=", Results["Nbe_unique_occ.",1]),
                        paste("Number of sub-populations (radius",Resol_sub_pop,"km)=",Results["Nbe_subPop",1]),
                        paste("Number of locations (grid res.:",round(Resolution/1000,1)," km)","=", Results["Nbe_loc",1]),
                        paste("Number of occupied protected areas=", Results["Nbe_loc_PA",1]),
                        paste("IUCN category according to criterion B:", Results["Category_CriteriaB",1]),
                        paste("Proportion of occurences within protected areas"), Results["Ratio_occ_within_PA",1]), cex=3.5, bg = grDevices::grey(0.9))
      }
      graphics::par(mar=c(4,1,1,1))
      sp::plot(full_poly_borders, lty=1, lwd=1,axes=FALSE)
      graphics::points(XY[,1],XY[,2], pch=8, cex=2, col="red") 
    }
    
    if (any(colnames(DATA) == "coly") & add.legend) {
      graphics::par(
        mar = c(12, 6, 1, 2),
        las = 2,
        yaxs = "r",
        xpd = FALSE
      )
      subdata <- DATA[which(DATA[, "tax"] == NamesSp), "coly"]
      if ((sum(subdata, na.rm = T)) > 0) {
        plot(
          table(subdata),
          col = "grey",
          ylab = " ",
          xlab = " ",
          cex.lab = 4,
          cex.axis = 4,
          axes = F
        )
        graphics::axis(
          1,
          outer = FALSE,
          cex.axis = 3,
          tick = FALSE,
          line = 1.5
        )  #pos=min(range(XY[,2]))-2)
        graphics::axis(
          1,
          outer = FALSE,
          labels = FALSE,
          cex.axis = 3,
          tick = TRUE,
          line = 0
        )  #pos=min(range(XY[,2]))-2)
        graphics::axis(
          2,
          outer = FALSE,
          cex.axis = 3,
          tick = FALSE,
          line = 1.5
        )  #pos=min(range(XY[,2]))-2)
        graphics::axis(
          2,
          outer = FALSE,
          labels = FALSE,
          cex.axis = 3,
          tick = TRUE,
          line = 0
        )  #pos=min(range(XY[,2]))-2)
      }
    }
    
    if(!map_pdf) grDevices::dev.off()
  } # end draw map
  
  
  
  
  if(write_shp) {
    
    dir.create(file.path(paste(getwd(), "/shapesIUCN", sep = "")), showWarnings = FALSE)
    
    if (!is.null(p1)) {
      
      if (length(list.files(paste(getwd(), "/shapesIUCN", sep = ""))) > 0) {
        
        if (length(grep(gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_EOO_poly")), unique(sub(
          "....$", '', list.files(paste(getwd(), "/shapesIUCN", sep = ""))
        )))) > 0) {
          
          FILES <-
            list.files(paste(getwd(), "/shapesIUCN", sep = ""), full.names = TRUE)
          file.remove(FILES[grep(gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_EOO_poly")), FILES)])
          
        }
        
      }
      
      
      # NAME <- names(p1)
      # p1@polygons[[1]]@ID <- "1"
      # 
      # ConvexHulls_poly_dataframe <-
      #   sp::SpatialPolygonsDataFrame(p1, data = as.data.frame(names(p1)))
      # 
      # colnames(ConvexHulls_poly_dataframe@data) <-
      #   paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = '')
      # 
      # rgdal::writeOGR(
      #   ConvexHulls_poly_dataframe,
      #   "shapesIUCN",
      #   gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_EOO_poly")),
      #   driver = "ESRI Shapefile"
      # )
      
      
      
      
      rgdal::writeOGR(
        p1,
        "shapesIUCN",
        gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_EOO_poly")),
        driver = "ESRI Shapefile"
      )
      
    }
    
    if(SubPop) {
      
      if (length(list.files(paste(getwd(), "/shapesIUCN", sep = ""))) > 0) {
        if (length(grep(gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_subpop_poly")), unique(sub(
          "....$", '', list.files(paste(getwd(), "/shapesIUCN", sep = ""))
        )))) > 0) {
          FILES <-
            list.files(paste(getwd(), "/shapesIUCN", sep = ""), full.names = TRUE)
          file.remove(FILES[grep(gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_subpop_poly")), FILES)])
          
        }
      }
      
      
      df <-
        data.frame(FID = row.names(SubPopPoly),
                   row.names = row.names(SubPopPoly))
      SubPopPoly <- SpatialPolygonsDataFrame(SubPopPoly, data=df)
      
      rgdal::writeOGR(
        obj = SubPopPoly,
        dsn = "shapesIUCN",
        layer = gsub(pattern = " ", replacement = "_" , x = paste0(NamesSp, "_subpop_poly")),
        driver = "ESRI Shapefile"
      )
      
      # NAME <- names(SubPopPoly)
      # SubPopPoly@polygons[[1]]@ID <- "1"
      # ConvexHulls_poly_dataframe <- sp::SpatialPolygonsDataFrame(SubPopPoly, data=as.data.frame(names(SubPopPoly)))
      # colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = '')
      # rgdal::writeOGR(ConvexHulls_poly_dataframe,"shapesIUCN",paste(NamesSp,"_subpop_poly", sep=""),driver="ESRI Shapefile")
      
    }
    
  }
  
  if (SubPop) {
    
    OUTPUT <- list(Results, p1, SubPopPoly)
    names(OUTPUT) <-
      c("Results", "spatialPoly_EOO", "spatialPoly_subpop")
    
  } else{
    
    OUTPUT <- list(Results, p1)
    
    names(OUTPUT) <- c("Results", "spatialPoly_EOO")
    
  }
  
  return(OUTPUT)
}