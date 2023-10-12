#' @title Mapping in grid cell results of IUCN.eval
#' 
#' @description Provides four maps showing in grid cells of a given resolution : number of
#' records, species richness, number of threatened species (CR+EN+VU) and
#' proportion of threatened species. Based on \code{\link[fields]{quilt.plot}}.
#' 
#' @details
#' \strong{Input} \code{Occurrences} as a \code{dataframe} should have the
#' following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not
#' matter}
#' 
#' \tabular{ccc}{ [,1] \tab ddlat \tab numeric, latitude (in decimal
#' degrees)\cr [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#' [,3] \tab tax \tab character or factor, taxa names\cr }
#' 
#' @param Results The default output of \code{\link{IUCN.eval}} applied to
#' multiple species
#' @param Occurrences A \code{dataframe}, see Details
#' @param country_map A \code{SpatialPolygonsDataFrame} or
#' \code{SpatialPolygons} showing for example countries or continent borders
#' @param Resol numeric , resolution in decimal degrees
#' @param threshold numeric, only grid cells with at least this number of
#' records will be shown
#' @param LatMin numeric, minimum latitude for the map
#' @param LatMax numeric, maximum latitude for the map
#' @param LongMin numeric, minimum longitude for the map
#' @param LongMax numeric, maximum longitude for the map
#' @param export_map logical, if TRUE, four maps in png will be created in the
#' working directory if FALSE, maps will be displayed in the R session
#' @param file_name character string. Name of the file
#' @param export_data logical. If TRUE, a \code{dataframe} containing all
#' information on the grid cell mapped is exported
#' @return Produce four maps either in the R session (if \code{export_map} is
#' FALSE) or in png format in the working directory (if \code{export_map} is
#' TRUE)
#' 
#' If \code{export_data} is TRUE
#' 
#' \strong{Output} \tabular{cccccccc}{ [,1] \tab X \tab numeric, x coordinates
#' of cell [,2] \tab Y \tab numeric, y coordinates of cell [,3] \tab meanLat
#' \tab numeric, mean latitude of occurrences within cell\cr [,4] \tab meanLat
#' \tab numeric, mean longitude of occurrences within cell\cr [,5] \tab NbeRec
#' \tab numeric, Number of records\cr [,6] \tab NbeEsp \tab numeric, Number of
#' species\cr [,7] \tab NbeThreatened \tab numeric, Number of threatened
#' species\cr [,8] \tab PropThreatened \tab numeric, Proportion of threatened
#' species\cr }
#' @author Gilles Dauby
#' @seealso package fields function quilt.plot
#' @examples
#' 
#' \dontrun{
#' data(land)
#' data(Malagasy_amphibian)
#' Results <- IUCN.eval(Malagasy_amphibian, DrawMap=FALSE, country_map=land, SubPop=FALSE)
#' ### This should run for 3 to 6 minutes depending of the computer.
#' 
#' ### Maps covering the whole dataset with a minimum of 5 records in each cell
#' map.res(Results=Results, Occurrences=Malagasy_amphibian, country_map=land, 
#' export_map=FALSE, threshold=5)
#' 
#' ## Maps focusing on Madagascar with a minimum of 5 records in each cell
#' map.res(Results=Results, Occurrences=Malagasy_amphibian, country_map=land, export_map=FALSE, 
#' 	threshold=5, LatMin=-25,LatMax=-12,LongMin=42, LongMax=52)
#' 
#' ## Maps focusing on Madagascar at half degree resolution with a minimum of 5 records in each cell
#' map.res(Results=Results, Occurrences=Malagasy_amphibian, country_map=land, 
#' export_map=FALSE,Resol=0.5, 
#' 	threshold=5, LatMin=-25,LatMax=-12,LongMin=42, LongMax=52)
#' 
#' ## Maps have been exported in the directory IUCN__results_map
#' map.res(Results=Results, Occurrences=Malagasy_amphibian, country_map=land, export_map=TRUE, 
#' 	threshold=5, LatMin=-25,LatMax=-12,LongMin=42, LongMax=52)
#' 
#' ## Install speciesgeocodeR package for an example with their lemurs dataset
#' library(speciesgeocodeR)
#' data("lemurs_in")
#' 
#' Results <- IUCN.eval(lemurs_in, DrawMap=FALSE, country_map=land, SubPop=FALSE)
#' 
#' map.res(Results=Results, Occurrences=lemurs_in, country_map=land, export_map=FALSE, threshold=3, 
#' 	LatMin=-25,LatMax=-12,LongMin=42, LongMax=52, Resol=1)
#' 
#' }
#' 
#' @importFrom fields quilt.plot two.colors image.plot
#' @importFrom rnaturalearth ne_countries
#' 
#' @importFrom grDevices dev.cur dev.off grey pdf png rgb
#' @importFrom graphics axis box layout legend mtext par points
#' @importFrom methods as slot
#' @importFrom utils installed.packages write.csv
#' 
#' @export
#' 
map.res <- function(Results,
                    Occurrences,
                    country_map = NULL,
                    Resol = 1,
                    threshold = 0,
                    LatMin = NULL,
                    LatMax = NULL,
                    LongMin = NULL,
                    LongMax = NULL,
                    export_map = FALSE,
                    file_name = NULL,
                    export_data = FALSE) {
  
  if (nrow(Results) != length(unique(as.character(Occurrences[[3]]))))
    stop("Results and Occurrences input files have different number of species")
  
  if (class(Occurrences) == "spgeoIN") {
    DATA_2 <-
      cbind(Occurrences$species_coordinates, Occurrences$identifier)
    Occurrences <- DATA_2[, c(2, 1, 3)]
    colnames(Occurrences)[1:3] <- c("ddlat", "ddlon", "tax")
  } else{
    colnames(Occurrences)[1:3] <- c("ddlat", "ddlon", "tax")
  }
  
  Results_full <- cbind(rownames(Results), Results)
  colnames(Results_full)[1] <- "tax"
  merged_data_criteriaB <-
    base::merge(Results_full, Occurrences, by.x = "tax", by.y = "tax")
  
  if (is.null(LatMin))
    LatMin = min(merged_data_criteriaB[, "ddlat"])
  if (is.null(LatMax))
    LatMax = max(merged_data_criteriaB[, "ddlat"])
  if (is.null(LongMin))
    LongMin = min(merged_data_criteriaB[, "ddlon"])
  if (is.null(LongMax))
    LongMax = max(merged_data_criteriaB[, "ddlon"])
  
  if (LatMin >= LatMax)
    stop("LatMin must be lower than LatMax")
  if (LongMin >= LongMax)
    stop("LongMin must be lower than LongMax")
  if (LongMin > 180 ||
      LongMin < -180 ||
      LatMin > 180 ||
      LatMin < -180 ||
      LatMax > 180 ||
      LatMax < -180 ||
      LongMax > 180 ||
      LongMax < -180)
    stop("Latitude and longitude must be within [-180; 180] intervall")
  
  EXTENT <- raster::extent(LongMin, LongMax, LatMin, LatMax)
  
  if(is.null(country_map))  {
    
    land <- 
      country_map <- 
      rnaturalearth::ne_countries(scale = 50, returnclass = "sp")
    
    # data('land', package='ConR', envir=environment()) 
    # land <- get("land", envir=environment()) 
    # country_map=land
  }
  
  if (!is.null(country_map)) {
    cropped_country_map <- raster::crop(country_map, EXTENT + 20)
  }
  
  DATA_SF <- merged_data_criteriaB
  sp::coordinates(DATA_SF) <-  ~ ddlon + ddlat
  
  sp::proj4string(DATA_SF) <- sp::CRS(SRS_string='EPSG:4326')
  
  # raster::crs(DATA_SF) <- raster::crs(country_map)
  
  DATA_SF$X <- floor(merged_data_criteriaB[, "ddlon"] / Resol)
  DATA_SF$Y <- floor(merged_data_criteriaB[, "ddlat"] / Resol)
  DATA_SF$Cell <- paste("M", DATA_SF$X, "x", DATA_SF$Y, sep = "")
  DATA_SF@data <-
    cbind(DATA_SF@data, merged_data_criteriaB[, c("ddlat", "ddlon")])
  
  print(paste("Number of cell with at least one occurrence is", length(unique(
    as.character(DATA_SF@data[, "Cell"])
  ))))
  print(
    paste(
      "Number of cell with number of occurrences higher or equal to",
      threshold,
      "is",
      length(which(table(DATA_SF@data[, "Cell"]) > threshold))
    )
  )
  if (length(which(table(DATA_SF@data[, "Cell"]) > threshold)) == 0)
    stop("No cell left after filtering")
  
  counts <-
    by(DATA_SF@data, DATA_SF@data$Cell, function(d)
      c(
        d$X[1],
        d$Y[1],
        mean(d$ddlat),
        mean(d$ddlon),
        .prop_threat(d, threshold)
      ))
  
  threatened_rec <- matrix(unlist(counts), nrow=8)
  rownames(threatened_rec) <- c("X", "Y","meanLat","meanLong", "NbeRec", "NbeEsp", "NbeThreatened", "PropThreatened")
  colnames(threatened_rec) <- names(counts)
  if(ncol(threatened_rec)<2) stop("All records are within one grid cell, decrease the resolution to have a relevant map")
  
  threatened_rec_cut <- as.data.frame(threatened_rec[, which(threatened_rec["Y",]*Resol> (EXTENT+Resol)[3])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  threatened_rec_cut <- as.data.frame(threatened_rec_cut[,which(threatened_rec_cut["Y",]*Resol< (EXTENT+Resol)[4])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  threatened_rec_cut <- as.data.frame(threatened_rec_cut[,which(threatened_rec_cut["X",]*Resol> (EXTENT+Resol)[1])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  threatened_rec_cut <- as.data.frame(threatened_rec_cut[,which(threatened_rec_cut["X",]*Resol< (EXTENT+Resol)[2])])
  if(ncol(threatened_rec_cut)<2) stop("All records are within one grid cell, decrease the resolution/threshold or modify extent")
  
  COORD <- t(rbind(threatened_rec_cut[1, ], threatened_rec_cut[2, ]))
  
  grid.list <-
    list(x = (Resol / 2 + Resol * seq(range(COORD[, 1])[1], range(COORD[, 1])[2], 1)),
         y = (Resol / 2 + Resol * seq(range(COORD[, 2])[1], range(COORD[, 2])[2], 1)))
  
  
  SelectedCells <- which(threatened_rec_cut["NbeRec", ] >= threshold)
  
  if (export_map) {
    DIRECTORY_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
    dir.create(file.path(paste(
      getwd(), paste("/", DIRECTORY_NAME, "results_map", sep = ""), sep = ""
    )), showWarnings = FALSE)
  }
  
  BG <- 'white'
  Border <- 'black'
  COlor <- 'grey97' # rgb(0.1, 0.3, 0.1, alpha=0.1)
  
  if (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  if (grDevices::dev.cur() > 1) {
    grDevices::dev.off()
  }
  
  if (!export_map)
    graphics::par(mfrow = c(2, 2))
  if (export_map)
    graphics::par(mfrow = c(1, 1))
  if (export_map)
    grDevices::png(
      paste(file.path(paste(
        getwd(), paste("/", DIRECTORY_NAME, "results_map", sep = ""), sep = ""
      )), "/", "number_of_records", ".png", sep = ""),
      width = 20,
      height = 20,
      units = "cm",
      res = 150
    )
  
  ################################### Number of records
  coltab <-
    fields::two.colors(256,
                       start = "lightblue",
                       end = "red",
                       middle = "yellow")
  
  if (!export_map)
    graphics::par(mar = c(4, 1, 1, 4),
                  las = 1,
                  omi = c(0.5, 1, 0.5, 0.3))
  if (export_map)
    graphics::par(
      mar = c(2, 2, 1, 5),
      las = 1,
      omi = c(0.3, 0.4, 0.3, 0.1)
    )
  
  if (export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1
    )
  
  if (!export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1,
      xaxt = 'n'
    )
  
  graphics::mtext(text = "Number of records", side = 3, cex = 1)
  VALUES <- as.numeric(threatened_rec_cut["NbeRec", SelectedCells])
  
  fields::quilt.plot(
    COORD[SelectedCells, 1] * Resol + Resol / 2 ,
    COORD[SelectedCells, 2] * Resol + Resol / 2,
    VALUES  ,
    grid = grid.list ,
    cex.axis = 1,
    cex.lab = 1,
    add.legend = FALSE,
    col = coltab,
    add = T
  )
  sp::plot(cropped_country_map, add = T)
  if (min(VALUES) == max(VALUES))
    Range <- c(min(VALUES), min(VALUES) + 1)
  if (min(VALUES) != max(VALUES))
    Range <- range(VALUES)
  fields::image.plot(
    zlim = Range,
    legend.only = TRUE,
    col = coltab,
    legend.shrink = 1 ,
    legend.width = 1,
    cex.lab = 2,
    axis.args = list(
      cex.axis = 1,
      col.lab = Border,
      col.axis = Border
    )
  )
  graphics::box()
  if (export_map)
    grDevices::dev.off()
  
  ################################### Number of species
  coltab <- fields::two.colors(256, start="cyan", end="darkorange4", middle="gold")
  
  if (export_map)
    grDevices::png(
      paste(file.path(paste(
        getwd(), paste("/", DIRECTORY_NAME, "results_map", sep = ""), sep = ""
      )), "/", "species_richness", ".png", sep = ""),
      width = 20,
      height = 20,
      units = "cm",
      res = 150
    )
  if (export_map)
    graphics::par(
      mar = c(2, 2, 1, 5),
      las = 1,
      omi = c(0.3, 0.4, 0.3, 0.1)
    )
  if (export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1
    )
  
  if (!export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1,
      xaxt = 'n',
      yaxt = 'n'
    )
  graphics::mtext(text = "Species richness", side = 3, cex = 1)
  VALUES <- as.numeric(threatened_rec_cut["NbeEsp", SelectedCells])
  
  fields::quilt.plot(
    COORD[SelectedCells, 1] * Resol + Resol / 2 ,
    COORD[SelectedCells, 2] * Resol + Resol / 2,
    VALUES  ,
    grid = grid.list ,
    cex.axis = 1,
    cex.lab = 1,
    add.legend = FALSE,
    col = coltab,
    add = T
  )
  
  sp::plot(cropped_country_map, add = T)
  if (min(VALUES) == max(VALUES))
    Range <- c(min(VALUES), min(VALUES) + 1)
  if (min(VALUES) != max(VALUES))
    Range <- range(VALUES)
  fields::image.plot(
    zlim = Range,
    legend.only = TRUE,
    col = coltab,
    legend.shrink = 1 ,
    legend.width = 1,
    cex.lab = 2,
    axis.args = list(
      cex.axis = 1,
      col.lab = Border,
      col.axis = Border
    )
  )
  graphics::box()
  if (export_map)
    grDevices::dev.off()
  
  ################################### Number of threatened species
  coltab<- fields::two.colors(256, start="slategray2", end="deeppink4", middle="burlywood2")
  
  if (export_map)
    grDevices::png(
      paste(file.path(paste(
        getwd(),
        paste("/", DIRECTORY_NAME, "results_map", sep = ""),
        sep = ""
      )), "/", "number_threatened_sp", ".png", sep = ""),
      width = 20,
      height = 20,
      units = "cm",
      res = 150
    )
  if (export_map)
    graphics::par(
      mar = c(2, 2, 1, 5),
      las = 1,
      omi = c(0.3, 0.4, 0.3, 0.1)
    )
  
  if (export_map) 
    sp::plot(cropped_country_map, axes=T, lty=1,border=Border, col=COlor, xlim=c(LongMin,LongMax), ylim=c(LatMin,LatMax), cex.axis=1,lwd=1)
  
  if (!export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1
    )
  if (export_map)
    graphics::par(
      mar = c(2, 2, 1, 5),
      las = 1,
      omi = c(0.3, 0.4, 0.3, 0.1)
    )
  graphics::mtext(text = "Number of threatened species", side = 3, cex =
                    1)
  VALUES <-
    as.numeric(threatened_rec_cut["NbeThreatened", SelectedCells])
  
  fields::quilt.plot(COORD[SelectedCells,1]*Resol+Resol/2 , COORD[SelectedCells,2]*Resol+Resol/2, VALUES  ,grid=grid.list , cex.axis=1, 
                     cex.lab=1, add.legend=FALSE, col=coltab, add=T)
  sp::plot(cropped_country_map, add=T)
  if(min(VALUES)==max(VALUES)) Range <- c(min(VALUES), min(VALUES)+1)
  if(min(VALUES)!=max(VALUES)) Range <- range(VALUES)
  fields::image.plot(zlim=Range,legend.only=TRUE, col=coltab, legend.shrink = 1 ,
                     legend.width=1, cex.lab=2, axis.args=list(cex.axis = 1, col.lab = Border, col.axis = Border))
  graphics::box()
  if (export_map) grDevices::dev.off()
  
  ################################### Proportion of threatened species
  coltab <- fields::two.colors(256, start="darkslategray1", end="hotpink2", middle="khaki")
  
  if (export_map)
    grDevices::png(
      paste(
        file.path(paste(
          getwd(),
          paste("/", DIRECTORY_NAME, "results_map", sep = ""),
          sep = ""
        )),
        "/",
        "proportion_threatened_sp",
        ".png",
        sep = ""
      ),
      width = 20,
      height = 20,
      units = "cm",
      res = 150
    )
  if (export_map)
    graphics::par(
      mar = c(2, 2, 1, 5),
      las = 1,
      omi = c(0.3, 0.4, 0.3, 0.1)
    )
  
  if (!export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1,
      yaxt = 'n'
    )
  if (export_map)
    sp::plot(
      cropped_country_map,
      axes = T,
      lty = 1,
      border = Border,
      col = COlor,
      xlim = c(LongMin, LongMax),
      ylim = c(LatMin, LatMax),
      cex.axis = 1,
      lwd = 1
    )
  
  mtext(text="Proportion of threatened species",side=3, cex=1)
  VALUES <- as.numeric(threatened_rec_cut["PropThreatened",SelectedCells])
  
  fields::quilt.plot(COORD[SelectedCells,1]*Resol+Resol/2 , COORD[SelectedCells,2]*Resol+Resol/2, VALUES  ,grid=grid.list , cex.axis=1, 
                     cex.lab=1, add.legend=FALSE, col=coltab, add=T)
  sp::plot(cropped_country_map, add=T)
  if(min(VALUES)==max(VALUES)) Range <- c(min(VALUES), min(VALUES)+1)
  if(min(VALUES)!=max(VALUES)) Range <- range(VALUES)
  fields::image.plot(zlim=Range,legend.only=TRUE, col=coltab, legend.shrink = 1 ,
                     legend.width=1, cex.lab=2, axis.args=list(cex.axis = 1, col.lab = Border, col.axis = Border))
  graphics::box()
  
  if (export_map) grDevices::dev.off()
  
  if(export_data) return(threatened_rec_cut)
  
}