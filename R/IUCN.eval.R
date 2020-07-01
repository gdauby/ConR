
#' @title Preliminary conservation status assessment following IUCN Criterion B
#' 
#' @description Given a dataframe of georeferenced occurrences of one, or more, taxa, this
#' function provide statistics values (Extent of Occurrence, Area of Occupancy,
#' number of locations, number of subpopulations) and provide a preliminary
#' conservation status following Criterion B of IUCN.  A graphical map output
#' is also available.
#' 
#' @details
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not
#' matter}
#' 
#' \tabular{ccccc}{ [,1] \tab ddlat \tab numeric, latitude (in decimal
#' degrees)\cr [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#' [,3] \tab tax \tab character or factor, taxa names\cr [,4] \tab family \tab
#' character, optional field indicating higher taxonomic rank\cr [,5] \tab coly
#' \tab numeric, optional field indicating collection year\cr }
#' 
#' \code{coly} and \code{family} are optional fields
#' 
#' If the optional field named 'family' is provided, indicating higher
#' taxonomic rank, this will be displayed in the title of the map if
#' \code{DrawMap} is 'TRUE'.
#' 
#' If the optional field named 'coly' is provided, indicating collection year,
#' a sub-graph in the map will be displayed (if \code{DrawMap} and
#' \code{add.legend} are both TRUE) showing a barplot of collection year
#' 
#' \strong{Starting position of the raster used for estimating the Area Of
#' Occupancy}\cr Different starting position of the raster used for estimate
#' the AOO may provide different number of occupied cells. Hence, by default, 4
#' different translations of the raster is done (fixed increment of 1/4
#' resolution north and east) and the minimum number of occupied cells is used
#' for estimating AOO. It is also possible to define a given number of random
#' starting position of the raster using the argument
#' \code{nbe.rep.rast.AOO}\cr
#' 
#' \strong{Estimating number of locations}\cr Locations are estimated by
#' overlaying a grid of a given resolution (see \code{Cell_size_locations} for
#' specifying the resolution). The number of locations is simply the number of
#' occupied locations. Note that the grid position is overlaid in order to
#' minimize the number of locations (several translation of the grid are
#' performed and the one providing the minimum number of occupied cells is
#' provided).
#' 
#' \strong{Taking into account protected area for estimating the number of
#' locations}\cr A location is defined by the IUCN as a "geographically or
#' ecologically distinct area in which a single threatening event can affect
#' all individuals of the taxon". A simple way to include threat level is to
#' rely on a map of protected areas and assume that populations within and
#' outside protected areas are under different threat level.\cr
#' 
#' If a map of protected area is provided, this one is used for estimating the
#' number of locations by the following procedure:\cr - if
#' \code{method_protected_area} is "no_more_than_one", all occurrences within a
#' given protected area will be considered as one location. Occurrences outside
#' protected area will be used for estimating the number of locations using
#' overlaying grid as descrived above. See the vignette for illustration. \cr -
#' if \code{method_protected_area} is NOT "no_more_than_one", number of
#' locations will be estimated by the overlaying grid as described above, but
#' by considering differently occurrences outside and inside protected area.
#' See the vignette for illustration. \cr
#' 
#' The protected areas layers should be given as as
#' \code{SpatialPolygonsDataFrame} in \code{protec.areas}. The
#' \code{ID_shape_PA} should also be given and should represent the unique ID
#' of each protected area in the provided shapefile. This can be checked by the
#' following code:
#' 
#' \code{colnames(ProtectedAreas@data)} Where ProtectedAreas is the name of
#' your shapefile.
#' 
#' \strong{Limitation in the estimations of EOO}\cr
#' 
#' For a species whose occurrences span more than 180 degrees, EOO is not
#' computed. This is the case for example for species whose distribution span
#' the 180th meridian.
#' 
#' @param DATA a \code{dataframe} or an object of class \code{spgeoIN} see
#' \url{https://github.com/azizka/speciesgeocodeR}. See Details
#' @param country_map a \code{SpatialPolygonsDataFrame} or
#' \code{SpatialPolygons} showing for example countries or continent borders.
#' This shapefile will be used for cropping the \code{SpatialPolygons} used for
#' EOO computation if \code{exclude.area} is TRUE. By default, it is
#' \code{land}
#' @param Cell_size_AOO a numeric, value indicating the grid size in kilometers
#' used for estimating Area of Occupancy.  By default, equal to 2
#' @param Cell_size_locations a numeric, value indicating the grid size in
#' kilometers used for estimating the number of location. By default, equal to
#' 10
#' @param Resol_sub_pop a numeric, value indicating the radius size in
#' kilometers used for estimating the number of sub-population. By default,
#' equal to 5
#' @param DrawMap a logical, if TRUE a map is produced for each species in png
#' format, unless map_pdf is TRUE. By default, it is FALSE
#' @param add.legend a logical, if TRUE a legend and a submap showing
#' distribution in 'country_map' are displayed for each map. By default, it is
#' TRUE
#' @param method_locations a character string, indicating the method used for
#' estimating the number of locations.  "fixed_grid" or "sliding scale". See
#' details. By default, it is "fixed_grid"
#' @param Rel_cell_size a numeric, if \code{method_locations="sliding scale"},
#' \code{Cell_size_locations} is ignored and the resolution is given by the
#' maximum distance separating two occurrences multiplied by
#' \code{Rel_cell_size}. By default, it is 0.05
#' @param file_name a character string. Name of the file. By default, it is
#' "IUCN_"
#' @param export_shp a logical, if TRUE, shapefiles of \code{SpatialPolygons}
#' used for EOO computation are exported. By default, it is FALSE
#' @param write_shp a logical, if TRUE, shapefiles of \code{SpatialPolygons}
#' used for EOO computation are written as ESRI shapefiles in a sub-directory
#' in the working directory. By default, it is FALSE
#' @param write_results a logical, if TRUE, results are exported in a file
#' which can csv or excel, see write_file_option. By default, it is TRUE
#' @param write_file_option a character, if "excel", results are exported in
#' excel file, if "csv", results are exported in csv. By default, it is "excel"
#' @param protec.areas a \code{SpatialPolygonsDataFrame}, shapefile with
#' protected areas.  If provided, this will be taken into account for
#' calculating number of location (see Details and
#' \code{method_protected_area}).  By default, no shapefile is provided
#' @param ID_shape_PA a character string, indicating the field name of
#' \code{protec.areas} with ID of the \code{SpatialPolygonsDataFrame} of
#' protected areas
#' @param map_pdf a logical, if TRUE, maps are exported in one pdf file.
#' Otherwise, each species map is exported in png. By default, it is FALSE
#' @param draw.poly.EOO a logical, if TRUE, the polygon used for estimating EOO
#' is drawn. By default, it is TRUE
#' @param exclude.area a logical, if TRUE, areas outside of \code{country_map}
#' are cropped of \code{SpatialPolygons} used for EOO computation. By default,
#' it is FALSE
#' @param method_protected_area a character string. By default is
#' "no_more_than_one"", which means occurrences within protected areas (if
#' provided) will not be taken into account for estimating the number of
#' locations following the grid system, see Details. By default, it is
#' "no_more_than_one"
#' @param SubPop a logical. If TRUE, sub-populations will be estimated. By
#' default, it is TRUE
#' @param buff_width a numeric. For a specific case where all points of a
#' species are on a straight line, EOO is computed by first drawing this
#' straight line and adding a buffer of \code{buff_width} decimal degrees
#' around this line. By default, it is 0.1
#' @param alpha a numeric, if \code{method.range} is "alpha.hull", alpha value
#' for the construction of alpha hull. By default, it is 1
#' @param buff.alpha a numeric, if \code{method.range} is "alpha.hull",
#' indicate the buffer added to the alpha hull in decimal degree. By default,
#' it is 0.1
#' @param method.range a character string, if "convex.hull", EOO is based on a
#' convex hull.  if "alpha.hull", EOO is based on alpha hull of \code{alpha}
#' value. By default, it is "convex.hull"
#' @param nbe.rep.rast.AOO a numeric, indicate the number of raster with random
#' starting position for estimating the AOO. By default, it is NULL but some
#' minimal translation of the raster are still done
#' @param showWarnings a logical. Whether R should report warnings
#' @param parallel a logical. Whether running in parallel. By default, it is
#' FALSE
#' @param NbeCores an integer. Register the number of cores for parallel
#' execution. By default, it is 2
#' @return A \code{dataframe} if 'export_shp' is FALSE. A \code{list} if
#' 'export_shp' is TRUE.
#' 
#' If a \code{list}, three elements are provided: 
#' 
#' \enumerate{
#'   \item a \code{dataframe} with results (see field description below)
#'   \item a list of \code{SpatialPolygons} used for EOO computation
#'   \item a list of \code{SpatialPolygons} used for subpopulations
#' }
#' 
#' The \code{dataframe} has as many rows as taxa and the following fields:
#' 
#' \tabular{cccccc}{ [,1] \tab EOO \tab numeric, EOO (square kilometers)\cr
#' [,2] \tab AOO \tab numeric, AOO (square kilometers)\cr [,3] \tab
#' Nbe_unique_occ. \tab numeric, Number of unique occurrences\cr [,4] \tab
#' Nbe_subPop \tab numeric, Number of subpopulations\cr [,5] \tab Nbe_loc \tab
#' numeric, Number of locations\cr [,6] \tab Category_CriteriaB \tab character,
#' IUCN threat category according to Criterion B\cr [,7] \tab Category_code
#' \tab character, IUCN annotation\cr [,8] \tab Category_AOO \tab character,
#' IUCN threat category according to Criterion B ignoring EOO\cr [,9] \tab
#' Category_EOO \tab character, IUCN threat category according to Criterion B
#' ignoring AOO\cr }
#' @author Gilles Dauby
#' 
#' \email{gildauby@@gmail.com}
#' @seealso \url{https://CRAN.R-project.org/package=biogeo}
#' 
#' \url{https://github.com/azizka/speciesgeocodeR}
#' @references Gaston KJ & Fuller AF, 2009, The sizes of species'geographic
#' ranges, Journal of Applied Ecology, 49 1-9
#' 
#' IUCN Standards and Petitions Subcommitte, 2010, Guidelines for Using the
#' IUCN Red List Categories and Criteria.
#' \url{https://www.iucnredlist.org/resources/categories-and-criteria}
#' 
#' Rivers CM, Bachman SP & Meagher TR, 2010, Subpopulations, locations and
#' fragmentation: applying IUCN red list criteria to herbarium specimen data,
#' Biodiversity Conservation 19:2071-2085
#' @examples
#' 
#' 
#' data(dataset.ex)
#' data(land)
#' \dontrun{
#' Results <- IUCN.eval(dataset.ex, country_map=land)
#' ## A directory has been created in your working directory 
#' and maps for each species has been produced
#' 
#' ### The method for computing locations is a sliding scale:
#' ## the grid resolution will be 0.05*the maximum distance separating occurrences
#' Results <- IUCN.eval(dataset.ex, 
#'                      country_map=land, Cell_size_locations=10,
#'                      Resol_sub_pop = 5, Cell_size_AOO = 4, method_locations="sliding scale")
#' 		}			 
#' \dontrun{
#' ## Install speciesgeocodeR package for an example with their lemurs dataset
#' library(speciesgeocodeR)
#' data("lemurs_in")
#' 
#' Results <- IUCN.eval(lemurs_in, DrawMap=FALSE, country_map=land, SubPop=FALSE)
#' 
#' 
#' 
#' }
#' 
#' @importFrom tibble is_tibble
#' @importFrom writexl write_xlsx
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom rnaturalearth ne_countries
#' @importFrom sf st_transform sf_project st_crs st_cast
#' 
#' @export
IUCN.eval <- function (DATA,
                       country_map = NULL,
                       Cell_size_AOO = 2,
                       Cell_size_locations = 10,
                       Resol_sub_pop = 5,
                       method_locations = "fixed_grid",
                       Rel_cell_size = 0.05,
                       DrawMap = FALSE,
                       add.legend = TRUE,
                       file_name = NULL,
                       export_shp = FALSE,
                       write_shp = FALSE,
                       write_results = TRUE,
                       protec.areas = NULL,
                       map_pdf = FALSE,
                       draw.poly.EOO = TRUE,
                       exclude.area = FALSE,
                       method_protected_area = "no_more_than_one",
                       ID_shape_PA = "WDPA_PID",
                       buff_width = 0.1,
                       SubPop = TRUE,
                       alpha = 1,
                       buff.alpha = 0.1,
                       method.range = "convex.hull",
                       nbe.rep.rast.AOO = 0,
                       showWarnings = TRUE,
                       # verbose=TRUE,
                       write_file_option = "excel",
                       parallel = FALSE,
                       NbeCores = 2) {
  
  
  if(class(DATA)[1]=="spgeoIN") {
    DATA_2 <- cbind(DATA$species_coordinates, DATA$identifier)
    DATA <- DATA_2[,c(2,1,3)]
  }
  colnames(DATA)[1:3] <- c("ddlat","ddlon","tax")
  
  if(tibble::is_tibble(DATA)) DATA <- as.data.frame(DATA)
  
  if (any(is.na(DATA[, 1:2]))) {
    length(which(rowMeans(is.na(DATA[, 1:2])) > 0))
    unique(DATA[which(rowMeans(is.na(DATA[, 1:2])) > 0), 3])
    print(paste(
      "Skipping",
      length(which(rowMeans(is.na(
        DATA[, 1:2]
      )) > 0)) ,
      "occurrences because of missing coordinates for",
      paste(as.character(unique(DATA[which(rowMeans(is.na(DATA[, 1:2])) >
                                             0), 3])), collapse = " AND ")
    ))
    DATA <- DATA[which(!is.na(DATA[, 1])), ]
    DATA <- DATA[which(!is.na(DATA[, 2])), ]
  }
  
  if(is.factor(DATA[,"tax"])) DATA[,"tax"] <- as.character(DATA[,"tax"])
  
  
  if(!is.numeric(DATA[,1]) || !is.numeric(DATA[,2]) ) {
    if(!is.double(DATA[,1]) || !is.double(DATA[,2])) stop("coordinates in DATA should be numeric")
  }
  
  if (any(DATA[, 1] > 180) ||
      any(DATA[, 1] < -180) ||
      any(DATA[, 2] < -180) ||
      any(DATA[, 2] > 180))
    stop("coordinates are outside of expected range")
  if (!is.null(country_map))
    if (!class(country_map) == "SpatialPolygonsDataFrame")
      stop("Country_map should be a spatialpolygondataframe")
  
  if (!is.null(protec.areas)) {
    if (!class(protec.areas) == "SpatialPolygonsDataFrame")
      stop("protec.areas should be a spatialpolygondataframe")
    
    if (!any(colnames(protec.areas@data) %in% ID_shape_PA))
      stop(
        "Check argument ID_shape_PA because selected ID field in the protected area shapefile does not exist"
      )
  }
  
  if (is.null(country_map)) {
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sp")
    
    # data('land', package='ConR', envir=environment())
    # land <- get("land", envir=environment())
    # country_map <- land
  }else{
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
  }
  
  if (!is.null(protec.areas)) {
    if (!sp::identicalCRS(protec.areas, country_map)) {
      
      sp::proj4string(protec.areas) <- 
        sp::CRS(SRS_string = 'EPSG:4326')
      sp::proj4string(country_map) <-
        sp::CRS(SRS_string = 'EPSG:4326')
      
    }
    # raster::crs(protec.areas) <- raster::crs(country_map)
  }
  
  if(length(grep("[?]", DATA[,3]))>0) DATA[,3] <- gsub("[?]", "_", DATA[,3])
  if(length(grep("[/]", DATA[,3]))>0) DATA[,3] <- gsub("[/]", "_", DATA[,3])
  
  #####
  
  list_data <- split(DATA, f = DATA$tax)
  
  if(map_pdf){
    if(!is.null(file_name)) {
      NAME_FILE <- file_name
    }else{
      NAME_FILE <- "IUCN_"
    }
    DIRECTORY_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
    dir.create(file.path(paste(
      getwd(), paste("/", DIRECTORY_NAME, "results_map", sep = ""), sep = ""
    )), showWarnings = FALSE)
    
    pdf(paste(paste(
      getwd(), paste("/", DIRECTORY_NAME, "results_map", sep = ""), sep = ""
    ), "/", "results.pdf", sep = ""),
    width = 25,
    height = 25)
  }
  
  
  
  # library(progressr)
  # library(doFuture)
  # 
  # registerDoFuture()
  # plan(multisession)
  # 
  # handlers("txtprogressbar")
  # 
  # xs <- 1:10
  # 
  # with_progress({
  #   p <- progressor(along = xs) ## create a 5-step progressor
  #   y <- foreach(x = xs) %dopar% {
  #     p()                       ## signal a progression update
  #     Sys.sleep(6.0-x)
  #     sqrt(x)
  #   }
  # })
  # 
  
  if(parallel) {
    
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    # registerDoParallel(NbeCores)
    message('Parallel running with ',
            NbeCores, ' cores')
    `%d%` <- foreach::`%dopar%`
  }else{
    `%d%` <- foreach::`%do%`
  }
  
  x <- NULL
  pb <- 
    utils::txtProgressBar(min = 0, max = length(list_data), style = 3)
  progress <- function(n)
    utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  Results <- 
    foreach::foreach(x = 1:length(list_data),
                     .options.snow = opts) %d% {
                       
                       if (!parallel)
                         utils::setTxtProgressBar(pb, x)
                       
                       res <- 
                         IUCN.comp(DATA = list_data[[x]],
                                    NamesSp = names(list_data)[x], 
                                    DrawMap = DrawMap, 
                                    exclude.area = exclude.area, 
                                    write_shp = write_shp, 
                                    method_protected_area = method_protected_area, 
                                    Cell_size_AOO = Cell_size_AOO, 
                                    Cell_size_locations = Cell_size_locations, 
                                    Resol_sub_pop = Resol_sub_pop, 
                                    method_locations = method_locations, 
                                    Rel_cell_size = Rel_cell_size, 
                                    poly_borders = country_map, 
                                    file_name = file_name, 
                                    buff_width = buff_width, 
                                    map_pdf = map_pdf, 
                                    ID_shape_PA = ID_shape_PA, 
                                    SubPop = SubPop, 
                                    protec.areas = protec.areas, 
                                    add.legend = add.legend,
                                    alpha = alpha, 
                                    buff.alpha = buff.alpha, 
                                    method.range = method.range, 
                                    nbe.rep.rast.AOO = nbe.rep.rast.AOO, 
                                    showWarnings = showWarnings,
                                    MinMax = c(min(list_data[[x]][,2]), max(list_data[[x]][,2]), min(list_data[[x]][,1]), max(list_data[[x]][,1])))
                       
                       
                       res
                     }
  
  if(parallel) snow::stopCluster(cl)
  close(pb)
  
  if(map_pdf) grDevices::dev.off()
  
  Results_short <- lapply(Results, `[`, 1)
  Results_short <-
    lapply(
      Results_short,
      FUN = function(x)
        t(x[[1]])
    )
  
  Results_short <-
    as.data.frame(do.call(rbind, Results_short), stringsAsFactors = FALSE)
  Results_short[, 1:4] <-
    apply(Results_short[, 1:4], MARGIN = 2, as.numeric)
  
  if(write_results) {
    if(!is.null(file_name)) {
      NAME_FILE <- file_name
    }else{
      NAME_FILE <- "IUCN_results"
    }
    
    if(write_file_option=="csv") 
      utils::write.csv(Results_short, paste(getwd(),"/", NAME_FILE, ".csv", sep=""))
    
    if(write_file_option=="excel") {
      Results_short <- data.frame(taxa=rownames(Results_short), Results_short)
      writexl::write_xlsx(Results_short, path = paste(getwd(),"/", NAME_FILE, ".xlsx", sep=""))
    }
    
  }
  
  if(!export_shp) {
    
    Results <- Results_short
    
    if(length(list_data)>5) {
      
      print("Number of species per category")
      print(table(Results[,"Category_CriteriaB"]))
      print("Ratio of species per category")
      print(round(table(Results[,"Category_CriteriaB"])/nrow(Results)*100,1))
      
    }
  }
  
  return(Results)
}
