#' @title Extent of Occurrences multi-taxa computation
#' 
#' @description Compute extent of occurrences (EOO) for multiple taxa in square kilometers
#' using \code{"geosphere"} package and provide
#' \code{SpatialPolygons} used for EOO computation
#' 
#' 
#' @details 
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not
#' matter}
#' 
#' \tabular{ccc}{ [,1] \tab ddlat \tab numeric, latitude (in decimal
#' degrees)\cr [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#' [,3] \tab tax \tab character or factor, taxa names\cr }
#' 
#' \strong{Important notes:}
#' 
#' EOO will only be computed if there is at least three unique occurrences
#' unless \code{method.less.than3} is put to "arbitrary". In that specific
#' case, EOO for species with two unique occurrences will be equal to
#' Dist*Dist*0.1 where Dist is the distance in kilometers separating the two
#' points.
#' 
#' For the very specific (and infrequent) case where all occurrences are
#' localized on a straight line (in which case EOO would be null), 
#' 'noises' are added to coordinates, using the `jitter` function, 
#' see \href{https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/jitter}. 
#' There is a warning when this happens. This means that EOO value will not be constant 
#' across multiple estimation (although the variation should be small)
#' 
#' \strong{Limitation}\cr
#' 
#' For a species whose occurrences span more than 180 degrees, EOO is not
#' computed. This is the case for example for species whose distribution span
#' the 180th meridian.
#' 
#' @param XY \code{dataframe} see Details
#' @param exclude.area a logical, if TRUE, areas outside of \code{country_map}
#' are cropped of \code{SpatialPolygons} used for calculating EOO. By default
#' is FALSE
#' @param country_map a \code{SpatialPolygonsDataFrame} or
#' \code{SpatialPolygons} showing for example countries or continent borders.
#' This shapefile will be used for cropping the \code{SpatialPolygons}l if
#' exclude.area is TRUE
#' @param export_shp a logical, whether shapefiles should be exported or not,
#' see Value. By default is FALSE
#' @param driver_shp a string, define the driver for exporting shapefiles, 
#' by default "ESRI Shapefile". See \code{\link[sf]{st_write}}
#' @param write_shp a logical, if TRUE, export \code{SpatialPolygons} used for
#' EOO computation as ESRI shapefiles in the working directory. By default is
#' FALSE
#' @param alpha a numeric, if \code{method.range} is "alpha.hull", value of
#' alpha of the alpha hull, see \code{\link[alphahull]{ahull}}. By default is 1
#' @param buff.alpha a numeric, if \code{method.range} is "alpha.hull", define
#' the buffer in decimal degree added to alpha hull. By default is 0.1
#' @param method.range a character string, "convex.hull" or "alpha.hull". By
#' default is "convex.hull"
#' @param Name_Sp a character string, if \code{XY} is for one taxon and field
#' containing taxon names is not provided, this item provide taxon name. By
#' default is "Species1"
#' @param method.less.than3 a character string. If equal to "arbitrary", will
#' give a value to species with two unique occurrences, see Details. By default
#' is "not comp"
#' @param write_results a logical. If TRUE, results will be exported in the
#' working environment as a csv file. By default is TRUE
#' @param file.name a character string. Name file for exported results in csv
#' file. By default is "EOO.results"
#' @param parallel a logical. Whether running in parallel. By default, it is
#' FALSE
#' @param NbeCores an integer. Register the number of cores for parallel
#' execution. By default, it is 2
#' @param show_progress logical. Whether a progress bar should displayed. TRUE by default
#' @param proj_type character string or numeric or object of CRS class, by default is "cea"
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'#'
#'
#' @return If \code{export_shp} is FALSE, a \code{dataframe} with one field
#' containing EOO in square kilometers.  \code{NA} is given when EOO could not
#' be computed because there is less than three unique occurrences (or two if
#' \code{method.less.than3} is put to "arbitrary").
#' 
#' If \code{export_shp} is TRUE, a \code{list} with: \enumerate{ \item EOO in
#' square kilometers \item \code{SpatialPolygons} used for EOO computation}
#' 
#' @author Gilles Dauby
#' 
#' \email{gildauby@@gmail.com}
#' @seealso \code{\link[alphahull]{ahull}}
#' 
#' \url{https://github.com/azizka/speciesgeocodeR}
#' 
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges,
#' Journal of Applied Ecology, 49 1-9
#' 
#' @examples
#' 
#' data(dataset.ex)
#' data(land)
#' \dontrun{
#' EOO <- EOO.computing(dataset.ex)
#' 
#' ## This exclude areas outside of land (i.e. ocean) for EOO computation
#' EOO <- EOO.computing(dataset.ex, 
#'   exclude.area=TRUE, country_map=land)
#' }
#' 
#' @import sf
#' 
#' @importFrom rnaturalearth ne_countries
#' @importFrom rgeos gBuffer
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export
EOO.computing <- function(XY,
                          exclude.area = FALSE,
                          country_map = NULL,
                          export_shp = FALSE,
                          driver_shp = "ESRI Shapefile",
                          write_shp = FALSE,
                          alpha = 1,
                          buff.alpha = 0.1,
                          method.range = "convex.hull",
                          Name_Sp = "species1",
                          method.less.than3 = "not comp",
                          write_results = TRUE,
                          file.name = "EOO.results",
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          proj_type = "cea",
                          mode = "spheroid"
) {
  
  
  list_data <- coord.check(XY = XY)
  

  ### Getting by default land map if poly_borders is not provided
  if (is.null(country_map)) {
    
    country_map <-
      rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    
  } else {
    
    if(any(grepl('sf', class(country_map))))
      country_map <- 
        as(country_map, "Spatial")
    
    country_map <-
      suppressWarnings(rgeos::gBuffer(country_map, byid = TRUE, width = 0))
    
    country_map <- 
      as(country_map, "sf")
  }
  
  # if (buff_width > 80)
  #   stop("buff_width has unrealistic value")
  

  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
  
  
  if (is.null(names(list_data))) {
    names_ <-
      rep(Name_Sp, length(list_data))
  } else {
    names_ <- names(list_data)
  }
  
  x <- NULL
  if (show_progress) {
    pb <-
      txtProgressBar(min = 0,
                            max = length(list_data),
                            style = 3)
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }
  
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      # source("./R/EOO.comp.R")
      # source("./R/Convex.Hull.Poly.R")
      # source("./R/proj_crs.R")
      # source("./R/ahull_to_SPLDF.R")
      # source("./R/coord.check.R")
      # library(sf)
      # library(sp)
      
      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      
      res <-
        EOO.comp(
          XY = list_data[[x]],
          exclude.area = exclude.area,
          country_map = country_map,
          Name_Sp = names_[x],
          method.range = method.range,
          alpha = alpha,
          buff.alpha = buff.alpha,
          method.less.than3 = method.less.than3, 
          mode = mode, 
          proj_type = proj_type
        )
      
      names(res)[1] <-
        paste0(names(res)[1], "_" , x)
      if (length(res) > 1)
        names(res)[2] <-
        paste0(names(res)[2], "_" , x)
      
      res
      
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  Results_short <-
    data.frame(EOO = unlist(output[grep("EOO", names(output))]))
  row.names(Results_short) <- names_
  
  if (length(output) == 1)
    names(output) <- Name_Sp
  
  
  if(export_shp) {
    
    output_spatial <- output[grep("spatial", names(output))]
    output_spatial <- output_spatial[!is.na(output_spatial)]
    
    id_spatial <-
      as.numeric(unlist(lapply(strsplit(
        names(output_spatial), "_"
      ), function(x)
        x[[2]])))
    
    if(length(output_spatial) > 1) {
      output_spatial <- 
        do.call("rbind", output_spatial)
    } else {
      output_spatial <- 
        output_spatial[[1]]
    }
    
    output_spatial <- 
      st_as_sf(data.frame(output_spatial[, -which(colnames(output_spatial) == 'a')], 
                          taxa = names_[id_spatial]))
    
  }
  
  if(write_shp) {
    
    message("Writing EOO shapefiles in shapesIUCN directory")
    
    dir.create(file.path(paste(getwd(), "/shapesIUCN", sep = "")), showWarnings = FALSE)
    

    
    # output_spatial <- 
    #   output_spatial[unlist(lapply(output_spatial, function(x) !is.vector(x)))]

    
    # exi_files <- 
    #   list.files(paste(getwd(), "/shapesIUCN", sep = ""))
    
    sf::write_sf(output_spatial,
                 dsn = "shapesIUCN",
                 layer = paste("EOO_poly", sep = ""),
                 driver = driver_shp,
                 overwrite = TRUE)
    
    # for (i in 1:length(output_spatial)) {
    #   
    #   NAME <- names_[id_spatial[i]]
    #   NAME <- gsub(" ", "_", NAME)
    #   
    #   sf::write_sf(output_spatial[[i]],
    #                dsn = "shapesIUCN",
    #                layer = paste(NAME, "_EOO_poly", sep = ""),
    #                driver = driver_shp,
    #                overwrite = TRUE)
    #   
    #   # output_spatial[[i]]@polygons[[1]]@ID <- "1"
    #   # ConvexHulls_poly_dataframe <-
    #   #   sp::SpatialPolygonsDataFrame(output_spatial[[i]], data = as.data.frame(names(output_spatial[[i]])))
    #   # colnames(ConvexHulls_poly_dataframe@data) <-
    #   #   paste(substr(names_[id_spatial[i]], 0, 3), collapse = '')
    #   # rgdal::writeOGR(
    #   #   ConvexHulls_poly_dataframe,
    #   #   "shapesIUCN",
    #   #   paste(names_[id_spatial[i]], "_EOO_poly", sep = ""),
    #   #   driver = "ESRI Shapefile",
    #   #   overwrite_layer = TRUE
    #   # )
    # }
  }
  
  if (write_results)
    write.csv(Results_short, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (!export_shp)
    output <- Results_short

  if (export_shp)
    output <- list(results = Results_short,
                   spatial = output_spatial)
  
  output
}
