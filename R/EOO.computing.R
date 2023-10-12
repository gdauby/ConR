#' @title Extent of Occurrences multi-taxa computation
#' 
#' @description
#' Compute extent of occurrences (EOO) for multiple taxa in square kilometers
#' and provide polygons used for EOO computation
#' 
#' 
#' @details
#' **Input** as a `dataframe` should have the following structure:
#' 
#' **It is mandatory to respect field positions, but field names do not
#' matter**
#' 
#' \tabular{lll}{
#' 1 \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#' 2 \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#' 3 \tab tax \tab character or factor, taxa names
#' }
#' 
#' 
#' **Important notes:**
#' 
#' EOO will only be computed if there is at least three unique occurrences
#' unless `method.less.than3` is put to "arbitrary". In that specific
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
#' **Limitation**\cr
#' 
#' For a species whose occurrences span more than 180 degrees, EOO should not be considered. 
#' This is the case for example for species whose distribution span the 180th meridian.
#' 
#' @param XY `dataframe` see Details
#' @param exclude.area a logical, if TRUE, areas outside of `country_map`
#' are cropped of `SpatialPolygons` used for calculating EOO. By default
#' is FALSE
#' @param country_map a `SpatialPolygonsDataFrame` or
#' `SpatialPolygons` showing for example countries or continent borders.
#' This shapefile will be used for cropping the `SpatialPolygons`l if
#' exclude.area is TRUE
#' @param export_shp a logical, whether shapefiles should be exported or not,
#' see Value. By default is FALSE
#' @param driver_shp a string, define the driver for exporting shapefiles, 
#' by default "ESRI Shapefile". See [sf::st_write()]
#' @param write_shp a logical, if TRUE, export `SpatialPolygons` used for
#' EOO computation as ESRI shapefiles in the working directory. By default is
#' FALSE
#' @param alpha a numeric, if `method.range` is "alpha.hull", value of
#' alpha of the alpha hull, see [alphahull::ahull()]. By default is 1
#' @param buff.alpha a numeric, if `method.range` is "alpha.hull", define
#' the buffer in decimal degree added to alpha hull. By default is 0.1
#' @param method.range a character string, "convex.hull" or "alpha.hull". By
#' default is "convex.hull"
#' @param method.less.than3 a character string. If equal to "arbitrary", will
#' give a value to species with two unique occurrences, see Details. By default
#' is "not comp"
#' @param write_results a logical. If TRUE, results will be exported in the working environment as a csv file. By default it is FALSE
#' @param file.name a character string. Name file for exported results in csv file. By default is "EOO.results"
#' @param parallel a logical. Whether running in parallel. By default, it is FALSE
#' @param NbeCores an integer. Register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical. Whether a progress bar should displayed. TRUE by default
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @inheritParams proj_crs
#'
#' @return
#' If `export_shp` is FALSE, a `dataframe` with one field
#' containing EOO in square kilometers. `NA` is given when EOO could not
#' be computed because there is less than three unique occurrences (or two if
#' `method.less.than3` is put to "arbitrary").
#' 
#' If `export_shp` is TRUE, a `list` with: \enumerate{ \item EOO in
#' square kilometers \item `SpatialPolygons` used for EOO computation}
#' 
#' @author Gilles Dauby
#' 
#' \email{gildauby@@gmail.com}
#' 
#' @seealso [alphahull::ahull()]
#' 
#' <https://github.com/azizka/speciesgeocodeR>
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
#' exclude.area=TRUE, country_map=land)
#' }
#' 
#' @import sf
#' 
#' @importFrom rnaturalearth ne_countries
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export EOO.computing
EOO.computing <- function(XY,
                          exclude.area = FALSE,
                          country_map = NULL,
                          export_shp = FALSE,
                          driver_shp = "ESRI Shapefile",
                          write_shp = FALSE,
                          alpha = 1,
                          buff.alpha = 0.1,
                          method.range = "convex.hull",
                          # Name_Sp = "species1",
                          method.less.than3 = "not comp",
                          write_results = FALSE,
                          file.name = "EOO.results",
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          proj_type = "cea",
                          mode = "spheroid"
) {
  
  list_data <- coord.check(XY = XY, check_eoo = TRUE)
  issue_long_span <- list_data$issue_long_span
  issue_nrow <- list_data$issue_nrow
  list_data <- list_data$list_data
  
  if (exclude.area) {
    ### Getting by default land map if country_map is not provided
    if (is.null(country_map)) {
      
      country_map <-
        rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
      
      country_map <- sf::st_make_valid(country_map)
      country_map <-
        suppressWarnings(sf::st_buffer(country_map, dist = 0))
      
    } else {
      
      if (any(!st_is_valid(country_map)))
        country_map <- sf::st_make_valid(country_map)
      
    }
    
  }
  
  res_df <-
    data.frame(eoo =  rep(NA, length(list_data)), 
               issue_eoo = rep(NA, length(list_data)))
  row.names(res_df) <- names(list_data)
  
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
  
  # if (is.null(names(list_data))) {
  #   names_ <-
  #     rep(Name_Sp, length(list_data))
  # } else {
  #   names_ <- names(list_data)
  # }
  
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
          # Name_Sp = names_[x],
          method.range = method.range,
          alpha = alpha,
          buff.alpha = buff.alpha,
          method.less.than3 = method.less.than3, 
          mode = mode, 
          proj_type = proj_type
        )
      
      names(res) <- c("eoo", "spatial")
      names(res)[1] <- list_data[[x]]$tax[1]
      
      # names(res)[1] <-
      #   paste0(names(res)[1], "_" , x)
      # if (length(res) > 1)
      #   names(res)[2] <-
      #   paste0(names(res)[2], "_" , x)
      
      res
      
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  res <- unlist(output[names(output) != "spatial"])
  
  res_df[which(row.names(res_df) %in% names(res)), 1] <-
    res
  
  if (length(issue_long_span) > 0)
    res_df[issue_long_span, 2] <-
    "Occurrences spans more than 180 degrees longitude, EOO unlikely reliable"
  
  if (length(issue_nrow) > 0)
    res_df[issue_nrow, 2] <-
    paste(res_df[issue_nrow, 2], "EOO cannot be estimated because less than 3 unique pair of coordinates", sep = "|")
    
  res_df$issue_eoo <- gsub("NA|", "", res_df$issue_eoo)
  
  if(export_shp) {
    
    output_spatial <- output[names(output) == "spatial"]
    output_spatial <- output_spatial[!is.na(output_spatial)]
    
    output_spatial <- 
        do.call("rbind", output_spatial)
      
    if (!is.null(output_spatial))
      row.names(output_spatial) <- 1:nrow(output_spatial)
  
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
    write.csv(res_df, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (!export_shp)
    output <- res_df

  if (export_shp)
    output <- list(results = res_df,
                   spatial = output_spatial)
  
  output
}
