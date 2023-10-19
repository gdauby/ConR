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
#' If `exclude.area` is TRUE and country_map is not provided, 
#' the world country polygons used comes from the package [rnaturalearth](https://www.rdocumentation.org/packages/rnaturalearth/versions/0.1.0/topics/ne_countries)
#' 
#' By default (`mode = "spheroid"`),the area of the polygon is based 
#' on a polygon in longitude/latitude coordinates considering the earth as an ellipsoid.  
#' 
#' To make a polygon more accurate, the function use the function `st_segmentize` 
#' from the [sf](https://CRAN.R-project.org/package=sf) package.
#' This adds vertices on the great circles (in order to make shortest distances between points, see example below) 
#' which can make difference for species with large distribution.  
#' 
#' An estimation of EOO based on projected data is also possible (`mode = "planar"`).
#' This allow the user to use its own projection.
#' By default, the projection (`proj_type = "cea"`) is a [global cylindrical equal-area projection](https://epsg.io/54034).
#' It can makes sense to use a planar projection to estimate the EOO. See example.
#' 
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
#' localized on a straight line (in which case EOO would be null), 'noises' are
#' added to coordinates, using the function
#' [jitter](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/jitter).
#' There is a warning when this happens. This means that EOO value will not be
#' constant across multiple estimation (although the variation should be small)
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
#' If `export_shp` is TRUE, a `list` with: \enumerate{ \item `results` a data.frame 
#' \item `spatial` used for EOO computation}
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
#' set.seed(2)
#' dataset <- dummy_dist(nsp = 3, max_occ = 30)
#' EOO.computing(XY = dataset, export_shp = TRUE)
#' 
#' EOO.computing(XY = dataset, mode = "planar")
#' 
#' EOO.computing(XY = dataset, method.range = "alpha")
#' 
#' res <- EOO.computing(XY = dataset, export_shp = TRUE)
#' res$spatial
#' 
#' 
#' data("land")
#' 
#' country_map <-
#'   rnaturalearth::ne_countries(scale = 50, returnclass = "sf", type = "map_units")
#' 
#' ### example Large distribution
#' # Spheroid estimation
#' par(mfrow = c(2,2))
#' XY <- 
#'   data.frame(x = c(-20, 5, 0, -20, -20), 
#'              y = c(-65, -45, -40, -55, -40),
#'              taxa = "species")
#' p1 <- 
#'   EOO.computing(XY =  XY, export_shp = TRUE)
#' p1$results
#' 
#' plot(st_geometry(country_map), extent = p1$spatial)
#' plot(p1$spatial, lwd = 2, col = "red",  add= TRUE)
#' 
#' p2 <- 
#'   EOO.computing(XY = XY, mode = "planar", export_shp = TRUE)
#' p2$results
#' plot(st_geometry(country_map), extent = p2$spatial)
#' plot(p2$spatial, lwd = 2, col = "red",  add= TRUE)
#' 
#' # World Mercartor projection
#' p2 <- 
#'   EOO.computing(XY = XY, mode = "planar", proj_type = 3395, export_shp = TRUE)
#' p2$results
#' plot(st_geometry(country_map), extent = p2$spatial)
#' plot(p2$spatial, lwd = 2, col = "red",  add= TRUE)
#' 
#' 
#' ### example Antartic distribution
#' XY <- 
#'   data.frame(x = c(-65, -62, -78, -65), 
#'              y = c(-150, 0, 120, 150),
#'              taxa = "species")
#' 
#' ## This throws an error
#' # p1 <- 
#' #   EOO.computing(XY = XY, mode = "planar", export_shp = TRUE) 
#' 
#' p1 <- 
#'   EOO.computing(XY = XY, mode = "planar", proj_type = "Antarctic", export_shp = TRUE)
#' 
#' 
#' p1 <- st_transform(p1$spatial, proj_crs(proj_type = "Antarctic"))
#' plot(st_geometry(st_transform(country_map, proj_crs(proj_type = "Antarctic"))), 
#'      extent = p1)
#' plot(p1, lwd = 2, col = 'red', add = TRUE)
#' 
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
                          # write_results = FALSE,
                          file.name = "EOO.results",
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          proj_type = "cea",
                          mode = "spheroid"
) {
  
  
  mode <- match.arg(mode, c("spheroid", "planar"))
  
  if (identical(mode, "planar"))
    proj_type <- proj_crs(proj_type = proj_type)
  
  list_data <- coord.check(XY = XY, 
                           check_eoo = TRUE, 
                           proj_type = if (!is.character(proj_type)) proj_type)
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
    data.frame(tax = names(list_data),
      eoo =  rep(NA, length(list_data)), 
               issue_eoo = rep(NA, length(list_data)))
  # row.names(res_df) <- names(list_data)
  
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
          # exclude.area = exclude.area,
          # country_map = country_map,
          # Name_Sp = names_[x],
          method.range = method.range,
          alpha = alpha,
          buff.alpha = buff.alpha,
          method.less.than3 = method.less.than3, 
          mode = mode, 
          proj_type = proj_type, 
          reproject = FALSE
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
  
  res_df[which(res_df$tax %in% names(res)), 2] <-
    res
  
  if (length(issue_long_span) > 0)
    res_df[issue_long_span, 3] <-
    "Occurrences spans more than 180 degrees longitude, EOO unlikely reliable"
  
  if (length(issue_nrow) > 0) {
    res_df[issue_nrow, 3] <-
      paste(res_df[issue_nrow, 3], "EOO cannot be estimated because less than 3 unique pair of coordinates", sep = "|")
    res_df[, 3] <- gsub("NA|", "", res_df[, 3])
  }
    
  if (export_shp | exclude.area) {
    
    output_spatial <- output[names(output) == "spatial"]
    output_spatial <- output_spatial[!is.na(output_spatial)]
    
    output_spatial <- 
        do.call("rbind", output_spatial)
      
    if (!is.null(output_spatial)) {
      
      row.names(output_spatial) <- 1:nrow(output_spatial)
      
      if (exclude.area) {
        
        if (identical(mode, "planar"))
          country_map <- sf::st_transform(country_map, crs = proj_crs(proj_type = proj_type))
        
        p1 <-
          suppressWarnings(suppressMessages(sf::st_intersection(output_spatial, 
                                                                st_make_valid(st_union(country_map)))))
        
        eoos <- data.frame(eoo = as.numeric(st_area(p1)) / 1000000,
                   tax = p1$tax)
        
        digits <-
          c(6, 5, 4, 3, 2, 1, 0)[findInterval(eoos$eoo, 
                                              c(0, 0.0001, 0.01, 0.1, 1, 10, 30000, Inf))]
        
        eoos$eoo <- round(eoos$eoo, digits)
        
        res_df[which(res_df$tax %in% eoos$tax), "eoo"] <- eoos$eoo
        output_spatial <- p1
        
      }
      
      if (identical(mode, "planar"))
        output_spatial <- sf::st_transform(output_spatial, crs = 4326)
      
    }
  }
  
  if(write_shp) {
    
    message("Writing EOO shapefiles in shapesIUCN directory")
    
    dir.create(file.path(paste(getwd(), "/shapesIUCN", sep = "")), showWarnings = FALSE)
    
    
    sf::write_sf(output_spatial,
                 dsn = "shapesIUCN",
                 layer = paste("EOO_poly", sep = ""),
                 driver = driver_shp,
                 overwrite = TRUE)
    
  }
  
  # if (write_results)
  #   utils::write.csv(res_df, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (!export_shp)
    output <- res_df

  if (export_shp)
    output <- list(results = res_df,
                   spatial = output_spatial)
  
  output
}
