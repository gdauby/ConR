
#' Internal function
#'
#' Build convex hull polygon
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @importFrom grDevices chull
#' @importFrom rgeos readWKT
#' @importFrom geosphere makePoly
#' 
.Convex.Hull.Poly <- function(XY) {
  hpts <- grDevices::chull(x =  XY[, 1], y = XY[, 2])
  hpts <- c(hpts, hpts[1])
  coord <- matrix(NA, length(hpts), 2)
  POLY <- "POLYGON(("
  for (i in 1:length(hpts)) {
    POLY <- paste(POLY, XY[hpts[i], 1], " ", XY[hpts[i], 2], sep = "")
    if (i != length(hpts))
      POLY <- paste(POLY, ", ", sep = "")
    if (i == length(hpts))
      POLY <- paste(POLY, "))", sep = "")
    
    coord[i, 1] <- XY[hpts[i], 2]
    coord[i, 2] <- XY[hpts[i], 1]
    
  }
  # p1 <- rgeos::readWKT(POLY)
  # raster::crs(p1) <- "+proj=longlat +datum=WGS84"
  # geosphere::makePoly(p1)
  
  p1 <- rgeos::readWKT(POLY)
  crs <- sp::CRS("+proj=longlat +datum=WGS84")
  raster::crs(p1) <- crs
  suppressWarnings(geosphere::makePoly(p1))
  
  
  return(p1)
}


#' Internal function
#'
#' Alpha hull processing
#'
#' @details 
#' The functions ahull_to_SPLDF and alpha.hull.poly were originally posted in the website https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
#' in a now broken link. It is also used in functions written by David Bucklin, see https://github.com/dnbucklin/r_movement_homerange 
#'
#'
#' @importFrom alphahull anglesArc
#' 
.ahull_to_SPLDF <- function(x, proj4string = NA)
{
  if (class(x) != 'ahull')
    stop('This function only works with `ahull` class objects')
  
  # convert ashape edges to DF
  x.ah.df <- as.data.frame(x$arcs)
  
  # convert each arc to a line segment
  l.list <- list()
  for (i in 1:nrow(x.ah.df))
  {
    # extract row i
    row_i <- x.ah.df[i, ]
    
    # extract elements for arc()
    v <- c(row_i$v.x, row_i$v.y)
    theta <- row_i$theta
    r <- row_i$r
    cc <- c(row_i$c1, row_i$c2)
    # from arc()
    angles <- alphahull::anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- cc[1] + r * cos(seqang)
    y <- cc[2] + r * sin(seqang)
    
    # convert to line segment
    l.list[[i]] <- sp::Line(cbind(x, y))
  }
  
  # promote to Lines class, then to SpatialLines class
  l <- sp::Lines(l.list, ID = 1)
  
  # copy over CRS data from original point data
  l.spl <-
    sp::SpatialLines(list(l), proj4string = sp::CRS(as.character(NA)))
  
  # promote to SpatialLinesDataFrame, required for export to GRASS / OGR
  l.spldf <-
    sp::SpatialLinesDataFrame(l.spl, data = data.frame(id = 1), match.ID =
                                FALSE)
  
  return(l.spldf)
}


#' Internal function
#'
#' Alpha hull process
#'
#' @details 
#' The functions ahull_to_SPLDF and alpha.hull.poly were originally posted in the website https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
#' in a now broken link. It is also used in functions written by David Bucklin, see https://github.com/dnbucklin/r_movement_homerange 
#'
#' @import sp raster
#' @importFrom rgeos gBuffer
#' @importFrom geosphere makePoly
#' @importFrom alphahull ahull
#' 
#' 
.alpha.hull.poly <- function(XY, alpha = 1, buff = 0.1) {
  
  Used_data = unique(XY)
  Used_data <- apply(Used_data, 2, jitter)
  if (any(rownames(utils::installed.packages()) == "alphahull")) {
    loadNamespace("alphahull")
    ahull.obj <-
      alphahull::ahull(Used_data[, c(1, 2)], alpha = alpha)
    y.as.spldf <- .ahull_to_SPLDF(ahull.obj)
    y.as.spldf_buff <- rgeos::gBuffer(y.as.spldf, width = buff)
    
    NZp <- methods::slot(y.as.spldf_buff, "polygons")
    holes <-
      lapply(NZp, function(x)
        sapply(methods::slot(x, "Polygons"), slot,
               "hole"))
    res <- lapply(1:length(NZp), function(i)
      methods::slot(NZp[[i]],
           "Polygons")[!holes[[i]]])
    IDs <- row.names(y.as.spldf_buff)
    NZfill <- sp::SpatialPolygons(lapply(1:length(res), function(i)
      sp::Polygons(res[[i]], ID = IDs[i])),
      proj4string = sp::CRS(sp::proj4string(y.as.spldf_buff)))
    
    
    crs <- sp::CRS("+proj=longlat +datum=WGS84")
    raster::crs(NZfill) <- crs
    
    return(NZfill)
    
    # raster::crs(NZfill) <- "+proj=longlat +datum=WGS84"
    
    
    
  } else{
    stop("The package alphahull is required for this procedure, please install it")
  }
}

#' Internal function
#'
#' Crop polygons
#'
#' @import sp raster
#' @importFrom geosphere areaPolygon
#' @importFrom sf st_intersection st_union
#' @importFrom methods as slot
#' 
.crop.poly <- function(poly, crop) {
  
  # @importFrom spatstat as.owin area.owin union.owin setminus.owin 
  
  poly_sf <- methods::as(poly, "sf")
  # crop_sf <- sf::st_combine(as(crop, "sf"))
  crop_sf <- methods::as(crop, "sf")
  
  suppressMessages(suppressWarnings(diff_croped <- 
                                      sf::st_intersection(crop_sf, poly_sf)))
  
  poly_masked <- methods::as(sf::st_union(diff_croped), "Spatial")
  
  # crs_crop <- raster::crs(crop)
  # 
  # raster::crs(poly) <- NA
  # raster::crs(crop) <- NA
  # 
  # p1_owin <- spatstat::as.owin(poly)
  # africa_owin <- spatstat::as.owin(crop)
  # 
  # if (round(spatstat::area.owin(spatstat::union.owin(p1_owin, africa_owin)), 3) != round(spatstat::area.owin(africa_owin), 3)) {
  #   w <- spatstat::setminus.owin(p1_owin, africa_owin)
  #   w2 <- spatstat::setminus.owin(p1_owin, w)
  #   poly_masked <- as(w2, "SpatialPolygons")
  #   
  #   raster::crs(poly_masked) <- crs_crop
  #   
  # } else{
  #   poly_masked <- poly
  #   
  # }
  
  EOO <-
    round(suppressWarnings(geosphere::areaPolygon(poly_masked)) / 1000000, 1)
  
  return(list(EOO, poly_masked))
}

#' Internal function
#'
#' EOO estimatiion
#'
#' @import sp raster
#' @importFrom  rgdal project
#' @importFrom  stats dist
#' @importFrom rgeos readWKT gBuffer
#' @importFrom geosphere makeLine areaPolygon
#' 
.EOO.comp <-  function(XY,
                       exclude.area = FALSE,
                       buff_width = 0.1,
                       country_map = NULL,
                       Name_Sp = "tax",
                       alpha.hull = FALSE,
                       convex.hull = TRUE,
                       alpha = 1,
                       buff.alpha = 0.1,
                       method.less.than3 = "not comp") {
  # , verbose=TRUE
  
  ### Checking if the method of calculating EOO has been chosen
  if (!convex.hull & !alpha.hull)
    stop("alpha.hull and convex.hull are both FALSE, choose one of them")
  
  if (nrow(unique(XY)) > 1)
    if (max(stats::dist(XY[, 2]), na.rm = T) >= 180)
      stop(
        paste(
          "EOO for species",
          as.character(Name_Sp),
          "cannot be computed because occurrences spans more than 180 degrees longitude"
        )
      )
  
  ## Check if there are less than 3 unique occurrences
  if (nrow(unique(XY)) < 3) {
    ## if there is only one occurrence, EOO is NA
    if (nrow(unique(XY)) < 2) {
      EOO <- NA
      message(
        paste(
          "EOO for",
          as.character(Name_Sp),
          "is not computed because there is only 1 unique occurrence"
        )
      )
      
    } else{
      if (method.less.than3 == "arbitrary") {
        
        projEAC <- .proj_crs()
        
        coordEAC <-
          as.data.frame(matrix(unlist(
            rgdal::project(
              as.matrix(unique(XY)[, 1:2]),
              proj = as.character(projEAC),
              inv = FALSE
            )
          ), ncol = 2))
        EOO <-
          as.numeric(stats::dist(coordEAC / 1000) * 0.1 * stats::dist(coordEAC / 1000))  #
      }
      
      if (method.less.than3 == "not comp") {
        ## if there are two unique occurences, EOO is not computed neither
        message(
          paste(
            "EOO for",
            as.character(Name_Sp),
            "is not computed because there is less than 3 unique occurrences"
          )
        )
        EOO <- NA
      }
    }
    
    OUTPUT <- round(EOO, 0)
    names(OUTPUT) <- c("EOO")
    
  } else{
    ### Checking if all occurrences are on a straight line
    if (length(unique(XY[, 1])) == 1 ||
        length(unique(XY[, 2])) == 1 ||
        round(abs(stats::cor(XY[, 1], XY[, 2])), 6) == 1) {
      ## If so, a straight line is built and a buffer of buff_width is added
      message(
        paste(
          "Occurrences of",
          as.character(Name_Sp),
          "follow a straight line, thus EOO is based on an artificial polygon using buff_width"
        )
      )
      hpts <- unique(XY[, c(2, 1)])
      POLY <- "LINESTRING("
      for (Z in 1:dim(hpts)[1]) {
        POLY <- paste(POLY, hpts[Z, 1], " ", hpts[Z, 2], sep = "")
        if (Z != dim(hpts)[1])
          POLY <- paste(POLY, ", ", sep = "")
        if (Z == dim(hpts)[1])
          POLY <- paste(POLY, ")", sep = "")
      }
      p1 <- rgeos::readWKT(POLY)
      
      p1 <- rgeos::readWKT(POLY)
      raster::crs(p1) <- sp::CRS("+proj=longlat +datum=WGS84")
      
      # crs <- CRS("+proj=longlat +datum=WGS84")
      # crs(p1) <- crs
      
      p1 <-
        suppressWarnings(geosphere::makeLine(p1)) ### Add vertices to line
      
      p1 <-
        suppressWarnings(rgeos::gBuffer(p1, width = buff_width)) ### Add buffer to line
      
      ## If exclude.area is TRUE
      if (exclude.area) {
        croped.EOO <- .crop.poly(poly = p1, crop = country_map)
        p1 <- croped.EOO[[2]]
        EOO <- croped.EOO[[1]]
      } else{
        EOO <- round(suppressWarnings(geosphere::areaPolygon(p1)) / 1000000,
                     1)
      }
      
    } else{
      
      if (alpha.hull)
        p1 <-
          .alpha.hull.poly(cbind(XY[, 2], XY[, 1]), alpha = alpha, buff = buff.alpha)
      
      if (convex.hull)
        p1 <- .Convex.Hull.Poly(cbind(XY[, 2], XY[, 1]))
      
      if (exclude.area) {
        croped.EOO <- 
          .crop.poly(poly = p1, crop = country_map)
        p1 <- croped.EOO[[2]]
      }
      
      ## If exclude.area is TRUE
      if (exclude.area) {
        EOO <-
          croped.EOO[[1]]
      } else{
        EOO <-
          round(suppressWarnings(geosphere::areaPolygon(p1)) / 1000000,
                0)
      }
    }
    
    OUTPUT <- list(EOO, p1)
    names(OUTPUT) <- c("EOO", "spatial.polygon")
  }
  
  # if(verbose) cat(" ",paste(Name_Sp,"EOO comp."))
  
  return(OUTPUT)
}


#' Extent of Occurrences
#' 
#' Compute extent of occurrences (EOO) for multiple taxa in square kilometres
#' using \code{\link[geosphere]{geosphere}} package and provide
#' \code{SpatialPolygons} used for EOO computation
#' 
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not
#' matter}
#' 
#' \tabular{ccc}{ [,1] \tab ddlat \tab numeric, latitude (in decimal
#' degrees)\cr [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#' [,3] \tab tax \tab character or factor, optinal field with taxa names\cr }
#' 
#' \strong{Important notes:}
#' 
#' EOO will only be computed if there is at least three unique occurrences
#' unless \code{method.less.than3} is put to "arbitrary". In that specific
#' case, EOO for species with two unique occurrences will be equal to
#' Dist*Dist*0.1 where Dist is the distance in kilometres separating the two
#' points.
#' 
#' For the very specific (and infrequent) case where all occurrences are
#' localized on a straight line (in which case EOO would be null), EOO is
#' estimated by the area of polygon surrounding this straight line with a
#' buffer of \code{buff.alpha} decimal degree. There is a warning when this
#' happen.
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
#' @param buff_width a numeric. For a specific case where all points of a taxa
#' are on a straight line, see Details. By default is 0.1
#' @param method.less.than3 a character string. If equal to "arbitrary", will
#' give a value to species with two unique occurrences, see Details. By default
#' is "not comp"
#' @param write_results a logical. If TRUE, results will be exported in the
#' working environment as a csv file. By default is TRUE
#' @param file.name a character string. Name file for exported results in csv
#' file. By default is "EOO.results"
#' @param parallel a logical. Wether running in parallel. By default, it is
#' FALSE
#' @param NbeCores an integer. Register the number of cores for parallel
#' execution. By default, it is 2
#' @param show_progress logical. Whether a progress bar should displayed. TRUE by default.
#' 
#' @return If \code{export_shp} is FALSE, a \code{dataframe} with one field
#' containing EOO in square kilometres.  \code{NA} is given when EOO could not
#' be computed because there is less than three unique occurrences (or two if
#' \code{method.less.than3} is put to "arbitrary").
#' 
#' If \code{export_shp} is TRUE, a \code{list} with: \enumerate{ \item EOO in
#' square kilometres \item \code{SpatialPolygons} used for EOO computation }
#' @author Gilles Dauby
#' 
#' \email{gildauby@@gmail.com}
#' @seealso \code{\link[alphahull]{ahull}}
#' 
#' \url{https://github.com/azizka/speciesgeocodeR}
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges,
#' Journal of Applied Ecology, 49 1-9
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
#' @import sp raster
#' @importFrom rgdal writeOGR
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
                          write_shp = FALSE,
                          alpha = 1,
                          buff.alpha = 0.1,
                          method.range = "convex.hull",
                          Name_Sp = "species1",
                          buff_width = 0.1,
                          method.less.than3 = "not comp",
                          write_results = TRUE,
                          file.name = "EOO.results",
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE
){ # , verbose=TRUE
  
  if (any(is.na(XY[, c(1:2)]))) {
    print(paste(
      "Skipping",
      length(which(rowMeans(is.na(
        XY[, 1:2]
      )) > 0)) ,
      "occurrences because of missing coordinates for",
      # if(verbose)
      paste(as.character(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                           0), 3])), collapse = " AND ")
    ))
    XY <- XY[which(!is.na(XY[, 1])),]
    XY <- XY[which(!is.na(XY[, 2])),]
  }
  
  XY <- as.data.frame(XY)
  
  if (exclude.area &
      is.null(country_map))
    stop("exclude.area is TRUE but no country_map is provided")
  if (buff_width > 80)
    stop("buff_width has unrealistic value")
  if (any(XY[, 2] > 180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 1] > 180))
    stop("coordinates are outside of expected range")
  
  if (method.range == "convex.hull") {
    convex.hull = TRUE
    alpha.hull = FALSE
  }
  
  if (method.range == "alpha.hull") {
    convex.hull = FALSE
    alpha.hull = TRUE
  }
  
  if (ncol(XY) > 2) {
    colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
    XY$tax <- as.character(XY$tax)
    list_data <- split(XY, f = XY$tax)
  } else{
    colnames(XY)[1:2] <- c("ddlat", "ddlon")
    list_data <- list(XY)
  }
  
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('doParallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else{
    `%d%` <- foreach::`%do%`
  }
  
  
  if (is.null(names(list_data))) {
    names_ <-
      rep(Name_Sp, length(list_data))    
  } else {
    names_ <- names(list_data)    
  }

  x <- NULL
  if(show_progress) {
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(list_data),
                            style = 3)
    
    progress <- function(n)
      utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{opts <- NULL}
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      # source("./R/IUCNeval.functionv11.R")
      
      if (!parallel & show_progress)
        utils::setTxtProgressBar(pb, x)
      
      res <-
        .EOO.comp(
          XY = list_data[[x]],
          exclude.area = exclude.area,
          buff_width = buff_width,
          country_map = country_map,
          Name_Sp = names_[x],
          alpha.hull = alpha.hull,
          convex.hull = convex.hull,
          alpha = alpha,
          buff.alpha = buff.alpha,
          method.less.than3 = method.less.than3
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
  
  if(write_shp) {
    
    dir.create(file.path(paste(getwd(), "/shapesIUCN", sep = "")), showWarnings = FALSE)
    output_spatial <- unlist(output[grep("spatial", names(output))])
    id_spatial <-
      as.numeric(unlist(lapply(strsplit(
        names(output_spatial), "_"
      ), function(x)
        x[[2]])))
    
    for (i in 1:length(output_spatial)) {
      
      if (length(list.files(paste(getwd(), "/shapesIUCN", sep = ""))) >
          0) {
        if (length(grep(paste(names(output)[i], "_EOO_poly", sep = ""), unique(sub(
          "....$", '', list.files(paste(getwd(), "/shapesIUCN", sep = ""))
        )))) > 0)
        {
          FILES <-
            list.files(paste(getwd(), "/shapesIUCN", sep = ""), full.names = TRUE)
          file.remove(FILES[grep(paste(names(output)[i], "_EOO_poly", sep =
                                         ""), FILES)])
        }
      }
      NAME <- names_[id_spatial[i]]
      output_spatial[[i]]@polygons[[1]]@ID <- "1"
      ConvexHulls_poly_dataframe <-
        sp::SpatialPolygonsDataFrame(output_spatial[[i]], data = as.data.frame(names(output_spatial[[i]])))
      colnames(ConvexHulls_poly_dataframe@data) <-
        paste(substr(names_[id_spatial[i]], 0, 3), collapse = '')
      rgdal::writeOGR(
        ConvexHulls_poly_dataframe,
        "shapesIUCN",
        paste(names_[id_spatial[i]], "_EOO_poly", sep = ""),
        driver = "ESRI Shapefile",
        overwrite_layer = TRUE
      )
    }
  }
  
  if (write_results)
    utils::write.csv(Results_short, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (!export_shp)
    output <- Results_short
  
  output
}

#' Internal function
#'
#' subpopulations estimatiion
#'
#' @importFrom  rgdal project
#' @importFrom rgeos readWKT gBuffer gUnion
#' 
.subpop.comp <- function(XY,
                         Resol_sub_pop) {
 
  projEAC <- .proj_crs()
  
  XY <- XY[, c(2, 1)]
  
  coordEAC <-
    as.data.frame(matrix(unlist(
      rgdal::project(as.matrix(XY), proj = as.character(projEAC), inv = FALSE)
    ), ncol = 2))
  rownames(coordEAC) <- seq(1, nrow(coordEAC), 1)
  
  p2 <-
    rgeos::readWKT(paste("POINT(", mean(unique(coordEAC)[1, 1]), " ", mean(unique(coordEAC)[1, 2]), ")", sep =
                           ""))
  
  
  p2_Buffered1 <-
    rgeos::gBuffer(p2, width = Resol_sub_pop * 1000, id = 1)
  if (nrow(unique(coordEAC)) > 1) {
    for (LL in 2:nrow(unique(coordEAC))) {
      p2 <-
        rgeos::readWKT(paste("POINT(", mean(unique(coordEAC)[LL, 1]), " ", mean(unique(coordEAC)[LL, 2]), ")", sep =
                               ""))
      p2_Buffered <-
        rgeos::gBuffer(p2, width = Resol_sub_pop * 1000, id = LL)
      p2_Buffered1 <- rgeos::gUnion(p2_Buffered1, p2_Buffered)
    }
  }
  
  splited_pol <-
    lapply(p2_Buffered1@polygons, slot, "Polygons")[[1]]
  
  NbeSubPop <- length(splited_pol)
  
  SubPopPoly <-
    sp::SpatialPolygons(
      Srl = list(p2_Buffered1@polygons[[1]]),
      pO = as.integer(1),
      proj4string = projEAC
    )
  
  SubPopPoly <-
    sp::spTransform(SubPopPoly,
                    sp::CRS("+proj=longlat +datum=WGS84"))
  
  OUTPUT <- list(NbeSubPop, SubPopPoly)
  names(OUTPUT) <- c("Number of subpopulation", "subpop.poly")
  return(OUTPUT)
}

#' Number of subpopulations
#'
#' Estimate the number of locations following the method **circular buffer method**
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY string, indicating the method used for estimating the number of locations. Either "fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"
#' @param Resol_sub_pop numeric. Defines in kilometres the radius of the circles around each occurrence
#' 
#' @details 
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not matter}
#' 
#' \tabular{ccc}{
#'   [,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#'   [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#'   [,3] \tab tax \tab character or factor, optinal field with taxa names\cr
#' }
#' 
#' @references Rivers MC, Bachman SP, Meagher TR, Lughadha EN, Brummitt NA (2010) Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation 19: 2071-2085. doi: 10.1007/s10531-010-9826-9
#'
#' @return A list with one list for each taxa containing [[1]]Number of subpopulation and [[2]]SpatialPolygons.
#' 
#' @examples 
#' data(dataset.ex)
#' \dontrun{
#'subpop <- subpop.comp(dataset.ex, Resol_sub_pop = 30)
#'}
#'
#' 
#' @export
subpop.comp <- function(XY, Resol_sub_pop = NULL) {
  if (is.null(Resol_sub_pop))
    stop("Resol_sub_pop is missing, please provide a value")
  
  if (any(is.na(XY[, c(1, 2)]))) {
    length(which(rowMeans(is.na(XY[, 1:2])) > 0))
    unique(XY[which(rowMeans(is.na(XY[, 1:2])) > 0), 3])
    print(
      paste(
        "Skipping",
        length(which(rowMeans(is.na(
          XY[, 1:2]
        )) > 0)) ,
        "occurrences because of missing coordinates for",
        paste(as.character(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                             0), 3])), collapse = " AND ")
      )
    )
    XY <- XY[which(!is.na(XY[, 1])), ]
    XY <- XY[which(!is.na(XY[, 2])), ]
  }
  
  if (any(XY[, 1] > 180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 2] > 180))
    stop("coordinates are outside of expected range")
  
  colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
  XY$tax <- as.character(XY$tax)
  list_data <- split(XY, f = XY$tax)
  
  OUTPUT <-
    lapply(list_data, function(x)
      .subpop.comp(XY = x, Resol_sub_pop = Resol_sub_pop))
  if (length(OUTPUT) == 1)
    OUTPUT <- OUTPUT[[1]]
  return(OUTPUT)
}

#' Internal function
#'
#' AOO estimatiion
#'
#' 
.AOO.estimation <- function(coordEAC, 
                            cell_size = 2, 
                            nbe_rep = 0, 
                            # poly_borders = NULL, 
                            export_shp=FALSE) {
  
  
  crs_proj <- 
    .proj_crs()
  
  res <-
    .cell.occupied(
      nbe_rep = nbe_rep,
      size = cell_size,
      coord = coordEAC,
      export_shp = export_shp
    )
  
  Corners <- rbind(c(min(coordEAC[, 1]),
                     max(coordEAC[, 1])),
                   c(min(coordEAC[, 2]),
                     max(coordEAC[, 2])))
  
  # if (nbe_rep == 0) {
  #   
  #   Occupied_cells <- vector(mode = "numeric", length = 4)
  #   decal <- c(0, 1, 2, 3)
  #   
  #   for (h in decal) {
  #     ext <-
  #       raster::extent(
  #         floor(Corners[1, 1]) - h * (cell_size * 1000 / 4) - 2 * cell_size * 1000,
  #         floor(Corners[1, 2]) + h * (cell_size * 1000 /
  #                                       4) + 2 * cell_size * 1000,
  #         floor(Corners[2, 1]) - h * (cell_size * 1000 /
  #                                       4) - 2 * cell_size * 1000,
  #         floor(Corners[2, 2]) + h * (cell_size * 1000 /
  #                                       4) + 2 * cell_size * 1000
  #       )
  #     
  #     r <-
  #       raster::raster(ext, resolution = cell_size * 1000, crs = crs_proj)
  #     
  #     r2_AOO <-
  #       raster::rasterize(coordEAC[, 1:2], r)
  #     
  #     OCC <-
  #       length(which(!is.na(raster::values(r2_AOO))))
  #     
  #     Occupied_cells[h + 1] <- OCC
  #     
  #     ### If only one occupied cell, stop the production of raster
  #     if (OCC == 1)
  #       break
  #   }
  #   # h <- decal[which.min(Occupied_cells)]
  #   # Occupied_cells <- min(Occupied_cells)
  # }
  # 
  # if (nbe_rep > 0) {
  #   Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
  #   
  #   for (h in 1:nbe_rep) {
  #     rd.1 <- runif(1) * cell_size * 1000
  #     rd.2 <- runif(1) * cell_size * 1000
  #     
  #     ext = raster::extent(
  #       floor(Corners[1, 1]) - rd.1 - 2 * cell_size * 1000,
  #       floor(Corners[1, 2]) + rd.1 + 2 * cell_size * 1000,
  #       floor(Corners[2, 1]) - rd.2 - 2 * cell_size *
  #         1000,
  #       floor(Corners[2, 2]) + rd.2 + 2 * cell_size * 1000
  #     )
  #     r = raster::raster(ext, resolution = cell_size * 1000, crs = crs_proj)
  #     # r
  #     r2_AOO <- raster::rasterize(coordEAC[, 1:2], r)
  #     OCC <- length(which(!is.na(raster::values(r2_AOO))))
  #     Occupied_cells[h] <- OCC
  #     # rd.1.vec <- c(rd.1.vec, rd.1)
  #     # rd.2.vec <- c(rd.2.vec, rd.2)
  #     if (OCC == 1)
  #       break
  #   }
  #   
  # }
  
  # Occupied_cells <- Occupied_cells[Occupied_cells>0]
  # Occupied_cells <- min(Occupied_cells)
  
  AOO <- res[[2]] * cell_size * cell_size  ### AOO
  if (export_shp)
    return(list(AOO, res[[1]]))
  if (!export_shp)
    return(AOO)
  
}


#' Internal function
#'
#' get proj CRS
#'
#' @importFrom utils packageVersion
#' @importFrom rgdal showWKT
#' 
.proj_crs <- function() {
  
  ## https://epsg.io/54032
  # World Azimuthal Equidistant
  proj <-
    "+proj=aeqd +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    
  ## https://epsg.io/102022
  ## Africa Albers Equal Area Conic
  # proj <- 
  #   "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  # "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  if (utils::packageVersion("sp") >= "1.3.3") {
    wkt_crs <-
      rgdal::showWKT(
        proj
      )
    crs_proj <- sp::CRS(projargs = proj, 
                        SRS_string = wkt_crs)
  }
  
  if (utils::packageVersion("sp") < "1.3.3")
    crs_proj <-
      sp::CRS(projargs = "+proj=eqc +lat_ts=60 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  
  return(crs_proj)
}


#' Area of occupancy
#'
#' Compute areas of occupancy (AOO) for multiple taxa in square kilometres
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param XY \code{"dataframe"} see Details
#' @param Cell_size_AOO numeric, value indicating the grid size in kilometres used for estimating Area of Occupancy.  By default, equal to 2
#' @param nbe.rep.rast.AOO numeric , indicate the number of raster with random starting position for estimating the AOO. By default, it is 0 but some minimal translation of the raster are still done
#' @param parallel logical, wether running in parallel. By default, it is FALSE
#' @param NbeCores string integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
#' @param export_shp logical, whether a shapefile of occupied cells should be exported. By default, it is FALSE
#'
#' @details 
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not matter}
#' 
#' \tabular{ccc}{
#'   [,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#'   [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#'   [,3] \tab tax \tab character or factor, optinal field with taxa names\cr
#' }
#' 
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9
#'
#' @return 
#' If \code{export_shp} if FALSE a vector of AOO estimates for each taxa
#' If \code{export_shp} if TRUE a list with two elements
#' \enumerate{
#'   \item a vector of AOO estimates for each taxa
#'   \item a list of SpatialPolygonsDataFrame for each taxa
#' }
#'   
#' @examples 
#' data(dataset.ex)
#' \dontrun{
#'AOO <- AOO.computing(dataset.ex)
#'}
#'
#'# This would estimate AOO for all taxa by overlaying randomly a 
#'# grid 100 times. For each taxa, the minimum value is kept
#' \dontrun{
#'AOO <- AOO.computing(dataset.ex, nbe.rep.rast.AO = 100)
#'}
#'
#' @importFrom rgdal project
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export
AOO.computing <- function(XY,
                          Cell_size_AOO = 2,
                          nbe.rep.rast.AOO = 0,
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          export_shp = FALSE
) {
  
  if (!any(class(XY) == "data.frame"))
    XY <- as.data.frame(XY)
  if (any(XY[, 2] > 180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 1] > 180))
    stop("coordinates are outside of expected range")
  
  projEAC <- .proj_crs()
  
  coordEAC <-
    data.frame(matrix(unlist(
      rgdal::project(as.matrix(XY[, c(2, 1)]),
                     proj = as.character(projEAC), inv =
                       FALSE)
    ),
    ncol = 2),
    tax = XY[, 3])
  
  ## if any missing coordinates
  if (any(is.na(coordEAC[, c(1:2)]))) {
    print(paste(
      "Skipping",
      length(which(rowMeans(is.na(
        coordEAC[, 1:2]
      )) > 0)),
      "occurrences because of missing coordinates for",
      paste(as.character(unique(coordEAC[which(rowMeans(is.na(coordEAC[, 1:2])) >
                                                 0), 3])), collapse = " AND ")
    ))
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 1])), ]
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 2])), ]
  }
  
  coordEAC$tax <- as.character(coordEAC$tax)
  list_data <- split(coordEAC, f = coordEAC$tax)
  
  # if(show_progress) prog. <- "text"
  # if(!show_progress) prog. <- "none"
  
  
  if(parallel) {
    # if("doParallel" %in% 
    #    rownames(installed.packages()) == FALSE) {stop("Please install doParallel package")}
    # 
    # library(doParallel)
    
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    # registerDoParallel(NbeCores)
    message('doParallel running with ',
            NbeCores, ' cores')
    `%d%` <- foreach::`%dopar%`
  }else{
    `%d%` <- foreach::`%do%`
  }
  
 
  
  x <- NULL
  if(show_progress) {
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(list_data),
                            style = 3)
    
    progress <- function(n)
      utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  }else{opts <- NULL}
  
  output <-
    foreach::foreach(
      x = 1:length(list_data),
      .combine = 'c', .options.snow = opts
    ) %d% {
      if (!parallel & show_progress)
        utils::setTxtProgressBar(pb, x)
      # source("./R/IUCNeval.functionv11.R")

      res <- .AOO.estimation(
        coordEAC = list_data[[x]],
        cell_size = Cell_size_AOO,
        nbe_rep = nbe.rep.rast.AOO,
        export_shp = export_shp
      )

      if (export_shp)
        names(res) <- c("aoo", "spatial")

      res
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  if(!export_shp) {
    
    res <- unlist(output)
    names(res) <- names(list_data)
    
  }
  
  if(export_shp) {
    
    res <- unlist(output[names(output) == "aoo"])
    names(res) <- names(list_data)
    shapes <-  unlist(output[names(output) == "spatial"])
    names(shapes) <- names(list_data)
    
  }
  
  if(!export_shp) return(res)
  if(export_shp) return(list(res, shapes))
  
}


#' Internal function
#'
#' Count number of occupied cells given resolution, projection
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' 
.cell.occupied <-
  function(nbe_rep = 0,
           size = 4,
           coord,
           export_shp = TRUE) {
    
    crs_proj <- 
      .proj_crs()
    
    Corners <- rbind(c(min(coord[, 1]),
                       max(coord[, 1])),
                     c(min(coord[, 2]),
                       max(coord[, 2])))
    
    if (nbe_rep == 0) {
      
      Occupied_cells <- vector(mode = "numeric", length = 4)
      decal <- c(0, 1, 2, 3)
      
      if (abs(Corners[1, 1] - Corners[1, 2]) > abs(Corners[2, 1] - Corners[2, 2])) {
        longer <-
          abs(Corners[1, 1] - Corners[1, 2])
      } else{
        longer <-
          abs(Corners[2, 1] - Corners[2, 2])
      }
      
      
      for (h in decal) {
        
        # xmin <- 
        #   floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000
        # xmax <- 
        #   xmin + longer + h * (size * 1000 / 4) + 2
        # 
        # ymin <- 
        #   floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *size * 1000
        # ymax <- 
        #   ymin + longer + h * (size * 1000 / 4) + 2
        # 
        # ext <-
        #   raster::extent(
        #     xmin,
        #     xmax,
        #     ymin,
        #     ymax
        #   )
        
        ext <-
          raster::extent(
            floor(Corners[1, 1]) - h * (size * 1000 / 4) - 2 * size * 1000,
            floor(Corners[1, 2]) + h * (size * 1000 / 4) + 2 *
              size * 1000,
            floor(Corners[2, 1]) - h * (size * 1000 / 4) - 2 *
              size * 1000,
            floor(Corners[2, 2]) + h * (size * 1000 / 4) + 2 *
              size * 1000
          )
        
        r <-
          raster::raster(ext, 
                         resolution = size * 1000, 
                         crs = crs_proj)
        
        r2_ <-
          raster::rasterize(coord[, 1:2], r)
        
        OCC <-
          length(which(!is.na(raster::values(r2_))))
        
        Occupied_cells[h + 1] <- OCC
        
        ### If only one occupied cell, stop the production of raster
        if (OCC == 1)
          break
      }
      # h <- decal[which.min(Occupied_cells)]
      # Occupied_cells <- min(Occupied_cells)
    }
    
    if (nbe_rep > 0) {
      Occupied_cells <- vector(mode = "numeric", length = nbe_rep)
      
      for (h in 1:nbe_rep) {
        rd.1 <- stats::runif(1) * size * 1000
        rd.2 <- stats::runif(1) * size * 1000
        
        ext = raster::extent(
          floor(Corners[1, 1]) - rd.1 - 2 * size * 1000,
          floor(Corners[1, 2]) + rd.1 + 2 * size * 1000,
          floor(Corners[2, 1]) - rd.2 - 2 * size * 1000,
          floor(Corners[2, 2]) + rd.2 + 2 * size * 1000
        )
        r = raster::raster(ext, resolution = size * 1000, crs = crs_proj)
        # r
        r2_ <- raster::rasterize(coord[, 1:2], r)
        OCC <- length(which(!is.na(raster::values(r2_))))
        Occupied_cells[h] <- OCC
        # rd.1.vec <- c(rd.1.vec, rd.1)
        # rd.2.vec <- c(rd.2.vec, rd.2)
        if (OCC == 1)
          break
      }
      
    }
    
    Occupied_cells <- Occupied_cells[Occupied_cells > 0]
    Occupied_cells <- min(Occupied_cells)
    
    if (export_shp)
      r2_ <-
      raster::projectRaster(from = r2_, 
                            crs = "+proj=longlat +datum=WGS84 +no_defs")
    
    if (export_shp)
      r2_pol <-
      raster::rasterToPolygons(
        r2_,
        fun = NULL,
        n = 4,
        na.rm = TRUE,
        digits = 6,
        dissolve = FALSE
      )
    
    if (export_shp)
      return(list(r2_pol, Occupied_cells))
    if (!export_shp)
      return(list(NA, Occupied_cells))
    
  }

#' Number of locations
#'
#' Estimate the number of locations for multiple taxa
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#'
#' @param method string, indicating the method used for estimating the number of locations. Either "fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"
#' @param nbe_rep numeric , indicate the number of raster with random starting position for estimating the number of locations By default, it is 0 but some minimal translation of the raster are still done
#' @param protec.areas \code{SpatialPolygonsDataFrame}, shapefile with protected areas. If provided, this will be taken into account for calculating number of location (see Details and \code{method_protected_area}). By default, no shapefile is provided
#' @param Cell_size_locations numeric, value indicating the grid size in kilometres used for estimating the number of location. By default, equal to 10
#' @param method_protected_area string, by default is "no_more_than_one"", which means occurrences within protected areas (if provided) will not be taken into account for estimating the number of locations following the grid system, see Details. By default, it is "no_more_than_one"
#' @param ID_shape_PA string, indicating the field name of \code{protec.areas} with ID of the \code{SpatialPolygonsDataFrame} of protected areas
#' @param Rel_cell_size numeric, if \code{method_locations="sliding scale"}, \code{Cell_size_locations} is ignored and the resolution is given by the maximum distance separating two occurrences multiplied by \code{Rel_cell_size}. By default, it is 0.05
#' @param parallel logical, wether running in parallel. By default, it is FALSE
#' @param NbeCores string integer, register the number of cores for parallel execution. By default, it is 2
#' @param show_progress logical, whether a bar showing progress in computation should be shown. By default, it is TRUE
#' 
#' @details 
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not matter}
#' 
#' \tabular{ccc}{
#'   [,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#'   [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#'   [,3] \tab tax \tab character or factor, optinal field with taxa names\cr
#' }
#' 
#' @references Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9
#'
#' @return A list with one list for each species containing [[1]]SpatialPolygonDataframe and [[2]]vector of the number of location.
#' 
#' @examples 
#' data(dataset.ex)
#' \dontrun{
#'locations <- locations.comp(dataset.ex)
#'}
#'
#'# This would estimate AOO for all taxa by overlaying 
#'# randomly a grid 100 times. For each taxa, the minimum value is kept
#' \dontrun{
#'AOO <- AOO.computing(dataset.ex, nbe.rep.rast.AO = 100)
#'}
#'
#' @importFrom rgdal project
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' 
#' @export
locations.comp <- function(XY,
                           method = "fixed_grid",
                           nbe_rep = 0,
                           protec.areas = NULL,
                           Cell_size_locations = 10,
                           method_protected_area = "no_more_than_one",
                           ID_shape_PA = "WDPA_PID",
                           Rel_cell_size = 0.05,
                           parallel = FALSE,
                           NbeCores = 2,
                           show_progress = TRUE) {
  
  if (!any(class(XY) == "data.frame"))
    XY <- as.data.frame(XY)
  if (any(XY[, 2] > 180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 1] > 180))
    stop("coordinates are outside of expected range")
  
  projEAC <- .proj_crs()
  
  coordEAC <-
    data.frame(matrix(unlist(
      rgdal::project(as.matrix(XY[, c(2, 1)]),
                     proj = as.character(projEAC), inv =
                       FALSE)
    ),
    ncol = 2),
    tax = XY[, 3])
  
  
  ## if any missing coordinates
  if (any(is.na(coordEAC[, c(1:2)]))) {
    print(
      paste(
        "Skipping",
        length(which(rowMeans(
          is.na(coordEAC[, 1:2])
        ) > 0)),
        "occurrences because of missing coordinates for",
        paste(as.character(unique(coordEAC[which(rowMeans(is.na(coordEAC[, 1:2])) >
                                                   0), 3])), collapse = " AND ")
      )
    )
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 1])),]
    coordEAC <- coordEAC[which(!is.na(coordEAC[, 2])),]
  }
  
  coordEAC$tax <- as.character(coordEAC$tax)
  list_data <- split(coordEAC, f = coordEAC$tax)
  
  
  crs_proj <- projEAC
  
  ## geographical distances for all pairs of occurrences
  
  if (is.null(protec.areas)) {
    if (nrow(coordEAC) > 1)
      pairwise_dist <- stats::dist(coordEAC[, 1:2],  upper = F)
    
    ## resolution definition
    if (any(method == "fixed_grid"))
      Resolution <- Cell_size_locations
    if (any(method == "sliding scale")) {
      if (nrow(coordEAC) > 1) {
        Resolution <- max(pairwise_dist) * Rel_cell_size / 1000
      } else{
        Resolution <- 10
      }
    }
    
    if (parallel) {
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      
      # registerDoParallel(NbeCores)
      message('doParallel running with ',
              NbeCores, ' cores')
      
      `%d%` <- foreach::`%dopar%`
    } else{
      `%d%` <- foreach::`%do%`
    }
    
    x <- NULL
    
    if(show_progress) {
      pb <-
        utils::txtProgressBar(min = 0,
                              max = length(list_data),
                              style = 3)
      
      progress <- function(n)
        utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    }else{opts <- NULL}
    
    output <-
      foreach::foreach(
        x = 1:length(list_data),
        .combine = 'c',
        .options.snow = opts
      ) %d% {
        
        if (!parallel & show_progress)
          utils::setTxtProgressBar(pb, x)
        
        res <- .cell.occupied(
          size = Resolution,
          coord = list_data[[x]],
          nbe_rep = nbe_rep
        )
        
        names(res) <- c("spatial", "nbe_occ")
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    Locations <- unlist(output[names(output) == "nbe_occ"])
    r2 <- unlist(output[names(output) == "spatial"])[[1]]
    names(Locations) <-
      names(r2) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
  }
  
  if (!is.null(protec.areas)) {
    ### Taking into account Protected Areas if provided
    DATA_SF <- as.data.frame(XY[, 1:2])
    colnames(DATA_SF) <- c("ddlat", "ddlon")
    sp::coordinates(DATA_SF) <-  ~ ddlon + ddlat
    raster::crs(DATA_SF) <- raster::crs(protec.areas)
    Links_NatParks <- sp::over(DATA_SF, protec.areas)
    
    coordEAC_pa <- coordEAC[!is.na(Links_NatParks[, 1]),]
    coordEAC_pa <-
      cbind(coordEAC_pa, id_pa = Links_NatParks[which(!is.na(Links_NatParks[, 1])), ID_shape_PA])
    
    LocNatParks <-
      vector(mode = "numeric", length = length(list_data))
    names(LocNatParks) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
    if (nrow(coordEAC_pa) > 0) {
      if (method_protected_area == "no_more_than_one") {
        ## if method is 'no_more_than_one' the number of location is the number of occupied protected areas
        
        loc_pa <-
          by(
            coordEAC_pa[, c("tax", "id_pa")],
            coordEAC_pa[, "tax"],
            FUN = function(x)
              length(unique(x$id_pa))
          )
        names(loc_pa) <-
          gsub(pattern = " ",
               replacement = "_",
               names(loc_pa))
        LocNatParks[names(LocNatParks) %in% names(loc_pa)] <-
          loc_pa
        r2_PA <- NA
        
      } else{
        coordEAC_pa$tax <- as.character(coordEAC_pa$tax)
        list_data_pa <- split(coordEAC_pa, f = coordEAC_pa$tax)
        
        ## geographical distances for all pairs of occurrences
        if (nrow(coordEAC_pa) > 1)
          pairwise_dist_pa <- stats::dist(coordEAC_pa[, 1:2],  upper = F)
        
        ## resolution definition
        if (any(method == "fixed_grid"))
          Resolution <- Cell_size_locations
        if (any(method == "sliding scale")) {
          if (nrow(coordEAC_pa) > 1) {
            Resolution <- max(pairwise_dist_pa) / 1000 * Rel_cell_size
          } else{
            Resolution <- 10
          }
        }
        
        if (parallel) {
          cl <- snow::makeSOCKcluster(NbeCores)
          doSNOW::registerDoSNOW(cl)
          
          # registerDoParallel(NbeCores)
          message('doParallel running with ',
                  NbeCores, ' cores')
          
          `%d%` <- foreach::`%dopar%`
        } else{
          `%d%` <- foreach::`%do%`
        }
        
        x <- NULL
        pb <- 
          utils::txtProgressBar(min = 0, max = length(list_data_pa), style = 3)
        progress <- function(n)
          utils::setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        output <-
          foreach::foreach(
            x = 1:length(list_data_pa),
            .combine = 'c',
            .options.snow = opts
          ) %d% {
            
            if (!parallel)
              utils::setTxtProgressBar(pb, x)
            
            res <- .cell.occupied(
              size = Resolution,
              coord = list_data_pa[[x]],
              nbe_rep = nbe_rep
            )
            
            names(res) <- c("spatial", "nbe_occ")
            res
          }
        
        if(parallel) snow::stopCluster(cl)
        if(show_progress) close(pb)
        
        loc_pa <- unlist(output[names(output) == "nbe_occ"])
        r2_PA <- unlist(output[names(output) == "spatial"])[[1]]
        names(loc_pa) <-
          names(r2_PA) <-
          gsub(pattern = " ",
               replacement = "_",
               names(list_data_pa))
        LocNatParks[names(LocNatParks) %in% names(loc_pa)] <-
          loc_pa
      }
    } else{
      r2_PA <- NA
    }
    
    coordEAC_not_pa <- coordEAC[is.na(Links_NatParks[, 1]),]
    LocOutNatParks <-
      vector(mode = "numeric", length = length(list_data))
    names(LocOutNatParks) <-
      gsub(pattern = " ",
           replacement = "_",
           names(list_data))
    
    if (nrow(coordEAC_not_pa) > 0) {
      coordEAC_not_pa$tax <- as.character(coordEAC_not_pa$tax)
      list_data_not_pa <-
        split(coordEAC_not_pa, f = coordEAC_not_pa$tax)
      
      ## geographical distances for all pairs of occurrences
      if (nrow(coordEAC_pa) > 1)
        pairwise_dist_not_pa <- stats::dist(coordEAC_not_pa[, 1:2],  upper = F)
      
      ## resolution definition
      if (any(method == "fixed_grid"))
        Resolution <- Cell_size_locations
      if (any(method == "sliding scale")) {
        if (nrow(coordEAC_pa) > 1) {
          Resolution <- max(pairwise_dist_not_pa) * Rel_cell_size
        } else{
          Resolution <- 10
        }
      }
      
      if (parallel) {
        cl <- snow::makeSOCKcluster(NbeCores)
        doSNOW::registerDoSNOW(cl)
        
        # registerDoParallel(NbeCores)
        message('doParallel running with ',
                NbeCores, ' cores')
        
        `%d%` <- foreach::`%dopar%`
      } else{
        `%d%` <- foreach::`%do%`
      }
      
      x <- NULL
      pb <- 
        utils::txtProgressBar(min = 0, max = length(list_data_not_pa), style = 3)
      progress <- function(n)
        utils::setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      output <-
        foreach::foreach(x = 1:length(list_data_not_pa),
                         .combine = 'c',
                         .options.snow = opts) %d% {
                           
                           if (!parallel)
                             utils::setTxtProgressBar(pb, x)
                           
                           res <- .cell.occupied(
                             size = Resolution,
                             coord = list_data_not_pa[[x]],
                             nbe_rep = nbe_rep
                           )
                           
                           names(res) <- c("spatial", "nbe_occ")
                           res
                         }
      
      if(parallel) snow::stopCluster(cl)
      if(show_progress) close(pb)
      
      loc_not_pa <- unlist(output[names(output) == "nbe_occ"])
      r2 <- unlist(output[names(output) == "spatial"])[[1]]
      names(loc_not_pa) <-
        names(r2) <-
        gsub(pattern = " ",
             replacement = "_",
             names(list_data_not_pa))
      LocOutNatParks[names(LocOutNatParks) %in% names(loc_not_pa)] <-
        loc_not_pa
      
    } else{
      r2 <- NA
    }
    
  }
  
  if (!is.null(protec.areas))
    return(list(r2, r2_PA, LocNatParks, LocOutNatParks))
  if (is.null(protec.areas))
    return(list(r2, Locations))
  
}



#' Internal function
#'
#' Compute IUCN eval
#'
#' @importFrom rgdal project writeOGR
#' @importFrom rnaturalearth ne_countries
#' 
.IUCN.comp <- function(DATA,
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
      ne_countries(scale = 50, returnclass = "sp")
    
    # data('land', package = 'ConR', envir = environment())
    # land <- get("land", envir = environment())
    # data(land, envir = environment())
    poly_borders <- land
  }
  
  ### cropping poly_borders according to range of occurrences shapefile for producing lighter maps
  if (DrawMap) {
    full_poly_borders <- poly_borders
    if (!is.null(poly_borders))
      poly_borders <- raster::crop(poly_borders, raster::extent(MinMax) + 30)
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
    locations.comp(XY = DATA, method = method_locations, 
                 protec.areas = protec.areas, 
                 Cell_size_locations = Cell_size_locations, 
                 ID_shape_PA=ID_shape_PA, show_progress = FALSE)
  
  if(is.null(protec.areas)) {
    r2 <- locations_res[[1]]
    Locations <- locations_res[[2]]
  }else{
    r2 <- locations_res[[1]]
    r2_PA <- locations_res[[2]]
    LocNatParks <- locations_res[[3]]
    LocOutNatParks <- locations_res[[4]]
  }
  
  if(!is.null(protec.areas)) {
    DATA_SF <- as.data.frame(unique(XY))
    colnames(DATA_SF) <- c("ddlon","ddlat")
    sp::coordinates(DATA_SF) <-  ~ddlon+ddlat
    raster::crs(DATA_SF) <- raster::crs(protec.areas)
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
    
    if(is.null(protec.areas)) Results["Nbe_loc",1] <- Locations
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
    
    Results["Category_code",1] <- paste(Results["Category_CriteriaB",1],"B2a")
    if(showWarnings) warning(paste("EOO statistic is not computed for", NamesSp,"because there is less than 3 records"))
  } ## End less than 3 records
  
  ############ Map ###########
  if(DrawMap) {
    
    ## pdf or png format initialization
    if(!map_pdf) {
      if(!is.null(file_name)) {
        NAME_FILE <- paste(file_name, gsub(" ",replacement = "_", as.character(NamesSp)) , sep="")
      }else{
        NAME_FILE <- paste("IUCN_", gsub(" ",replacement = "_", as.character(NamesSp)), sep="")
      }
      FILE_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
      dir.create(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")), showWarnings = FALSE)
    }
    
    if (!map_pdf)
      grDevices::png(
        paste(file.path(paste(
          getwd(), paste("/", FILE_NAME, "_results_map", sep = ""), sep = ""
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
    if(!is.null(protec.areas)){
      if(LocOutNatParks==0){
        graphics::plot(poly_borders, xlim=c(range(XY[,1])[1]-1, range(XY[,1])[2]+1), ylim=c(range(XY[,2])[1]-1, range(XY[,2])[2]+1), axes=FALSE, xlab="", ylab="")
      }else{
        # r2_pol <- r2
        if(LocOutNatParks==1){
          
          graphics::plot(
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
          
          graphics::plot(
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
      graphics::plot(r2, col=rgb(red=1, green=0, blue=0, alpha=0.2), 
           xlim=c(range(XY[,1])[1]-1, range(XY[,1])[2]+1), 
           ylim=c(range(XY[,2])[1]-1, range(XY[,2])[2]+1))
    }
    
    if(SubPop) graphics::plot(SubPopPoly, add=T, border="black", lwd=2, lty=1)
    
    if(!is.null(protec.areas)){
      if(LocNatParks>0){
        if(method_protected_area!="no_more_than_one"){
          # r2_PA_pol <- rasterToPolygons(r2_PA, fun=NULL, n=4, na.rm=TRUE, digits=6, dissolve=FALSE)
          graphics::plot(r2_PA, add=T, col=rgb(red=0, green=0, blue=1, alpha=0.2))
        }
      }
    }
    
    if (!is.null(p1) &
        draw.poly.EOO)
      graphics::plot(p1,
           add = T,
           col = rgb(
             red = 0.2,
             green = 0.2,
             blue = 0.2,
             alpha = 0.1
           ))
    
    graphics::plot(
      poly_borders,
      axes = FALSE,
      lty = 1,
      add = T,
      lwd = 1
    )
    
    if (!is.null(protec.areas))
      graphics::plot(
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
        graphics::plot(
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
        graphics::plot(
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
      graphics::plot(
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
    raster::scalebar(scaleBAR, type="bar", below="kilometres", cex=2.5)
    
    mtext(NamesSp, side=3, cex=3, line=3)
    if(any(colnames(DATA)=="higher.tax.rank")) mtext(DATA[which(DATA[,3]==NamesSp),"higher.tax.rank"][1], side=3, cex=3, line=0.4)
    
    if(add.legend) {
      graphics::par(mar=c(1,1,1,1), xpd=T)
      graphics::plot(1:10, 1:10, type="n", bty='n', xaxt='n', yaxt='n')
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
      graphics::plot(full_poly_borders, lty=1, lwd=1,axes=FALSE)
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
        graphics::plot(
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
    
    dir.create(file.path(paste(getwd(),"/shapesIUCN", sep="")), showWarnings = FALSE)
    
    if (!is.null(p1)) {
      if (length(list.files(paste(getwd(), "/shapesIUCN", sep = ""))) > 0) {
        if (length(grep(paste(NamesSp, "_EOO_poly", sep = ""), unique(sub(
          "....$", '', list.files(paste(getwd(), "/shapesIUCN", sep = ""))
        )))) > 0) {
          FILES <-
            list.files(paste(getwd(), "/shapesIUCN", sep = ""), full.names = TRUE)
          file.remove(FILES[grep(paste(NamesSp, "_EOO_poly", sep = ""), FILES)])
        }
      }
      NAME <- names(p1)
      p1@polygons[[1]]@ID <- "1"
      ConvexHulls_poly_dataframe <- sp::SpatialPolygonsDataFrame(p1, data=as.data.frame(names(p1)))
      colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = '')
      rgdal::writeOGR(ConvexHulls_poly_dataframe,"shapesIUCN",paste(NamesSp,"_EOO_poly", sep=""),driver="ESRI Shapefile")
    }
    
    if(SubPop) {
      if(length(list.files(paste(getwd(),"/shapesIUCN", sep="")))>0){
        if(length(grep(paste(NamesSp,"_subpop_poly", sep=""), unique(sub("....$", '', list.files(paste(getwd(),"/shapesIUCN", sep=""))))))>0) {
          FILES <- list.files(paste(getwd(),"/shapesIUCN", sep=""), full.names = TRUE)
          file.remove(FILES[grep(paste(NamesSp,"_subpop_poly", sep=""), FILES)])
        }
      }
      NAME <- names(SubPopPoly)
      SubPopPoly@polygons[[1]]@ID <- "1"
      ConvexHulls_poly_dataframe <- sp::SpatialPolygonsDataFrame(SubPopPoly, data=as.data.frame(names(SubPopPoly)))
      colnames(ConvexHulls_poly_dataframe@data) <- paste(substr(unlist(strsplit(NamesSp, "[ ]")), 0, 3), collapse = '')
      rgdal::writeOGR(ConvexHulls_poly_dataframe,"shapesIUCN",paste(NamesSp,"_subpop_poly", sep=""),driver="ESRI Shapefile")      
    }
  }
  
  if(SubPop) {
    OUTPUT <- list(Results, p1, SubPopPoly)
    names(OUTPUT) <- c("Results","spatialPoly_EOO","spatialPoly_subpop")
  } else{
    OUTPUT <- list(Results, p1)
    names(OUTPUT) <- c("Results","spatialPoly_EOO")    
  }

  return(OUTPUT)
}

#' Preliminary conservation status assessment following IUCN Criterion B
#' 
#' Given a dataframe of georeferenced occurrences of one, or more, taxa, this
#' function provide statistics values (Extent of Occurrence, Area of Occupancy,
#' number of locations, number of subpopulations) and provide a preliminary
#' conservation status following Criterion B of IUCN.  A graphical map output
#' is also available.
#' 
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
#' \code{coly} and \code{family} are optinal fields
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
#' @param Cell_size_AOO a numeric, value indicating the grid size in kilometres
#' used for estimating Area of Occupancy.  By default, equal to 2
#' @param Cell_size_locations a numeric, value indicating the grid size in
#' kilometres used for estimating the number of location. By default, equal to
#' 10
#' @param Resol_sub_pop a numeric, value indicating the radius size in
#' kilometres used for estimating the number of sub-population. By default,
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
#' @param showWarnings a logical. Wether R should report warnings
#' @param parallel a logical. Wether running in parallel. By default, it is
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
#' \tabular{cccccc}{ [,1] \tab EOO \tab numeric, EOO (square kilometres)\cr
#' [,2] \tab AOO \tab numeric, AOO (square kilometres)\cr [,3] \tab
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
    
    country_map <- rgeos::gBuffer(country_map, byid=TRUE, width=0)
    
  }
  
  if (!is.null(protec.areas)) {
    if (!sp::identicalCRS(protec.areas, country_map))
      raster::crs(protec.areas) <- raster::crs(country_map)
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
    FILE_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
    dir.create(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")), showWarnings = FALSE)
    
    pdf(paste(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep=""),"/","results.pdf", sep=""), width=25, height=25)
  }
  
  if(parallel) {
    # if("doParallel" %in% 
    #    rownames(installed.packages()) == FALSE) {stop("Please install doParallel package")}
    # 
    # library(doParallel)
    
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    # registerDoParallel(NbeCores)
    message('doParallel running with ',
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
                         .IUCN.comp(DATA = list_data[[x]],
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


#' Internal function
#'
#' Compute prop and nbr taxa per cell
#'
#' 
.prop_threat <- function(Cell_count, threshold) {
  NbeRec <- nrow(Cell_count)
  if(NbeRec >= threshold) {
    NbeEsp <- length(unique(Cell_count$tax))
    NbeThreatened <- length(unique(Cell_count[which(Cell_count$Category_CriteriaB %in% c("CR","EN","VU")),"tax"]))
    PropThreatened <- round(NbeThreatened/NbeEsp*100,1)
  }else{
    NbeEsp <- NbeThreatened <- PropThreatened <- NA
  }
  c(NbeRec, NbeEsp, NbeThreatened, PropThreatened)
}



#' Mapping in grid cell results of IUCN.eval
#' 
#' Provides four maps showing in grid cells of a given resolution : number of
#' records, species richness, number of threatened species (CR+EN+VU) and
#' proportion of threatened species. Based on \code{\link[fields]{quilt.plot}}.
#' 
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
#' @importFrom graphics axis box layout legend mtext par plot points
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
  raster::crs(DATA_SF) <- raster::crs(country_map)
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
    FILE_NAME <- ifelse(!is.null(file_name), file_name, "IUCN_")
    dir.create(file.path(paste(
      getwd(), paste("/", FILE_NAME, "_results_map", sep = ""), sep = ""
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
  
  if (!export_map) graphics::par(mfrow=c(2,2))
  if (export_map) graphics::par(mfrow=c(1,1))
  if (export_map) grDevices::png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","number_of_records",".png", sep=""),width=20, height=20, units="cm",res=150)
  
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
        getwd(), paste("/", FILE_NAME, "_results_map", sep = ""), sep = ""
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
  
  if (export_map) grDevices::png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","number_threatened_sp",".png", sep=""),width=20, height=20, units="cm",res=150)
  if (export_map) graphics::par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))

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
  
  if (export_map) grDevices::png(paste(file.path(paste(getwd(),paste("/",FILE_NAME,"_results_map", sep=""), sep="")),"/","proportion_threatened_sp",".png", sep=""),width=20, height=20, units="cm",res=150)
  if (export_map) graphics::par(mar=c(2,2,1,5), las=1, omi=c(0.3,0.4,0.3,0.1))

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


