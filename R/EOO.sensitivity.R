#' @title Sensitiity of Extent of Occurrences
#' 
#' @description Compute changes in the extent of occurrences (EOO) in square
#'   kilometers using occurrence data with different confidence level and
#'   quantifies indirectly the most influentional occurrences related to EOO
#'   estimation.
#' 
#' @details 
#'
#' \strong{Input} as a \code{dataframe} should have the following structure:
#' 
#' \strong{It is mandatory to respect field positions, but field names do not
#' matter}
#' 
#' \tabular{ccc}{ [,1] \tab ddlat \tab numeric, latitude (in decimal
#' degrees)\cr [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#' [,3] \tab tax \tab character or factor, taxa names\cr 
#' [,4] \tab valid \tab character or factor, classes of confidence level}
#' 
#' The function works with a minimum of two classes, that can be levels of
#' confidence in the geographical coordinate, the species determination or a
#' period of collection year. For each confidence level, the EOO is computed for
#' a more restricted and reliable set of records. It is essential that the
#' classes of of confidence level provided using the argument `levels.order`.
#' For instance, these classes can simply be `c("low", "high")` or `c(FALSE,
#' TRUE)`. But they should match the classes provided in \code{XY}.  
#' 
#' If argument `occ.based` is `TRUE` (the default), the function calculates and
#' return a measure of influence of each occurrence on the estimation of EOO.
#' This measure simply is the distance of each occurrence to the convex hull
#' polygon (i.e. EOO) obtained using only the occurrences in the highest level
#' of confidence declared, divided by the 95% quantile of the pair-wise
#' distances of the occurrences within the high confidence EOO. Thus, this is
#' only an indirect measure of how much the EOO should increase if a given
#' occurrence was included, beacuse in practice the difference in EOO is not
#' calculated for subsets with and without each occurrence (i.e leave-one-out
#' cross validation).   
#' 
#' The measure is thus a proportion and the highest the value, the more distant
#' the occurrence is from the high confidence EOO. Zeros means that the
#' occurrence is within the high confidence EOO and `NA` means that the high
#' confidence EOO could not be obtained (i.e. less then 2 or 3 high confidence
#' occurrences) or that latitude or longitude was missing for the occurrence.
#' Note that even if the EOO cannot not be calculated or if coordinates were
#' missing, the function returns a zero for all occurrences with the highest
#' confidence level, which would be by definition within the high confidence
#' EOO.
#' 
#' It is up to the user to define the most appropriate threshold to include or
#' exclude occurrences. A preliminary assessment assuming a circular hgh
#' confidence EOO suggest that values of 0.5, 1 and 2 would lead to increases of
#' about 25, 50 and 100% in EOO.
#' 
#' As mentioned above, EOO will only be computed if there is at least three
#' unique occurrences in the high confidence level class, unless
#' \code{method.less.than3} is put to "arbitrary". In that specific case, EOO
#' for species with two unique occurrences will be equal to Dist*Dist*0.1 where
#' Dist is the distance in kilometers separating the two points.
#' 
#' For the very specific (and infrequent) case where all occurrences are
#' localized on a straight line (in which case EOO would be null), EOO is
#' estimated by the area of polygon surrounding this straight line with a
#' buffer of \code{buff.alpha} decimal degree. There is a warning when this
#' happen.
#' 
#' \strong{Notes on computational time} The processing time depends on several
#' factors, including the total number of occurrences, number of confidence
#' levels rpovided and user's computer specifications. Using argument `parellel`
#' equals `TRUE`, greatly increase the processing time, but the processing of
#' large data sets (milions of occurrences) may take hours. On a Intel Core i5,
#' CPU 1.70GHz, 64-bit OS and 16 GB RAM it took 20 min to process about 800
#' thousand records from ~5100 species using 5 cores.
#' 
#' @return A data frame containing, for each taxon, the EOO in square kilometers
#'   for each level of confidence or a list containing this same data frame and
#'   the input data with a new column for each confidence level with a measure 
#'   of influence of the occurrences on the estimation of EOO.
#' 
#' @param XY \code{dataframe} see Details
#' @param file.name a character string. Name file for exported results in csv
#' file. By default is "EOO.sensitivity.results"
#' @param levels.order a character vector with at least two classes ordered from
#'   the least confident to the more confident class of records. See Details.
#' @param occ.based logical. Should the measure of influence of each record be returned? Default to TRUE.
#' @inheritParams EOO.computing
#' @inheritParams .over.valid.poly
#' 
#' @author Renato A. Ferreira de Lima & Gilles Dauby
#' 
#' @references 
#' Nic Lughadha, Staggemeier, Vasconcelos, Walker, Canteiro & Lucas (2019).
#' Harnessing the potential of integrated systematics for conservation of
#' taxonomically complex, megadiverse plant groups. Conservation Biology, 33(3):
#' 511-522.
#' 
#' Gaston & Fuller (2009). The sizes of species' geographic ranges.
#' Journal of Applied Ecology, 46(1): 1-9.
#' 
#' @examples 
#' 
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-28.0,-44.0,-26.0,-34.0,-22.0,-40.0,-46.0,-34.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-17.2,-22.0,-46.0,-23.0,-25.0,-43.0,-32.0,-32.0),
#'                    tax = rep(c("a", "b"), each=7),
#'                    valid = c(c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE),c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE)),
#'                    stringsAsFactors = FALSE)
#' plot(mydf[,2:1], col = as.factor(mydf$tax), 
#'            pch = as.double(as.factor(mydf$tax)))
#' points(mydf[mydf$valid,2:1], col = as.factor(mydf$tax)[mydf$valid], 
#'            pch = 15 + as.double(as.factor(mydf$tax))[mydf$valid])                                 
#' EOO.sensitivity(mydf, levels.order = c(FALSE, TRUE), 
#'            proj_user = 5641)
#'                    
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom dplyr left_join
#' @importFrom rgeos gBuffer
#' @importFrom methods slot
#' @importFrom sp SpatialPolygonsDataFrame
#' 
#' @export EOO.sensitivity
#' 
EOO.sensitivity <- function(XY,
                          levels.order = NULL,
                          occ.based = TRUE,
                          value = "dist",
                          proj_user = NULL,
                          exclude.area = FALSE,
                          country_map = NULL,
                          alpha = 1,
                          buff.alpha = 0.1,
                          method.range = "convex.hull",
                          #Name_Sp = "species1",
                          buff_width = 0.1,
                          method.less.than3 = "not comp",
                          write_results = FALSE,
                          file.name = "EOO.sensitivity.results",
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE
){ 
  
  XY$recordID = 1:dim(XY)[1]
  
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
    XY1 <- XY[which(!is.na(XY[, 1])), ]
    XY1 <- XY1[which(!is.na(XY1[, 2])), ]
    XY1 <- as.data.frame(XY1)
  } else {
    XY1 <- as.data.frame(XY)
  }
  
  if (exclude.area & is.null(country_map))
    stop("exclude.area is TRUE but no country_map is provided")
  
  if (buff_width > 80)
    stop("buff_width has an unrealistic value")
  
  if (any(XY1[, 2] > 180) ||
      any(XY1[, 2] < -180) ||
      any(XY1[, 1] < -180) ||
      any(XY1[, 1] > 180))
    stop("coordinates are outside of the expected range")
  
  if(length(unique(XY1[,4])) <2)
    stop("there is only one class of confidence level")
  
  if(!is.null(country_map))
    country_map <- suppressWarnings(rgeos::gBuffer(country_map, byid=TRUE, width=0))
  
  if (method.range == "convex.hull") {
    convex.hull = TRUE
    alpha.hull = FALSE
  }
  
  if (method.range == "alpha.hull") {
    convex.hull = FALSE
    alpha.hull = TRUE
  }

  if (ncol(XY1) > 4) {
    colnames(XY1)[1:4] <- c("ddlat", "ddlon", "tax", "valid")
    XY1$valid <- as.character(XY1$valid)
    XY1$tax <- as.character(XY1$tax)
    XY1 <- XY1[, c("ddlat", "ddlon", "tax", "valid", "recordID")]
  } else{
    colnames(XY1)[1:3] <- c("ddlat", "ddlon", "valid")
    XY1$valid <- as.character(XY1$valid)
    XY1$tax <- "Species 1"
    XY1 <- XY1[, c("ddlat", "ddlon", "tax", "valid", "recordID")]
  }
  
  levels.order <- levels.order[levels.order %in% unique(XY1$valid)]
  XY1 <- XY1[XY1$valid %in% levels.order, ]
  n.levels <- length(levels.order)
  
  if(occ.based) {
    export_shp = c(rep(FALSE, n.levels - 1), TRUE)
  } else {
    export_shp = c(rep(FALSE, n.levels - 1), FALSE) 
  }
  
  ## Obtaining the data for each class of confidence level
  XY1$classes <- as.double(factor(XY1$valid, levels = levels.order, labels = 1:n.levels))

  XY.list <- sp_names <- vector("list", n.levels)
  for(i in 1:n.levels) {
    tmp <- XY1[XY1$classes >= i, ]
    XY.list[[i]] <- tmp
    sp_names[[i]] <- sort(unique(tmp$tax))
    rm(tmp)
  }
  names(XY.list) <- names(sp_names) <- paste0("level.", 1:n.levels)
  
  ## Obtaining EOO for each taxon and class of confidence level
  result <- vector("list", n.levels) 
  names(result) <- paste0("level.", 1:n.levels)  
  cat("Starting the EOO analysis for each species and confidence levels...", sep= "\n")
  for(i in 1:length(result)) {
    result[[i]] <- EOO.computing(XY = XY.list[[i]],
                               exclude.area = exclude.area,
                               buff_width = buff_width,
                               country_map = country_map,
                               write_shp = FALSE,
                               #Name_Sp = names_list[[i]],
                               method.range = method.range,
                               alpha = alpha,
                               buff.alpha = buff.alpha,
                               method.less.than3 = method.less.than3,
                               #alpha.hull = alpha.hull,
                               #convex.hull = convex.hull,
                               write_results = write_results,
                               export_shp = export_shp[i],
                               parallel = parallel,
                               NbeCores = NbeCores,
                               show_progress = show_progress)
  }
  
  if(occ.based) {
    eoos <- do.call(rbind.data.frame, result[[n.levels]][grepl("EOO", names(result[[n.levels]]))])
    dimnames(eoos) <- list(sp_names[[n.levels]], "EOO")
    shps <- result[[n.levels]][grepl("spatial.polygon", names(result[[n.levels]]))]
    for (i in 1:length(shps))
      methods::slot(methods::slot(shps[[i]], "polygons")[[1]], "ID") <- sp_names[[n.levels]][!is.na(eoos$EOO)][i]
    shps <- do.call(rbind, shps)
    shps_df <- sp::SpatialPolygonsDataFrame(shps, data.frame(tax = names(shps), row.names = names(shps)))
    shps_df$tax <- as.character(shps_df$tax) 
    result[[n.levels]] <- eoos
    rm(eoos, shps)
  }

  for (i in 1:length(result))
    result[[i]]$n.occs <- as.double(table(XY.list[[i]]$tax)) 
  
  Results_short <- merge(result[[1]], result[[2]], by="row.names", 
                         all = TRUE, suffixes = c(".conf1",".conf2"))
  if(n.levels > 2) {
    for(i in 3:n.levels)
      Results_short <- merge(Results_short, result[[i]], 
                             all = TRUE, by.x="Row.names", by.y="row.names", suffixes = c("", paste0(".conf",i)))
  }
  names(Results_short) <- c("Species", paste0(rep(c("EOO.level.","Occs.level."),n.levels), rep(1:n.levels, each = 2)))

  result1 <- Results_short[, grepl("EOO", names(Results_short))]
  result1 <- do.call(rbind.data.frame,
                     lapply(1:nrow(result1),  
                        function(x) round((result1[x, 1:(n.levels-1)] - result1[x, n.levels]) / 
                            result1[x, n.levels], 5)))
  names(result1) <- paste0("EOO.increase.", 1:(n.levels-1))
  Results_short <- cbind.data.frame(Results_short,
                                    result1, stringsAsFactors = FALSE)
  Results_short <- Results_short[,c(1, grep("Occs", names(Results_short)),
                                    grep("EOO.level", names(Results_short)),
                                    grep("EOO.increase", names(Results_short)))]
  
  if(occ.based) {
    cat("Starting the occurrence-based analysis...", sep= "\n")

    list_data <- split(XY.list[[1]], f = XY.list[[1]]$tax)
    
    if (parallel) {
      cl <- snow::makeSOCKcluster(NbeCores)
      doSNOW::registerDoSNOW(cl)
      
      message('Parallel running with ',
              NbeCores, ' cores')
      
      `%d%` <- foreach::`%dopar%`
    } else {
      `%d%` <- foreach::`%do%`
    }
    
    if(show_progress) {
      pb <- txtProgressBar(min = 0,
                       max = length(list_data),
                       style = 3)
      
      progress <- function(n)
        setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    } else {
      opts <- NULL
    }
    
    x <- NULL
    output <-
      foreach::foreach(
        x = 1:length(list_data),
        .combine = 'c',
        .options.snow = opts
      ) %d% {
        #source("C://Users//renato//Documents//raflima//R_packages//ConR//R//over.valid.poly.R")
        
        if (!parallel & show_progress)
         setTxtProgressBar(pb, x)
        
        res <- .over.valid.poly(shps_df, list_data[[x]], 
                                proj_user = proj_user, value = value)
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    XY.list[[1]]$prop.dist.eoo <- output
    Results_long <- dplyr::left_join(XY,
                                     XY.list[[1]][, c("recordID", "classes","prop.dist.eoo")],
                                     by = "recordID")
    Results_long$prop.dist.eoo[is.na(Results_long$prop.dist.eoo) &
                                 Results_long$classes >= max(Results_long$classes, na.rm = TRUE)] <- 0
    Results_long <- Results_long[order(Results_long$recordID), ]
    Results_long <- 
      Results_long[, -which(names(Results_long) %in% c("recordID", "classes"))]
  }
  
  if (write_results)
    write.csv(Results_short, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (occ.based) {
    output <- list(Results_short,
                   Results_long)
    names(output) = c("EOO.change", "Occ.influence")
  } else {
    output <- Results_short
  }
  
  cat("Returning the results.", sep= "\n") 
  output
}

