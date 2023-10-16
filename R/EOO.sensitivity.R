#' @title Sensitivity of Extent of Occurrences
#' 
#' @description Compute changes in the extent of occurrences (EOO) in square
#'   kilometers using occurrence data with different confidence level and
#'   quantifies indirectly the most influential occurrences related to EOO
#'   estimation.
#' 
#' @details 
#'
#' **Input** as a [data.frame][base::data.frame()] should have the following structure:
#' 
#' **It is mandatory to respect field positions, but field names do not
#' matter**
#' 
#' \enumerate{
#'   \item The first column is contains numeric value i.e. latitude in decimal degrees
#'   \item The second column is contains numeric value i.e. longitude in decimal degrees
#'   \item The third column is contains character value i.e. the names of the species
#' }
#' 
#' The function works with a minimum of two classes, that can be, for example, levels of
#' confidence in the geographical coordinates, the taxonomic determination or a
#' period of collection year. For each confidence level, the EOO is computed for
#' a more restricted and reliable set of records. It is essential that the
#' classes of confidence level provided using the argument `levels.order`.
#' For instance, these classes can simply be `c("low", "high")` or `c(FALSE,
#' TRUE)`. But they should match the classes provided in `XY`.  
#' 
#' If argument `occ.based` is `TRUE` (default), the function calculates and
#' return a measure of influence of each occurrence on the estimation of EOO.
#' This measure is simply the distance of each occurrence to the convex hull
#' polygon (i.e. EOO) obtained using only with the occurrences in the highest level
#' of confidence declared, divided by the 95% quantile of the pair-wise
#' distances of the occurrences within the high confidence EOO. Thus, this is
#' only an indirect measure of how much the EOO should increase if a given
#' occurrence was included, because in practice the difference in EOO is not
#' calculated for subsets with and without each occurrence (i.e leave-one-out
#' cross validation).   
#' 
#' The measure is thus a proportion and the highest the value, the more distant
#' the occurrence is from the high confidence EOO. Zeros means that the
#' occurrence is within the high confidence EOO and `NA` means that the high
#' confidence EOO could not be obtained (e.g. less than 2 or 3 high confidence
#' occurrences).
#' Note that even if the EOO cannot be calculated or if coordinates were
#' missing, the function returns a zero for all occurrences with the highest
#' confidence level, which would be by definition within the high confidence
#' EOO.
#' 
#' It is up to the user to define the most appropriate threshold to include or
#' exclude occurrences. A preliminary assessment assuming a circular high
#' confidence EOO suggest that values of 0.5, 1 and 2 would lead to increases of
#' about 25, 50 and 100% in EOO.
#' 
#' As mentioned above, EOO will only be computed if there is at least three
#' unique occurrences in the high confidence level class, unless
#' `method.less.than3` is put to "arbitrary". In that specific case, EOO
#' for species with two unique occurrences will be equal to Dist*Dist*0.1 where
#' Dist is the distance in kilometers separating the two points.
#' 
#' For the very specific (and infrequent) case where all occurrences are
#' localized on a straight line (in which case EOO would be null), 
#' 'noises' are added to coordinates. There is a warning when this happens.
#' 
#' **Notes on computational time** The processing time depends on several
#' factors, including the total number of occurrences, number of confidence
#' levels rpovided and user's computer specifications. Using argument `parallel`
#' equals `TRUE`, greatly increase the processing time, but the processing of
#' large data sets (millions of occurrences) may take hours. On a Intel Core i5,
#' CPU 1.70GHz, 64-bit OS and 16 GB RAM it took 20 min to process about 800
#' thousand records from ~5100 species using 5 cores.
#' 
#' @return A data frame containing, for each taxon, the EOO in square kilometres
#'   for each level of confidence or a list containing this same data frame and
#'   the input data with a new column for each confidence level with a measure 
#'   of influence of the occurrences on the estimation of EOO.
#' 
#' @param XY [data.frame][base::data.frame()] see Details
#' @param file.name a character string. Name file for exported results in csv
#' file. By default is "EOO.sensitivity.results"
#' @param levels.order a character vector with at least two classes ordered from
#'   the least confident to the more confident class of records. See Details.
#' @param occ.based logical. Should the measure of influence of each record be returned? Default to TRUE.
#' @inheritParams EOO.computing
#' @inheritParams over.valid.poly
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
#' mydf <- data.frame(ddlat = c(-44.6,-46.2,-45.4,-42.2,-43.7,-45.0,-28.0,-44.0,
#'                              -26.0,-34.0,-22.0,-40.0,-46.0,-34.0),
#'                    ddlon = c(-42.2,-42.6,-45.3,-42.5,-42.3,-39.0,-17.2,-22.0,
#'                              -46.0,-23.0,-25.0,-43.0,-32.0,-32.0),
#'                    tax = rep(c("a", "b"), each=7),
#'                    valid = c(c(TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,FALSE),
#'                               c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE)))
#' plot(mydf[,2:1], col = as.factor(mydf$tax), 
#'            pch = as.double(as.factor(mydf$tax)))
#' points(mydf[mydf$valid,2:1], col = as.factor(mydf$tax)[mydf$valid], 
#'            pch = 15 + as.double(as.factor(mydf$tax))[mydf$valid])                                 
#' EOO.sensitivity(mydf, levels.order = c(FALSE, TRUE))
#'                    
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom data.table setDT data.table
#' 
#' @export EOO.sensitivity
#' 
EOO.sensitivity <- function(XY,
                          levels.order = NULL,
                          occ.based = TRUE,
                          value = "dist",
                          exclude.area = FALSE,
                          country_map = NULL,
                          alpha = 1,
                          buff.alpha = 0.1,
                          method.range = "convex.hull",
                          #Name_Sp = "species1",
                          method.less.than3 = "not comp",
                          file.name = "EOO.sensitivity.results",
                          parallel = FALSE,
                          NbeCores = 2,
                          show_progress = TRUE,
                          proj_type = "cea",
                          mode = "spheroid"
){ 
  
  XY$recordID <- 1:dim(XY)[1]
  
  if(length(unique(as.data.frame(XY)[,4])) < 2)
    stop("there is only one class of confidence level")
  
  if (ncol(XY) > 4) {
    colnames(XY)[1:4] <- c("ddlat", "ddlon", "tax", "valid")
    XY$valid <- as.character(XY$valid)
    XY$tax <- as.character(XY$tax)
    XY <- XY[, c("ddlat", "ddlon", "tax", "valid", "recordID")]
  } else{
    colnames(XY)[1:3] <- c("ddlat", "ddlon", "valid")
    XY$valid <- as.character(XY$valid)
    XY$tax <- "Species 1"
    XY <- XY[, c("ddlat", "ddlon", "tax", "valid", "recordID")]
  }
  
  
  # levels.order_in_data <- levels.order[levels.order %in% unique(XY$valid)]
  levels.order_in_data <- unique(XY$valid)[unique(XY$valid) %in% levels.order]
  
  if(length(levels.order_in_data) != length(levels.order))
    stop("The number of element of levels.order does not match the number classes of confidence level")
  
  levels.order <- levels.order_in_data
  XY <- XY[XY$valid %in% levels.order, ]
  n.levels <- length(levels.order)
  
  
  ## Obtaining the data for each class of confidence level
  XY$classes <- 
    # as.double(factor(XY$valid, levels = levels.order, labels = 1:n.levels))
    as.double(factor(XY$valid, levels = levels.order_in_data, labels = 1:n.levels))
  
  XY.orig <- XY

  XY.list <- 
    coord.check(XY = XY, 
                listing = TRUE, listing_by_valid = TRUE)$list_data
  
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
  
  names_ <- names(XY.list)
  
  x <- NULL
  if (show_progress) {
    pb <-
      txtProgressBar(min = 0,
                     max = length(XY.list),
                     style = 3)
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
  }
  
  output <-
    foreach::foreach(
      x = 1:length(XY.list),
      .combine = 'c',
      .options.snow = opts
    ) %d% {

      if (!parallel & show_progress)
        setTxtProgressBar(pb, x)
      
      res <-
        ConR::EOO.comp(
          XY = XY.list[[x]],
          exclude.area = exclude.area,
          country_map = country_map,
          #Name_Sp = names_[x],
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
    data.frame(EOO = unlist(output[grep("EOO", names(output))]),
               nbe_occ = unlist(lapply(XY.list, nrow)),
               tax = unlist(lapply(strsplit(names_, split = "___"), function(x) x[1])),
               classes = as.numeric(unlist(lapply(strsplit(names_, split = "_"), function(x) x[length(x)]))))
  row.names(Results_short) <- names_
  
  cat("EOO differences...", sep= "\n")

  data.table::setDT(Results_short)
  increases <-
    Results_short[, lapply(.SD, function(x)
      (x[2:length(x)] - x[1]) / x[1] * 100),
      by = .(tax), .SDcols = "EOO"]
  
  comp_names <- 
    Results_short[, lapply(.SD, function(x) x[2:length(x)]), 
                  by = .(tax), .SDcols = "classes"]
  
  increases <- 
    data.table::data.table(increases, classes = comp_names$classes)
  
  colnames(increases)[colnames(increases) == "EOO"] <- 
    "EOO.increase"
  
  Results_short <- 
    merge(Results_short, increases, by = c("tax", "classes"), all.x = TRUE)
  
  #   colnames(increases)[2] <-
  #     paste("EOO", 2, sep = "_")
  #
  #
  #   Results_short_subst <-
  #     Results_short[Results_short$classes == 1, c("EOO", "nbe_occ", "tax")]
  #
  #   Results_short_subst <-
  #     merge(Results_short_subst, increases, by = 'tax')
  #
  # print(increases)
  # print(comp_names)
  
  # if (length(output) == 1)
  #   names(output) <- Name_Sp


  if(occ.based) {
    # cat("extract spatial")
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
    
    ### convert into sf, may be slow for large dataset, see if necessary step
    output_spatial <- 
      st_as_sf(data.frame(output_spatial, name = names_[id_spatial]))
    
  }
  
  if(occ.based) {
    cat("Starting the occurrence-based analysis...", sep= "\n")

    XY.list.taxa <- 
      coord.check(XY = XY, listing = TRUE)$list_data
    
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
                       max = length(XY.list.taxa),
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
        x = 1:length(XY.list.taxa),
        .combine = 'rbind',
        .options.snow = opts
      ) %d% {

        if (!parallel & show_progress)
         setTxtProgressBar(pb, x)
        
        res <-
          over.valid.poly(
            poly = output_spatial,
            points = XY.list.taxa[[x]],
            names_taxa = names(XY.list.taxa)[x],
            value = value,
            names_poly = names_[id_spatial],
            min.dist = 0.1
          )
        res
      }
    
    if(parallel) snow::stopCluster(cl)
    if(show_progress) close(pb)
    
    output <- 
      output[!is.na(output[,1]),]
    
    # XY.list[[1]]$prop.dist.eoo <- output
    # Results_long <- dplyr::left_join(XY,
    #                                  XY.list[[1]][, c("recordID", "classes","prop.dist.eoo")],
    #                                  by = "recordID")
    # Results_long$prop.dist.eoo[is.na(Results_long$prop.dist.eoo) &
    #                              Results_long$classes >= max(Results_long$classes, na.rm = TRUE)] <- 0
    # Results_long <- Results_long[order(Results_long$recordID), ]
    # Results_long <- 
    #   Results_long[, -which(names(Results_long) %in% c("recordID", "classes"))]
  }
  
  # if (write_results)
  #   write.csv(Results_short, paste(getwd(), "/", file.name, ".csv", sep = ""))
  
  if (occ.based) {
    
    output <- merge(XY.orig, output[ ,c("recordID", "rel_dist")], 
                    by = "recordID", all.x = TRUE)
    output <- output[order(output$recordID), ]
    
    output <- list(results = as.data.frame(Results_short),
                   results_occ = output,
                   spatial = output_spatial)
  } else {
    output <- Results_short
  }
  
  cat("Returning the results.", sep= "\n")
  
  return(output)
}

