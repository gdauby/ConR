
#' @title Number of Subpopulations
#'
#' @description Estimate the number of subpopulations following the method
#'   **circular buffer method** (overlapping buffered circles form a single subpopulation)
#'
#' @author Gilles Dauby & Renato A. Ferreira de Lima
#'
#' @param XY a data frame containing the geographical coordinates for each taxon
#'   (see Details).
#' @param Resol_sub_pop a value defining the radius of the circles around each
#'   occurrence (in kilometres) or data frame vector containing a column 'tax' 
#'   with the taxa names and a column 'radius' with the species-specific radius
#'   (in kilometre as well). Typically, this data frame is the output of
#'   ```ConR``` function ```subpop.radius```.
#' @param export_shp logical. Whether the resulting shapefiles should be
#'   exported. FALSE by default.
#' @param parallel logical. Whether compute should run in parallel. FALSE by default.
#' @param NbeCores integer. Number of cores for parallel computation. Default to
#'   2.
#' @param show_progress logical. Whether a bar showing progress in computation
#'   should be shown. TRUE by default.
#' @param proj_type character string or numeric or object of CRS class, by
#'   default is "cea"
#' 
#' @details 
#' **Input** as a `dataframe` should have the following structure:
#' 
#' **It is mandatory to respect field positions, but field names do not matter**
#' 
#' \tabular{ccc}{
#'   [,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
#'   [,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
#'   [,3] \tab tax \tab character or factor, taxa names\cr
#' }
#' 
#' @references Rivers MC, Bachman SP, Meagher TR, Lughadha EN, Brummitt NA
#'   (2010) Subpopulations, locations and fragmentation: applying IUCN red list
#'   criteria to herbarium specimen data. Biodiversity and Conservation 19:
#'   2071-2085. doi: 10.1007/s10531-010-9826-9
#'
#' @return 
#' If `export_shp` is TRUE, a list with [[1]]number_subpop and [[2]]poly_subpop a `Simple feature collection` with as many MULTIPOLYGON as taxa
#' If `export_shp` is FALSE, a vector with estimated number of subpopulation per taxa
#' 
#' @examples 
#' data(dataset.ex)
#'
#' subpop.comp(dataset.ex, Resol_sub_pop = 5)
#' rad.df <- data.frame(
#'     tax = unique(dataset.ex$tax),
#'     radius = seq(3,13, by=2),
#'     stringsAsFactors = FALSE
#'   )
#' subpop.comp(dataset.ex, Resol_sub_pop = rad.df)
#' subpop.comp(dataset.ex, Resol_sub_pop = rad.df, export_shp = TRUE)
#' 
#' @importFrom snow makeSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom foreach foreach
#' 
#' 
#' @export subpop.comp
#' 
subpop.comp <- function(XY, 
                        Resol_sub_pop = NULL, 
                        proj_type = "cea",
                        export_shp = FALSE,
                        parallel = FALSE,
                        show_progress = TRUE,
                        NbeCores = 2) {
  
  ### ADDED BY RENATO: SPC RECOMMEND NOT USING ANY DEFAULTS TO FORCE ASSESSORS TO THINK ###
  if (is.null(Resol_sub_pop)) 
    stop("Radius is missing, please provide a value for all species or a data frame with species-specific values")
  
  proj_type <- 
    proj_crs(proj_type = proj_type)
  
  ### PART INCLUDED BY RENATO ###
  if ("data.frame" %in% class(Resol_sub_pop)) {
    XY <- merge(XY, Resol_sub_pop, 
                by = "tax", all.X = TRUE, sort = FALSE)
    XY <- XY[,c("ddlat", "ddlon", "tax", "radius")]
  } else {
    XY$radius <- Resol_sub_pop   
  }
  ### END OF PART INCLUDED BY RENATO ###
  
  list_data <-
    coord.check(XY = XY, listing = TRUE, proj_type = proj_type)
  
  if (parallel) {
    cl <- snow::makeSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else{
    `%d%` <- foreach::`%do%`
  }
  
  x <- NULL
  
  if (show_progress) {
    pb <-
      utils::txtProgressBar(min = 0,
                            max = length(list_data[[1]]),
                            style = 3)
    
    progress <- function(n)
      utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else{
    opts <- NULL
  }
  
  output <-
    foreach::foreach(
      x = 1:length(list_data[[1]]),
      .combine = 'c',
      .options.snow = opts
    ) %d% {
      
      if (!parallel & show_progress)
        utils::setTxtProgressBar(pb, x)
      
      res <- 
        subpop.estimation(
          XY = list_data[[1]][[x]], 
          Resol_sub_pop = unique(list_data[[1]][[x]]$radius), #### PART EDITED BY RENATO #### 
          proj_type = proj_type,
          export_shp = export_shp #### NEW RGUMENT ADDED BY RENATO ####
        )
      
      if (export_shp) {
        names(res) <- c("subpop", "spatial")
        names(res)[1] <- list_data[[1]][[x]]$tax[1]
        res$spatial <- cbind(res$spatial, tax = list_data[[1]][[x]]$tax[1])
        
      } else {
        
        names(res)[1] <- list_data[[1]][[x]]$tax[1]
        
      }
      
      res
    }
  
  if(parallel) snow::stopCluster(cl)
  if(show_progress) close(pb)
  
  if (export_shp) {
    
    number_subpop <-
      data.frame(subpop =  unlist(output[names(output) != "spatial"]))

    shapes <- output[names(output) == "spatial"]
    shapes <- do.call('rbind', shapes)
    row.names(shapes) <- 1:nrow(shapes)
    
  } else {
    
    number_subpop <-
      data.frame(subpop =  unlist(output))
  }
  

  
  ### GILLES, NOT SURE WHY IT IS NECESSARY TO TRANFORM SPECIES ORIGINAL NAMES...
  #SO I CHANGED IT, BUT LEFT THE PREVIOUS CODE IF YOU WANT TO TAKE IT BACK
  # SpNames <- gsub(pattern = " ", 
  #                 replacement = "_", 
  #                 names(list_data))
  # SpNames <- names(list_data)
  # names(number_subpop) <- SpNames
  # 
  # if (export_shp) { ## IF/ELSE ADDED BY RENATO
  #   poly <- 
  #     output[names(output) == "poly_subpop"]
  #   names(poly) <- SpNames

    ### GILLES: I INCLUDE THIS PART FROM ANOTHER FUNCTION, SINCE NOW
    #THE OUTPUT 'subpop.estimation' ARE sf OBJECTS WITH MULTIPLE POLYGONS/CIRCLES
  #   if(length(poly) > 1) {
  #     poly <-
  #       do.call("rbind", poly)
  #     row.names(poly) <- NULL
  #     #### GILLES: MAYBE RETURN THE POLYGONS IN THE SAME CRS OF THE OCCURRENCES: WSG84?   
  #     # poly <- 
  #     #   sf::st_transform(poly, crs = 4326)
  #     
  #   } else {
  #     poly <-
  #       poly[[1]]
  #     poly <- 
  #       sf::st_as_sf(data.frame(poly, tax = SpNames[1]))
  #   }
  # 
  #   # if (length(OUTPUT) == 1)
  #   #   OUTPUT <- OUTPUT[[1]]
  # 
  #   OUTPUT <- list(number_subpop = number_subpop, poly_subpop = poly)
  # 
  # } else {
  #   
  #   OUTPUT <- number_subpop
  #   
  # }  
  
  if (!export_shp) return(number_subpop)
  
  if (export_shp) return(list(number_subpop = number_subpop,
                              poly_subpop = shapes))
  
}
