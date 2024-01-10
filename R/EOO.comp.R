#' @title Compute Extent of Occurrence
#'
#' @description Compute EOO given a data frame of coordinate in decimals degrees
#'
#' @param XY data.frame
#' @param exclude.area logical. Default if FALSE
#' @param country_map SpatialPolygonDataframe. Default if NULL
#' @param method.range string, by default "convex.hull", can also take
#'   "alpha.hull"
#' @param alpha integer
#' @param buff.alpha numeric
#' @param method.less.than3 string
#' @param mode character string either 'spheroid' or 'planar'. By default
#'   'spheroid'
#' @param proj_type crs
#' @param reproject logical FALSE whether the polygon should be converted to geographic coordinates if mode is `planar`
#'
#' @author Gilles Dauby & Renato A. Ferreira de Lima
#' 
#' @return A list
#' \enumerate{
#'   \item EOO a numeric vector of AOO estimates for each taxa
#'   \item spatial.polygon a simple feature collection
#' }
#' 
#' @keywords internal
#' 
#' @import sf
#' @importFrom stats dist cor
#' @importFrom rnaturalearth ne_countries
#' @export
EOO.comp <-  function(XY,
                      # exclude.area = FALSE,
                      # country_map = NULL,
                      # Name_Sp = "tax",
                      method.range = "convex.hull",
                      alpha = 1,
                      buff.alpha = 0.1,
                      method.less.than3 = "not comp",
                      mode = "spheroid",
                      proj_type = NULL,
                      reproject = FALSE
) {
  
  if (!requireNamespace("lwgeom", quietly = TRUE))
    stop(
      "The 'lwgeom' package is required to run this function. ",
      "Please install it first."
    )
  
  method.range <- match.arg(method.range, c("convex.hull", "alpha.hull"))
  
  # XY <- 
  #   coord.check(XY = XY, listing = FALSE, check_eoo = FALSE)
  
  # XY <- XY$list_data
  
  Name_Sp <- XY[1,3]
  
  # if (exclude.area) {
  #   ### Getting by default land map if poly_borders is not provided
  #   if (is.null(country_map)) {
  #     
  #     country_map <-
  #       rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
  #     
  #   } 
  # }
  
  
  
  ## Check if there are less than 3 unique occurrences
  if (nrow(XY) < 3) {
    ## if there is only one occurrence, EOO is NA
    if (nrow(XY) < 2) {
      
      EOO <- NA
      # message(
      #   paste(
      #     "EOO parameter cannot be estimated for",
      #     as.character(Name_Sp),
      #     "because there is only 1 unique occurrence"
      #   )
      # )
      
    } else {
      if (method.less.than3 == "arbitrary") {
        
        # projEAC <- proj_type
        # 
        # # coordEAC <-
        # #   as.data.frame(matrix(unlist(
        # #     rgdal::project(
        # #       as.matrix(unique(XY)[, 1:2]),
        # #       proj = as.character(projEAC),
        # #       inv = FALSE
        # #     )
        # #   ), ncol = 2))
        # 
        # coordEAC <-
        #   sf_project(
        #   from = sf::st_crs(4326),
        #   to =
        #     sf::st_crs(projEAC),
        #   pts = XY[, c(1, 2)]
        # )
        
        EOO <-
          as.numeric(dist(XY[, c(1, 2)] / 1000) * 0.1 * dist(XY[, c(1, 2)] / 1000))
      }
      
      if (method.less.than3 == "not comp") {
        ## if there are two unique occurences, EOO is not computed neither
        # message(
        #   paste(
        #     "EOO parameter cannot be estimated for",
        #     as.character(Name_Sp),
        #     "because there is less than 3 unique occurrences"
        #   )
        # )
        EOO <- NA
      }
    }
    
    OUTPUT <- list(EOO = EOO, spatial.polygon = NA)
    
  } else {
    
    ### Checking if all occurrences are on a straight line
    if (length(XY[, 1]) == 1 ||
        length(XY[, 2]) == 1 ||
        round(abs(suppressWarnings(cor(XY[, 1], XY[, 2]))), 6) == 1 ||
        is.na(round(abs(suppressWarnings(cor(XY[, 1], XY[, 2]))), 6))) {
      
      message(
        paste(
          "Occurrences of",
          as.character(Name_Sp),
          "follow a straight line, 'noise' to coordinates is added"
        )
      )
      
      check_line <- TRUE
      while(check_line) {
        
        XY[,c(1, 2)] <- 
          apply(XY[,c(1, 2)], 2, function(x) jitter(x, factor = 0.1))
        
        if (round(abs(cor(XY[, 1], XY[, 2])), 6) != 1)
          check_line <- FALSE
        
      }
      
    }
    # else {
      
      if (method.range == "alpha.hull") {
        
        ### work around to avoid bug appearing randomly
        # cont_try <- TRUE
        # while(cont_try) {
          
          p1 <-
            alpha.hull.poly(
              XY = XY[, c(2, 1)],
              alpha = alpha,
              buff = buff.alpha,
              mode = mode,
              proj_type = proj_type)

        #   if (!grepl("trye-error", class(p1)))
        #     cont_try <- FALSE
        # }
        
      }
        
      if (method.range == "convex.hull")
        p1 <-
          Convex.Hull.Poly(XY = XY[, c(2, 1)],
                            mode = mode,
                            proj_type = proj_type)
      
      
      if(any(class(p1) == "SpatialPolygons") | any(class(p1) == "sfc") | any(class(p1) == "sf")) {
        
        if (!sf::st_is_valid(p1)) {
          
          p1 <- sf::st_make_valid(p1)
          
          warning(
            paste(
              "Failed to build a valid polygon to estimate EOO for _",
              as.character(Name_Sp),
              "_ This is probably because occurrences spans more than 180 degrees longitude."
            )
          )
          
          EOO <- NA
          p1 <- NA
          
        } else {
          
          EOO <-
            as.numeric(st_area(p1)) / 1000000
          
          if (mode == "planar" & reproject)
            p1 <-
              sf::st_transform(p1, 4326)
          
          p1$tax <- 
            Name_Sp
          
        }
        
      } else  {
        
        EOO <- NA
        
      }
    # }
    
    OUTPUT <- list(EOO = EOO, spatial.polygon = p1)
  }
  
  digits <-
    c(6, 5, 4, 3, 2, 1, 0)[findInterval(OUTPUT$EOO, 
                                        c(0, 0.0001, 0.01, 0.1, 1, 10, 30000, Inf))]
  
  OUTPUT$EOO <- round(OUTPUT$EOO, digits)
  
  # if(verbose) cat(" ",paste(Name_Sp,"EOO comp."))
  
  return(OUTPUT)
}

