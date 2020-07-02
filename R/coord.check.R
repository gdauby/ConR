
#' Coordinates check
#'
#' @param XY 
#' @param listing 
#' @param proj_type
#'
#' @return
#'
coord.check <- function(XY, listing = TRUE, proj_type = NULL) {
  XY <- as.data.frame(XY)
  
  if (any(is.na(XY[, c(1:2)]))) {
    print(
      paste(
        "Skipping",
        length(which(rowMeans(is.na(
          XY[, 1:2]
        )) > 0)) ,
        "occurrences because of missing coordinates for",
        # if(verbose)
        length(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                 0), 3])), "taxa"
      )
    )
    XY <- XY[which(!is.na(XY[, 1])), ]
    XY <- XY[which(!is.na(XY[, 2])), ]
  }
  
  if (any(XY[, 2] > 180) ||
      any(XY[, 2] < -180) ||
      any(XY[, 1] < -180) ||
      any(XY[, 1] > 180))
    stop("coordinates are outside of expected range")
  
  if (!is.null(proj_type)) {
    
    XY_proj <-
      sf::sf_project(
        from = sf::st_crs(4326),
        to =
          sf::st_crs(proj_type),
        pts = XY[, c(2, 1)]
      )[, c(2, 1)]
    
    XY[, c(1, 2)] <- 
      XY_proj[, c(1, 2)]
    
  }
  
  
  if (listing) {
    if (ncol(XY) > 2) {
      colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
      XY$tax <- as.character(XY$tax)
      
      if(length(grep("[?]", XY[,3]))>0) XY[,3] <- gsub("[?]", "_", XY[,3])
      if(length(grep("[/]", XY[,3]))>0) XY[,3] <- gsub("[/]", "_", XY[,3])
      
      list_data <- split(XY, f = XY$tax)
    } else{
      colnames(XY)[1:2] <- c("ddlat", "ddlon")
      list_data <- list(XY)
    }
  } else{
    list_data <-
      XY
  }
  
  return(list_data)
}