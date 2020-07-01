
#' Coordinates check
#'
#' @param XY 
#' @param listing 
#'
#' @return
#'
coord.check <- function(XY, listing = TRUE) {
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
        paste(as.character(unique(XY[which(rowMeans(is.na(XY[, 1:2])) >
                                             0), 3])), collapse = " AND ")
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
  
  if (listing) {
    if (ncol(XY) > 2) {
      colnames(XY)[1:3] <- c("ddlat", "ddlon", "tax")
      XY$tax <- as.character(XY$tax)
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