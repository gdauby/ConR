#' @title Internal function
#'
#' @description Alpha hull processing
#'
#' @param x ahull class object
#'
#' @details 
#' The functions ahull_to_SPLDF and alpha.hull.poly were originally posted in the website https://casoilresource.lawr.ucdavis.edu/software/r-advanced-statistical-package/working-spatial-data/converting-alpha-shapes-sp-objects/
#' in a now broken link. It is also used in functions written by David Bucklin, see https://github.com/dnbucklin/r_movement_homerange 
#'
#'
#' 
ahull_to_SPLDF <- function(x)
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
    sp::SpatialLines(list(l), proj4string = sp::CRS(as.character(NA), doCheckCRSArgs=TRUE))
  
  # promote to SpatialLinesDataFrame, required for export to GRASS / OGR
  l.spldf <-
    sp::SpatialLinesDataFrame(l.spl, data = data.frame(id = 1), match.ID =
                                FALSE)
  
  return(l.spldf)
}