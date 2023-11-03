
#' @title Internal function
#'
#' @description Display progress
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#' 
#' @param show_progress logical
#' @param max_pb integer
#' 
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
display_progress_bar <- function(show_progress = TRUE, max_pb) {
  if (show_progress) {
    pb <-
      txtProgressBar(min = 0,
                     max = max_pb,
                     style = 3)
    
    progress <- function(n)
      setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
  } else {
    opts <- NULL
    pb <- NA
  }
  return(list(opts = opts, pb = pb))
}