#' @title Internal function
#'
#' @description Activate paralle processing
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#' 
#' @param XY data.frame
#' @param mode character string either 'spheroid' or 'planar'. By default 'spheroid'
#' @param proj_type crs
#' 
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel makePSOCKcluster
#' @importFrom foreach %dopar% %do% foreach
#' @keywords internal
#' @export
activate_parallel <- function(parallel = FALSE) {
  if (parallel) {
    cl <- parallel::makePSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
}