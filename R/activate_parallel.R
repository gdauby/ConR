#' @title Internal function
#'
#' @description Activate paralle processing
#'
#' @author Gilles Dauby, \email{gildauby@gmail.com}
#' 
#' @param parallel a logical. Whether running in parallel. By default, it is FALSE
#' @param NbeCores an integer. Register the number of cores for parallel execution. By default, it is 2
#' 
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel makePSOCKcluster
#' @importFrom foreach %dopar% %do% foreach
#' @keywords internal
#' @export
activate_parallel <- function(parallel = FALSE, NbeCores = 2) {
  if (parallel) {
    cl <- parallel::makePSOCKcluster(NbeCores)
    doSNOW::registerDoSNOW(cl)
    
    message('Parallel running with ',
            NbeCores, ' cores')
    
  } else {
    cl <- NULL
  }
  return(cl)
}


c_par <- function(parallel = FALSE) {
  
  if (parallel) {
    `%d%` <- foreach::`%dopar%`
  } else {
    `%d%` <- foreach::`%do%`
  }
  return(`%d%`)
}


