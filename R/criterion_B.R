#' @title Assess IUCN Criterion B
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion B, which is based on species geographic distribution (i.e. extent
#'  of occurrence - EOO, and area of occupancy, AOO)
#'
#' @param x the character string.
#'
#' @return 
#' 
#' @details The function ... 
#' 
#' @author Dauby, G. & Lima, R.A.F.
#'
#' @references 
#'
#' @export criterion_B
#'
#' @examples
#'
#'
criterion_B = function(x, 
                       protec.areas = NULL, 
                       NamesSp = "species1", 
                       file_name = NULL, 
                       #add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                       EOO.threshold = c(20000, 5000, 100), 
                       AOO.threshold = c(2000, 500, 10), 
                       Loc.threshold = c(10, 5, 1)) {
  
  if (is.null(protec.areas)) {
    
    tmp <- as.data.frame(matrix(NA, 4, dim(DATA)[2]))
    rownames(tmp) <- c("Category_CriteriaB", "Category_code","Category_AOO", "Category_EOO")
    colnames(tmp) <- colnames(DATA)
    Results <- rbind.data.frame(DATA, tmp, make.row.names = TRUE, deparse.level = 2)
    
  } else {
    
    tmp <- as.data.frame(matrix(NA, 5, dim(DATA)[2]))
    rownames(tmp) <- c("Category_CriteriaB", "Category_code","Ratio_occ_within_PA", "Category_AOO", "Category_EOO")
    colnames(tmp) <- colnames(DATA)
    Results = rbind.data.frame(DATA, tmp, make.row.names = TRUE, deparse.level = 2)
    
  }
  
  # "If EOO is less than AOO, EOO should be changed to make it equal to AOO to ensure
  # consistency with the definition of AOO as an area within EOO."
  if (any(Results["EOO",] < Results["AOO",])) 
    Results["EOO",][Results["EOO",] < Results["AOO",]] <- 
      Results["AOO",][Results["EOO",] < Results["AOO",]]
  
  #### CHECK HERE ####
  if (!is.null(protec.areas)) 
    Results["Ratio_occ_within_PA", 1] <- round(length(which(!is.na(Links_NatParks[, 1])))/nrow(Links_NatParks) * 100, 1)
  #### END CHECK ####
  
  
  cat_criterion_b(
    EOO = ,
    AOO = ,
    locations = ,
    protected = ,
    EOO.threshold = ,
    AOO.threshold = ,
    Loc.threshold =
  )
  
  #Results["Category_code", ] <- paste(Results["Category_CriteriaB", ], "B2a")
  
  Results1 = as.data.frame(t(Results), stringsAsFactors = FALSE)  
  return(Results1)
}






