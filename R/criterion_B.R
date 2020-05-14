#' @title Assess IUCN Criterion B
#'
#' @description Preliminary assessment of species conservation status following
#'  IUCN Criterion B, which is based on species geographic distribution (extent
#'  of occurrence)
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
criterion_B = function(DATA, protec.areas = NULL, 
                       NamesSp = "species1", 
                       file_name = NULL, 
                       #add.legend = FALSE, DrawMap = FALSE, map_pdf = FALSE, draw.poly.EOO = FALSE, 
                       EOO.threshold = c(20000, 5000, 100), AOO.threshold = c(2000, 500, 10), Loc.threshold = c(10, 5, 1)) {
  if (is.null(protec.areas)) {
    tmp <- as.data.frame(matrix(NA, 4, dim(DATA)[2]))
    rownames(tmp) <- c("Category_CriteriaB", "Category_code","Category_AOO", "Category_EOO")
    colnames(tmp) <- colnames(DATA)
    Results = rbind.data.frame(DATA, tmp, make.row.names = TRUE, deparse.level = 2)
  } else {
    tmp <- as.data.frame(matrix(NA, 5, dim(DATA)[2]))
    rownames(tmp) <- c("Category_CriteriaB", "Category_code","Ratio_occ_within_PA", "Category_AOO", "Category_EOO")
    colnames(tmp) <- colnames(DATA)
    Results = rbind.data.frame(DATA, tmp, make.row.names = TRUE, deparse.level = 2)
  }
  
  if (any(Results["EOO",] < Results["AOO",])) 
    Results["EOO",][Results["EOO",] < Results["AOO",]] <- Results["AOO",][Results["EOO",] < Results["AOO",]]
  
  #### CHECK HERE ####
  if (!is.null(protec.areas)) 
    Results["Ratio_occ_within_PA", 1] <- round(length(which(!is.na(Links_NatParks[, 1])))/nrow(Links_NatParks) * 100, 1)
  #### END CHECK ####
  
  Rank_EOO = rep(4, dim(DATA)[2])
  Rank_EOO[Results["EOO",] < EOO.threshold[1]] <-3
  Rank_EOO[Results["EOO",] < EOO.threshold[2]] <-2
  Rank_EOO[Results["EOO",] < EOO.threshold[3]] <-1
  
  Rank_AOO = rep(4, dim(DATA)[2])
  Rank_AOO[Results["AOO",] < AOO.threshold[1]] <-3
  Rank_AOO[Results["AOO",] < AOO.threshold[2]] <-2
  Rank_AOO[Results["AOO",] < AOO.threshold[3]] <-1
  
  Rank_Loc = rep(4, dim(DATA)[2])
  Rank_Loc[Results["Nbe_loc",] < Loc.threshold[1]] <-3
  Rank_Loc[Results["Nbe_loc",] < Loc.threshold[2]] <-2
  Rank_Loc[Results["Nbe_loc",] < Loc.threshold[3]] <-1
  
  Rank_B1a <- apply(rbind(Rank_EOO, Rank_Loc), 2, max, na.rm=TRUE)
  Rank_B2a <- apply(rbind(Rank_AOO, Rank_Loc), 2, max, na.rm=TRUE)
  Rank_CriteriaB <- apply(rbind(Rank_B1a, Rank_B2a), 2, min, na.rm=TRUE)
  
  Nbe_Loc <- as.numeric(Results["Nbe_loc", ])
  Cat = rep("LC or NT", dim(DATA)[2])
  Cat[Rank_CriteriaB == 1] <- "CR"
  Cat[Rank_CriteriaB == 2] <- "EN"
  Cat[Rank_CriteriaB == 3 & Nbe_Loc > 0 & Nbe_Loc < 11] <- "VU"
  
  Cat_Code <- rep(NA, dim(DATA)[2])
  if(any(Rank_B1a > Rank_B2a))    
    Cat_Code[Rank_B1a > Rank_B2a] <- paste(Cat[Rank_B1a > Rank_B2a], "B2a")
  if(any(Rank_B1a < Rank_B2a)) 
    Cat_Code[Rank_B1a < Rank_B2a] <- paste(Cat[Rank_B1a < Rank_B2a], "B1a")
  if (any(Rank_B1a == Rank_B2a)) 
    Cat_Code[Rank_B1a == Rank_B2a & Rank_B1a != 4] <- paste(Cat[Rank_B1a == Rank_B2a & Rank_B1a != 4], "B1a+B2a") 
  
  if (!is.null(protec.areas)) {
    if (any(as.numeric(Results["Ratio_occ_within_PA", ]) == 100)) {
      Results["Category_CriteriaB", as.numeric(Results["Ratio_occ_within_PA", ]) == 100] <- "LC or NT"
      Results["Category_code", as.numeric(Results["Ratio_occ_within_PA", ]) == 100] <- Cat_Code
    } else {
      Results["Category_CriteriaB", ] <- Cat
      Results["Category_code", ] <- Cat_Code
    }
  } else {
    Results["Category_CriteriaB", ] <- Cat
    Results["Category_code", ] <- Cat_Code
  }
  
  if (any(Rank_B2a == 1)) 
    Results["Category_AOO", Rank_B2a == 1] <- "CR"
  if (any(Rank_B2a == 2)) 
    Results["Category_AOO", Rank_B2a == 2] <- "EN"
  if (any(Rank_B2a == 3)) 
    Results["Category_AOO", Rank_B2a == 3] <- "VU"
  if (any(Rank_B2a > 3)) 
    Results["Category_AOO", Rank_B2a > 3] <- "LC or NT"
  if (any(Rank_B1a == 1)) 
    Results["Category_EOO", Rank_B1a == 1] <- "CR"
  if (any(Rank_B1a == 2)) 
    Results["Category_EOO", Rank_B1a == 2] <- "EN"
  if (any(Rank_B1a == 3)) 
    Results["Category_EOO", Rank_B1a == 3] <- "VU"
  if (any(Rank_B1a > 3)) 
    Results["Category_EOO", Rank_B1a > 3] <- "LC or NT"

  #### CHECK HERE ####
  if (!is.null(protec.areas)) 
    Results["Ratio_occ_within_PA", 1] <- round(length(which(!is.na(Links_NatParks[, 1])))/nrow(Links_NatParks) * 100, 2)
  #### END CHECK ####
  
  if (!is.null(protec.areas)) {
    if (any(as.numeric(Results["Ratio_occ_within_PA", ]) == 100)) {
      Results["Category_CriteriaB", as.numeric(Results["Ratio_occ_within_PA", ]) == 100] <- "LC or NT"
    } else {
      if (any(as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3])) {
        Results["Category_AOO", as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc",
                                                                                                     ]) <= Loc.threshold[3]] <- Results["Category_CriteriaB", as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3]] <- "CR"
      } else {
        if (any(as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2])) {
          Results["Category_AOO", as.numeric(Results["AOO", 1]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", 
                                                                                                        ]) <= Loc.threshold[2]] <- Results["Category_CriteriaB", as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2]] <- "EN"
        } else {
          if (any(as.numeric(Results["AOO", ]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1])) {
            Results["Category_AOO", as.numeric(Results["AOO", 1]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", 
                                                                                                          ]) <= Loc.threshold[1]] <- Results["Category_CriteriaB", as.numeric(Results["AOO", ]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
          } else {
            Results["Category_AOO", ] <- Results["Category_CriteriaB", ] <- "LC or NT"
          }
        }
      }
    }
  } else {
    if (any(as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3])) {
      Results["Category_AOO", as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3]] <- "CR"
      Results["Category_CriteriaB", as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[3]] <- "CR"
    } else {
      if (any(as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2])) {
        Results["Category_AOO", as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2]] <- "EN"
        Results["Category_CriteriaB", as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[2]] <- "EN"
      } else {
        if (any(as.numeric(Results["AOO", ]) < AOO.threshold[1] &  as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1])) {
          Results["Category_AOO", as.numeric(Results["AOO", ]) < AOO.threshold[1] &  as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
          Results["Category_CriteriaB", as.numeric(Results["AOO", ]) < AOO.threshold[1] &  as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
        } else {
          #Results["Category_AOO", ] <- "LC or NT"
          #Results["Category_CriteriaB", ] <- "LC or NT"
        }
      }
    }
  }
  if (any(is.na(Results["Category_AOO", ]))) {
    if (any(as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[3])) {
      Results["Category_AOO", as.numeric(Results["AOO", ]) < AOO.threshold[3] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[3]] <- "CR"
    } else {
      if (any(as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[2])) {
        Results["Category_AOO", as.numeric(Results["AOO", ]) < AOO.threshold[2] & as.numeric(Results["Nbe_loc",  ]) <= Loc.threshold[2]] <- "EN"
      } else {
        if (any(as.numeric(Results["AOO", 1]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1])) {
          Results["Category_AOO", as.numeric(Results["AOO", 1]) < AOO.threshold[1] & as.numeric(Results["Nbe_loc", ]) <= Loc.threshold[1]] <- "VU"
        } else {
          #Results["Category_AOO", ] <- "LC or NT"
        }
      }
    }
  }
  #Results["Category_code", ] <- paste(Results["Category_CriteriaB", ], "B2a")
  
  Results1 = as.data.frame(t(Results), stringsAsFactors = FALSE)  
  return(Results1)
}
