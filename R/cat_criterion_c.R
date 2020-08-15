#' @title Categorize taxa according to IUCN criterion C
#'
#' @description Provide the consensus IUCN category based on the sub-criteria of
#'   IUCN criterion C (C1, C2ai, C2aii and C2b) and the thresholds recommended
#'   by IUCN.
#'
#' @param C1_df data frame with the parameters necessary to assess IUCN
#'   sub-criterion C1.
#' @param C2_df data frame with the parameters necessary to assess IUCN
#'   sub-criteria C2.
#' @param C.threshold numeric vector with the criterion C thresholds to define
#'   small population sizes (e.g. number of mature individuals). Default values
#'   are the thresholds recommended by the IUCN.    
#' @param C1.threshold numeric vector with the C1 thresholds to convert
#'   continuing decline estimates into categories. Default is the thresholds
#'   recommended by IUCN.
#' @param C2ai.threshold numeric vector with the C2ai thresholds to assess the
#'   number of mature individuals in each subpopulation. Default is the
#'   thresholds recommended by IUCN.
#' @param C2aii.threshold numeric vector with the C2aii thresholds to assess the
#'   proportion of mature individuals in one subpopulation. Default is the
#'   thresholds recommended by IUCN.
#' @param mag.fluct numerical. Threshold of mean order of magnitude of the
#'   differences between population minima and maxima to classify populations
#'   with extreme fluctuations. Default to 10 as recommended by IUCN (2019).
#' @param high.alter numerical. Threshold of proportion of changes that are
#'   followed by a change in the opposite direction. Default to 80%, but
#'   currently not implemented
#' @param all.cats logical. Should the categories from all criteria be returned
#'   and not just the consensus categories? Default to TRUE.
#'
#' @return A list containing a vector of the consensus category from all
#'   sub-criteria evaluated for each taxon (`ranks_C`) and the sub-criteria used
#'   to obtain the consensus category (`cat_codes`). If `all.cats  == TRUE` the
#'   function also returns a data frame containing the categories classified by
#'   each sub-criteria individually (`all.cats`).
#'
#' @details By default, the function provides the consensus category, following
#'   the recommendations of IUCN (2019) that states "Only the criteria for the
#'   highest category of threat that the taxon qualifies for should be listed".
#'   Therefore, the consensus category is the highest category of threat among
#'   the sub-criteria evaluated.
#'
#'   The function assumes that the order of the values in C1_df and C2_df are
#'   from the same taxa (i.e. first element from C1_df and C2_df is always the
#'   same species i). Therefore, the order of the estimates of population size
#'   and continuing decline for each sub-criterion *must* be the same.
#'
#'
#' @author Lima, R.A.F. & Dauby, G.
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
# @importFrom stringr str_replace_all
#'
#' @importFrom stringr str_replace_all
#'
#' @export cat_criterion_c
#'
#'
cat_criterion_c <- function(C1_df = NULL,
                            C2_df = NULL,
                            C.threshold = c(10000, 2500, 250),
                            C1.threshold = c(10, 20, 25),
                            C2ai.threshold = c(1000, 250, 50), 
                            C2aii.threshold = c(90, 95, 100),
                            mag.fluct = 10,
                            high.alter = 80,
                            all.cats = TRUE) {
  
  L <- list(C1 = C1_df,
            C2 = C2_df)
  rank <- L[!unlist(lapply(L, is.null))]
  rank <- sapply(rank, function(x) rep(NA, dim(x)[1]))
  cat <- code <- rank
  
  if((is.data.frame(C1_df) & is.data.frame(C2_df)) & length(unique(lapply(L, function(x) dim(x)[1]))) > 1)
    stop("Numbers of data frame rows provided for each criterion should be identical")
  
  if(any(unlist(lapply(L, is.null))))
    warning(paste("The following subcriteria were not used in the assessment: ", paste(names(L)[unlist(lapply(L, is.null))], collapse = ", ")),
            call. = FALSE)
  
  #### CHECK HERE: Should we return this warning? ####
  # if(any(unlist(lapply(L, function(x) any(apply(x[,-c(1,2)], 1, function(x) any(x<0)))))))
  #   warning(paste("The following subcriteria had population increase for one or more species: ", paste(names(L)[unlist(lapply(L, function(x) any(x<0)))], collapse = ", ")),
  #           call. = FALSE)
  
  rpl.cds <- c("CR", "EN", "VU", "LC or NT")
  names(rpl.cds) <- c("3", "2", "1", "0")
  
  ## Criteria C1: under criterion C1, the decline must be observed or estimated (thus removing projections of future decline)
  if(!is.null(C1_df)) { 
    
    tmp <- L[["C1"]]
    tmp$check <- tmp$assess.pop.size < max(C.threshold) &
                    grepl("Decreasing", tmp$cont.decline)
    
    tmp$red_3gen <- sapply(1:length(tmp$reduction_3gen), 
                                 function(i) if (tmp$assess.pop.size[i] < C.threshold[1] & tmp$reduction_3gen[i] >= C1.threshold[1]) 1 else 0)
    tmp$red_2gen <- sapply(1:length(tmp$reduction_2gen), 
                                 function(i) if (tmp$assess.pop.size[i] < C.threshold[2] & tmp$reduction_2gen[i] >= C1.threshold[2]) 2 else 0)
    tmp$red_1gen <- sapply(1:length(tmp$reduction_1gen), 
                                 function(i) if (tmp$assess.pop.size[i] < C.threshold[3] & tmp$reduction_1gen[i] >= C1.threshold[3]) 3 else 0)

    rank[,"C1"] <- apply(tmp[,c("red_3gen", "red_2gen", "red_1gen")], 1, max)
    rank[,"C1"][!tmp$check] <- 0
    cat[,"C1"] <- stringr::str_replace_all(rank[,"C1"], rpl.cds)
    code[,"C1"][rank[,"C1"]>0] <- "C1"

  }
  
  ## Criteria C2
  if(!is.null(C2_df)) {  
    
    tmp <- L[["C2"]]
    tmp$check <- tmp$assess.pop.size < max(C.threshold) &
                    grepl("Decreasing", tmp$any.decline)
    
    tmp$pop <- 3 - findInterval(tmp$assess.pop.size, sort(C.threshold))
    tmp$each <- 3 - findInterval(tmp$max.subpop.size, sort(C2ai.threshold))
    tmp$prop <- 0
    tmp$prop[tmp$pop == 1 & tmp$prop.subpop.size >= C2aii.threshold[3]] <- 2
    tmp$prop[tmp$pop == 2 & tmp$prop.subpop.size >= C2aii.threshold[2]] <- 2
    tmp$prop[tmp$pop == 3 & tmp$prop.subpop.size >= C2aii.threshold[1]] <- 3
    
    tmp$C2ai = rep(0, dim(tmp)[1])
    tmp$C2ai[tmp$check & tmp$pop >= 1 & tmp$each >= 1] <- 1
    tmp$C2ai[tmp$check & tmp$pop >= 2 & tmp$each >= 2] <- 2
    tmp$C2ai[tmp$check & tmp$pop >= 3 & tmp$each >= 3] <- 3
    
    tmp$C2aii = rep(0, dim(tmp)[1])
    tmp$C2aii[tmp$check & tmp$pop == 1 & tmp$prop == 1] <- 1
    tmp$C2aii[tmp$check & tmp$pop == 2 & tmp$prop == 2] <- 2
    tmp$C2aii[tmp$check & tmp$pop == 3 & tmp$prop == 3] <- 3
    
    #### CHECK HERE: INCLUDE ANY CRITERIA RELATED TO HIGH ALTERNANCE IN THE FLUCTUATIONS AS WELL? ####
    tmp$C2b <- tmp$pop
    tmp$C2b[tmp$mean.fluctuation < mag.fluct] <- 0
    tmp$C2b[!tmp$check] <- 0
    
    rank[,"C2"] <- apply(tmp[,c("C2ai", "C2aii", "C2b")], 1, max)
    cat[,"C2"] <- stringr::str_replace_all(rank[,"C2"], rpl.cds)
    code[,"C2"] <- suppressWarnings(apply(tmp[,c("C2ai", "C2aii", "C2b")], 1,
                                          FUN = function(x) {
                                            y = names(x[x == max(x[x > 0], na.rm = T)])
                                            paste(y[!is.na(y)], collapse = "+")
                                          }))
    
  }    
  
  ranks_C <- as.character(apply(rank, 1, max, na.rm=TRUE))
  ranks_C <- stringr::str_replace_all(ranks_C,  rpl.cds)
  
  code[code == ""] <- NA
  cats_code <- sapply(1:dim(code)[1], 
                      function(i) {
                        y <- code[i,][which(rank[i,] == max(rank[i,], na.rm = TRUE))]
                        paste(y[!is.na(y)], collapse = "+")
                      })
  
  if(all.cats & dim(cat)[1] > 1) {
    
    return(list(ranks_C = ranks_C, cats_code = cats_code, all_cats = as.data.frame(cat)))
    
  } else {
    
    return(list(ranks_C = ranks_C, cats_code = cats_code))
    
  }
}