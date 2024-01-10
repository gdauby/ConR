#' @title Categorize Taxa Using Multiple IUCN Criteria
#' 
#' @description Provide the consensus IUCN category based on multiples IUCN
#'   sub-criteria.
#'
#' @param assess.df a data frame containing the taxon name in the first columns
#'   and the assessments under each criterion in the subsequent columns.
#' @param evidence.df a data frame (not currently implemented)
#' 
#' @return The same data frame as ```assess.df``` with three new columns: the
#'   consensus category ('category'), the main criteria that lead to this
#'   category ('main.criteria'). It also returns the auxiliary category provided
#'   by other criteria ('aux.criteria'), separated by a ';'
#' 
#' @details The definiton of the main category of threat, follows the
#'   recommendations of IUCN (2019) that states "Only the criteria for the
#'   highest category of threat that the taxon qualifies for should be listed".
#'   Therefore, the consensus category is the highest category of threat among
#'   the sub-criteria evaluated. Nevertheless, the function also returns the
#'   categories and sub-criteria related to lower categories of threat.
#' 
#' @author Renato A. Ferreira de Lima
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'
#' @examples
#' df <- data.frame(tax = c("sp1","sp2","sp3","sp4"),
#'                     A2 = c(NA, "VU", "VU", NA),
#'                     B1 = c("LC", "VU", "LC", "LC"),
#'                     B2 = c("LC", "EN", "LC", "VU"),
#'                     D = c(NA, "LC", "LC", NA))
#' cat_mult_criteria(df)                     
#' 
#'
#'
#' @importFrom stringr str_replace_all
#' 
#' @export cat_mult_criteria
cat_mult_criteria <- function(assess.df = NULL, evidence.df = NULL){
  
  assess.df$tmp.order.plantr <- 1:dim(assess.df)[1]
  tmp <- assess.df
  tax <- names(tmp)[1]
  
  possible_crits <- c(LETTERS[1:5], paste0(LETTERS[1], 1:4),
                     paste0(rep(LETTERS[2:4], each = 2), 1:2),
                     paste0("category_", LETTERS[1:5]), "D2.Loc", "D.AOO")
  if (!any(names(tmp) %in% possible_crits))
    stop("Please provide a data frame with at least one of the following column names: ", 
         paste0(possible_crits, collapse = ", "))
  col.name.df <- names(tmp)[names(tmp) %in% possible_crits]
  tmp0 <- tmp[, names(tmp) %in% possible_crits]
  n.crits <- dim(tmp0)[2]
  
  tmp <- cbind.data.frame(tax = tmp[[1]], tmp0)
  
  ## Establishing an hierarchy of data availability
  hier <- c("Observed","Estimated","Projected","Inferred","Suspected")
  names(hier) <- c("2", "2", "1", "1", "0")
  
  ## Replacing the categories by ordered numbers
  tmp[] <- lapply(tmp, gsub, pattern = "^LC or NT$", replacement = 0, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^LC$", replacement = 0, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^DD$", replacement = 0.5, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^NT$", replacement = 1, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^VU$", replacement = 2, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^EN$", replacement = 3, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^CR$", replacement = 4, perl = TRUE)
  tmp[] <- lapply(tmp, gsub, pattern = "^EW$|^EX$|^RE$", replacement = 5, perl = TRUE)
  
  ## Getting the main (highest) criteria for each species
  tmp[,2:n.crits] <- apply(tmp[,2:n.crits], 2, as.double)
  tmp$category <- suppressWarnings(as.character(apply(tmp[,2:n.crits], 1, max, na.rm=TRUE)))
  tmp$category[tmp$category %in% c("Inf", "-Inf")] <- ""
  
  tmp$main.criteria <- 
    apply(tmp[,2:n.crits], 1, 
          function(x) paste(names(x)[which(x == suppressWarnings(max(x, na.rm = TRUE)))], 
                            collapse = "+"))
  rpl.cds <- c("EW","EX","RE","CR", "EN", "VU", "NT", "DD","LC", "LC or NT")
  names(rpl.cds) <- c("5", "5", "5", "4", "3", "2", "1", "0.5", "0", "0")
  
  ## Getting the auxiliary criteria for each species
  tmp$aux.criteria <- 
    apply(tmp[,2:n.crits], 1, function(x) 
      if(length(unique(x)) > 1) {
        crits <- sort(unique(as.double(x)), decreasing = TRUE)[-1]
        cats <- sapply(crits, function(y) paste(names(tmp[,2:n.crits])[x==y], collapse = "+"))
        paste(paste(stringr::str_replace_all(crits, rpl.cds),": ", cats, sep=""), collapse = "; ")
      } else {
        ""
      }  
    )
  
  #Final edits on the results
  tmp$aux.criteria <- gsub("^: $", "", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("\\+NA\\+", "+", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("NA\\+", "+", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("\\+NA", "+", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("\\+\\+", "+", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub(": \\+", ": ", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("\\+$", "", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("\\+;", ";", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("category_B", "B1+B2", tmp$aux.criteria, perl = TRUE)
  tmp$aux.criteria <- gsub("category_C", "C1+C2", tmp$aux.criteria, perl = TRUE)
  
  tmp$main.criteria[is.infinite(as.double(tmp$category))] <- ""
  tmp$category[is.infinite(as.double(tmp$category))] <- ""
  tmp$category <- stringr::str_replace_all(tmp$category,  rpl.cds)
  tmp$category <- gsub("category_B", "B1+B2", tmp$category, perl = TRUE)
  tmp$category <- gsub("category_C", "C1+C2", tmp$category, perl = TRUE)
  
  #Merging with the entry data.frame
  res <- merge(assess.df, tmp[,c(tax,"category","main.criteria","aux.criteria")], 
               by = tax)
  res <- res[order(res$tmp.order.plantr),]
  res <- res[,-which(names(res) %in% "tmp.order.plantr")]
  
  return(res)
}
