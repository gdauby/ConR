#' @title Identify Near Threatened Species
#'
#' @description Implementation of IUCN list of cases that qualify species as 
#' Near Threatened (NT) instead of Least Concerned (LC) or more generally
#' "LC or NT" (the standard output string in ConR).
#'
#' @param cats a character string containing the Red List categories.  
#' @param EOO numeric vector with species extent of occurrence - EOO (i.e.
#'   sub-criterion B1)
#' @param AOO numeric vector with species area of occupancy - AOO (i.e.
#'   sub-criterion B2)
#' @param decline string vector providing the status of the species continuing
#'   decline in EOO, AOO, habitat, locations or subpopulations or population
#'   size (i.e. condition 'b').
#' @param pop.reduction numeric vector with the estimated percentage of
#'   population reduction in the last three generations.
#' @param pop.size numeric vector with the mean estimate of population size in
#'   number of mature individuals.
#' @param pop.size.low numeric vector with the lower bound of the confidence
#'   interval of the population size in number of mature individuals.
#' @param locations numeric vector with the number of locations where the
#'   species occur (i.e. condition 'a')
#' @param sever.frag numeric vector with the proportion AOO which is in patches
#'   that are separated from other patches by a large distance.
#' @param subpop numeric vector with the number of sub-populationsfor the
#'   species
#' @param ext.fluct numeric vector with the mean order of magnitude of the
#'   differences between population minima and maxima.
#' @param subcriteria  character string with the sub-criteria used to perform
#'   the assessments (e.g. "A1")
#' @param many.more numeric value to numerically express what "many more" means.
#' @param extra.case logical. Should the extra case to detected probable "LC" be
#'   used? Default to TRUE.
#' 
#' @return A vector of the same length as ```cats``` with NT separated from the
#'   LC category. If the category provided is different than "LC" or "LC or NT"
#'   the function returns the same category.
#' 
#' @details 
#' 
#'  This function automatically identify the Near Threatened (NT) category among
#'  the "LC or NT" general category. A species qualify as NT if it is close to
#'  qualifying for the Vulnerable category (IUCN, 2019).
#'  
#'  The function try to translate the list of cases where NT applies and does
#'  not applies (IUCN 2019, pp 76-77). Not all the cases listed were translated
#'  here, particularly those based on uncertainties of the estimates. To perform
#'  this translation some interpretation, adaptation or generalization was
#'  carried out, namely:
#'  - "Population has declined by an estimated 20 - 25% in the last three generations: pop. reduction between 20 and 30%;
#'  - "Population has declined by an estimated 10%": pop. reduction >10%;
#'  - "many more": by default we use triple to numerically express many more;
#'  - "(...) has about 15,000/1,500 mature individuals": mature individuals between 10000 and 16000/1000 and 1600, respectively; 
#'  - "The population has more than 2,000 mature individuals": pop. size <2000 & 
#'  
#'  For many of the cases listed by IUCN (2019, pp 76 - 77), there is a minimum
#'  number of arguments that should be provided so that the NT criteria can be
#'  assigned. For instance, if only ```EOO``` and ```AOO``` is given no species
#'  will be classified as NT and the function will return "LC" for all "LC" or
#'  "LC or NT" provided as an input. On the other hand, if only ```pop.size```
#'  species can be classified as "NT". Please see all the IUCN cases for more
#'  details on this (IUCN 2019, pp 76-77).
#'  
#'  We added one extra case to the list of cases from IUCN (2019) where species
#'  should be classified as "LC" and not "NT". This case comes from an
#'  interpretation of the IUCN (2019) guidelines and it is unoficially used by
#'  the IUCN/SSC Global Tree Specialist Group as an indication of low
#'  probability of triggering a threatened or NT categories, if the species is
#'  not cited as threatened based on other criteria (e.g. uses and threats):
#'  - EOO > 30000 km2 or AOO > 3000 km2, and number of locations > 30. This 
#'  extra case can be excluded from the separation between "LC" and "NT" by 
#'  setting the argument ```extra.case``` to FALSE.
#'  
#' @author Renato A. Ferreira de Lima
#'
#' @references IUCN 2019. Guidelines for Using the IUCN Red List Categories and
#'   Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
#'   http://www.iucnredlist.org/documents/RedListGuidelines.pdf.
#'   
#' @examples
#' cats <- c("LC", "LC or NT", "LC or NT", "LC")
#' EOO <- c(35000, 10000, 10000, NA)
#' AOO <- c(3200, 800, 800, 32)
#' decline <- c("Decreasing", "Decreasing", "Increasing", "Increasing")
#' pop.reduction <- c(23.1, 23.1, 11, 6)
#' pop.size <- c(25000, 12000, 12000, 1200)
#' locations <- c(75, 25, 25, 9)
#' 
#' ## Example with different species metrics
#' near.threatened(cats = cats, EOO = EOO, AOO = AOO, decline = decline, 
#' pop.reduction = pop.reduction, pop.size = pop.size, locations = locations)
#'
#' ## Example with only EOO, AOO and number of locations (not enough metrics)
#' near.threatened(cats = cats, EOO = EOO, AOO = AOO)
#'
#' @export near.threatened
#' 
near.threatened <- function(cats = NULL,
                            EOO = NULL,
                            AOO = NULL,
                            decline = NULL,
                            pop.reduction = NULL,
                            pop.size = NULL,
                            pop.size.low = NULL,
                            locations = NULL,
                            sever.frag = NULL,
                            ext.fluct = NULL,
                            subpop = NULL,
                            subcriteria = NULL,
                            many.more = 3,
                            extra.case = TRUE) 
{
  
  output <- output1 <- data.frame(order = 1:length(cats),
                                  cats = cats, 
                                  stringsAsFactors = FALSE)
  output$lc.nt <- output1$lc.nt <- 1*cats %in% c("LC or NT","LC")
  
  ##1- Population has declined by an estimated 20 - 25% in the last three generations.
  if (!is.null(pop.reduction))
    output$case1 <- pop.reduction >= 20 & pop.reduction < 30
  
  ##2- The taxon meets the area requirements under criterion B for threatened 
  #(EOO <20,000 km2 and/or AOO < 2,000 km2) and is declining, but the population 
  #is not severely fragmented, occurs at many more than 10 locations, and there are no extreme fluctuations.
  if ((!is.null(EOO) | !is.null(AOO)) & !is.null(sever.frag) & !is.null(locations) & !is.null(ext.fluct))
    if (!is.null(EOO) & !is.null(AOO))
      output$case2 <- (EOO < 20000 | AOO < 2000) & sever.frag <= 0.5 & 
                          locations > 10*many.more & ext.fluct < 10
    if (!is.null(EOO) & is.null(AOO))
      output$case2 <- EOO < 20000 & sever.frag <= 0.5 & 
                          locations > 10*many.more & ext.fluct < 10
    if(is.null(EOO) & !is.null(AOO) )
      output$case2 <- AOO < 2000 & sever.frag <= 0.5 & 
                                    locations > 10*many.more & ext.fluct < 10

  ##3- The taxon meets the area requirements under criterion B for threatened (EOO <20,000 km2 and/or AOO <2,000 km2) and is severely fragmented, 
  #but the population is not declining, occurs at many more than 10 locations, and there are no extreme fluctuations.
    if((!is.null(EOO) | !is.null(AOO)) & !is.null(decline) & !is.null(sever.frag) & !is.null(locations) & !is.null(ext.fluct))
      if(!is.null(EOO) & !is.null(AOO))
        output$case3 <- (EOO < 20000 | AOO < 2000) & decline != "Decreasing" & 
                            sever.frag > 0.5 & locations > 10*many.more & ext.fluct < 10
      if(!is.null(EOO) & is.null(AOO))
        output$case3 <- EOO < 20000 & decline != "Decreasing" & 
                            sever.frag > 0.5 & locations > 10*many.more & ext.fluct < 10
      if(is.null(EOO) & !is.null(AOO))
        output$case3 <- AOO < 2000 & decline != "Decreasing" & 
                            sever.frag > 0.5 & locations > 10*many.more & ext.fluct < 10

  ##4- The taxon is declining and occurs at ten locations, but has an EOO of 30,000 km2 and/or an AOO of 3,000 km2, which are uncertain estimates.
    #not implemented
  ##5- The taxon is declining and severely fragmented, but has an EOO of 30,000 km2 and/or an AOO of 3,000 km2, which are uncertain estimates.
    #not implemented  
  ##6- The taxon is declining and severely fragmented, but has an EOO of 22,000 km2 and/or an AOO of 3,000 km2, which are highly certain estimates.
    #not implemented
  
  ##7- Population has declined by an estimated 10% in the last three generations, and is continuing to decline, and has about 15,000 mature individuals.
  if (!is.null(pop.reduction) & !is.null(decline) & !is.null(pop.size))
    output$case7 <- pop.reduction >= 9.9 & 
                    decline == "Decreasing" & pop.size >= 10000 & pop.size < 16000
  
  ##8- The taxon exists in a single subpopulation of about 15,000 individuals and is declining.
  if (!is.null(subpop) & !is.null(pop.size) & !is.null(decline))
    output$case8 <- subpop == 1 & pop.size >= 10000 & pop.size < 16000 & decline == "Decreasing"
  
  ##9- The population has about 1,500 mature individuals.
  if (!is.null(pop.size))
    output$case9 <- pop.size >= 1000 & pop.size < 2000
  
  ##10- The best estimate of population size is 2,000 mature individuals, but this estimate is very uncertain, and as low as 1,000 mature individuals cannot be ruled out.
  if (!is.null(pop.size) & !is.null(pop.size.low))
    output$case10 <- pop.size >= 2000 & pop.size.low < 1600
  
  ##11- The taxon exists at three sites, occupying an area of 12 km2; the population is being harvested but is not declining; there are no current threats, 
  #but there are plausible events that may cause the species to decline, but these are unlikely to make the species Extinct or Critically Endangered in a short time.
  if (!is.null(locations) & !is.null(AOO) & !is.null(decline))
    output$case11 <- locations == 3 & AOO == 12 & decline != "Decreasing"
  
  ##12- Population has declined by 40% in the last three generations, but the decline has stopped, and the causes of the decline have been understood.
  if (!is.null(pop.reduction) & !is.null(decline) & "A2" %in% subcriteria)
    output$case11 <- pop.reduction >= 40 & decline != "Decreasing" & "A2" %in% subcriteria
  
  ## The following are examples of species that should not be listed as NT (or any of the categories of threat), 
  #unless other criteria apply:
    
  ##1- Population has declined by an estimated 10% in the last three generations, and there are more than 20,000 mature individuals.
  if (!is.null(pop.reduction) & !is.null(pop.size))
    output1$case1 <- pop.reduction < 10.5 & pop.size > 20000
  
  ##2- Population has declined by an estimated 30% as part of fluctuations.
  if (!is.null(pop.reduction) & !is.null(ext.fluct))
    output1$case2 <- pop.reduction < 30 & ext.fluct > 10
  
    
  ##3- The taxon meets the area requirements under criterion B for CR (EOO <100 km2 and/or AOO <10 km2), but is not declining, not severely fragmented, 
  #there are no extreme fluctuations, and there are no obvious threats.
  if ((!is.null(EOO) | !is.null(AOO)) & !is.null(decline) & !is.null(ext.fluct))
    if (!is.null(EOO) & !is.null(AOO))
      output1$case3 <- (EOO < 100 | AOO < 10) & decline != "Decreasing" & ext.fluct < 10
  if (!is.null(EOO) & is.null(AOO))
    output1$case3 <- EOO < 100 & decline != "Decreasing" & ext.fluct < 10
  if (is.null(EOO) & !is.null(AOO))
    output1$case3 <- AOO < 10 & decline != "Decreasing" & ext.fluct < 10
  
  ##4- The taxon is long-lived and slow growing, but does not meet any criteria A-E.
    #Not implemented
      
  ##5- The population has more than 2,000 mature individuals.
  # if (!is.null(pop.size) & !is.null(subpop))
  #   output1$case5 <- pop.size > 2000 & subpop == 1
  
  ##6- The taxon exists at three sites, occupying an area of 30 km2; the population is not declining; there are no current threats, 
  #and the species is very unlikely to become Extinct or Critically Endangered in a short time.   
  
  if (!is.null(locations) & !is.null(AOO) & !is.null(decline))
    output1$case6 <- locations >= 3 & AOO >= 30 & decline != "Decreasing"
  
  ##7 - Additional case from GTA: EOO >30,000 km2; AOO >3,000 km2; and locations >30 (unique spatial occurrences 10 km appart)
  if (extra.case) {
    if ((!is.null(EOO) | !is.null(AOO)) & !is.null(locations))
      if (!is.null(EOO) & !is.null(AOO))
        output1$case7 <- (EOO > 30000 | AOO > 3000) & locations > 30
    if (!is.null(EOO) & is.null(AOO))
      output1$case7 <- EOO > 30000 & locations > 30
    if (is.null(EOO) & !is.null(AOO))
      output1$case7 <- AOO > 3000 & locations > 30
  }
  
  #### Merging the results #### 
  if (dim(output)[2] > 3) {
    output$nt <- 1*apply(output[,-c(1:3)], 1, any)
    output <- output[,c("order","cats","lc.nt","nt")]
  } else { 
    output$nt <- 0
    output <- output[,c("order","cats","lc.nt","nt")]
    
  }
  
  if (dim(output1)[2] > 3) {
    output1$non.nt <- (-1)*apply(output1[,-c(1:3)], 1, any)
    output1 <- output1[,c("order","non.nt")]
  } else { 
    output1$non.nt <- 0
    output1 <- output1[,c("order","non.nt")]
  }  
  
  output2 <- merge(output, output1, by = "order", sort = FALSE)
  output2$final <- apply(output2[, c("lc.nt", "nt", "non.nt")], 1, sum, na.rm = TRUE)

  # Separating NT and LC apart
  final.cats <- cats
  final.cats[output2$lc.nt %in% 1 & output2$final %in% 2] <- "NT"
  final.cats[output2$lc.nt %in% 1 & output2$final %in% c(0, 1)] <- "LC"
  
  return(final.cats)
}
