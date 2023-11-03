#' Dataset of plant species distribution
#'
#' A dataframe of three columns
#'
#'
#' @keywords datasets
#' @name dataset.ex
#' @usage data(dataset.ex)
#' @format A dataframe
NULL

#' Dataset to illustrate the assessment of IUCN criterion A
#'
#' A dataset containing the population sizes of six taxa and from 1970 until
#' 2032 at every two years. First two species are taken from the same example
#' used by IUCN (2019) to illustrate the assessment of this criterion.
#'
#' @keywords datasets
#' @name example_criterionA
#' @usage data(example_criterionA)
#' @format A data frame with six rows and 33 variables
#' @source <https://www.iucnredlist.org/resources/criterion-a>
NULL

#' Dataset to illustrate the assessment of IUCN criterion C
#'
#' A dataset containing the population sizes of three subpopulations for nine
#' different taxa and from 1970 until 2000 at every five years. 
#'
#' @keywords datasets
#' @name example_criterionC_subpops
#' @usage data(example_criterionC_subpops)
#' @format A data frame with 27 rows and 8 variables
NULL

#' Dataset to illustrate the assessment of IUCN criterion C
#'
#' A dataset containing the population sizes for nine
#' different taxa and from 1970 until 2000 at every five years,
#' without information on subpopulation sizes. 
#'
#' @keywords datasets
#' @name example_criterionC
#' @usage data(example_criterionC)
#' @format A data frame with 9 rows and 7 variables
NULL

#' Dataset to illustrate the assessment of populations with extreme fluctuation
#'
#' A dataset containing the population sizes adapted from some of the panels 
#' provided in Figure 4.4. of IUCN (2019) to illustrate different examples
#' of how population can fluctuate.
#'
#' @keywords datasets
#' @name example_fluctuation
#' @usage data(example_fluctuation)
#' @format A data frame with seven rows and 28 variables
#' @source <http://www.iucnredlist.org/documents/RedListGuidelines.pdf>
NULL


#' Dataset to run the package introductory vignette
#'
#' A dataset containing occurrence records and population size estimates from
#' three arborescent species from the Atlantic Forest obtained from different
#' sources, but mostly from GBIF and TreeCo.
#' 
#' Only occurrence records considered to have a high confidence level in terms
#' of georefferencing and species identifications are included.
#' 
#' In data frame 'occurrences' columns are: ddlat, ddlon, tax, coly, UC, source
#' 
#' In data frame 'population.sizes' columns are: species, and all years with 
#' population estimates
#'
#' @keywords datasets
#' @name example_tutorial
#' @usage data(example_tutorial)
#' @format A named list with two data frames: 'occurrences' and 'population.sizes'
#' @source <https://doi.org/10.15468/dl.mzmat2>
NULL
