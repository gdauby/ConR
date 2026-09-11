
<!-- README.md is generated from README.Rmd. Please edit README.Rmd and then
     regenerate README.md with devtools::build_readme(). Do not edit README.md
     by hand: your changes will be overwritten. -->

<img src="https://raw.githubusercontent.com/gdauby/ConR/devel/inst/figures/conr_sticker4.png" align="right" alt="" width="120" />

# ConR package

<br/><br/>

The `ConR` package aims at assisting the preliminary assessment of
species conservation status based on the International Union for
Conservation of Nature (IUCN) Red List Categories and Criteria. More
specifically, it helps users to calculate the population metrics related
to the IUCN criteria A, B, C and D and to assign one of the IUCN
categories (i.e. EX, EW, CR, EN, VU, NT, LC, DD). It was developed speed
up the assessment of hundreds to thousands of species at the same time.

The ideas behind this package, and its testing involved many people and
institutes.

<!-- See the original paper published in [Ecology and Evolution](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3704). -->

The package it not available anymore on CRAN because of the recent
retirement of various packages (rgdal, rgeos). The new version is
nevertheless available on github and will soon be available on CRAN.

To install the github version :

``` r
install.packages("devtools")
devtools::install_github("gdauby/ConR")
```

## Install R, ConR and dependent packages

**First step**. Install [R](https://cran.r-project.org/).

<!-- **Second step**. A proper way to work with R is to define a [working directory](https://bookdown.org/ndphillips/YaRrr/the-working-directory.html). If you are working with Rstudio, you can create the first time a [project](https://bookdown.org/ndphillips/YaRrr/projects-in-rstudio.html), which much simplify handling of scripts and data. -->

## Loading ConR

**Attach ConR package** This should be done everytime you open an R
session.

``` r
library(ConR)
```

**Help files** Any function in R is documented by a help file which can
be obtained by the following code:

``` r
?EOO.computing
?AOO.computing
?subpop.comp
?locations.comp
?criterion_A
?EOO.sensitivity
```

A mode detailed manual on how to use this package is available
[here](https://raw.githubusercontent.com/gdauby/ConR/devel/vignettes/articles/ConR.pdf).

## A quick example

Occurrence data are given as a `data.frame` whose first three columns
are, in this order, latitude, longitude and taxon name. Column names do
not matter, but their positions do. The package ships a small example
dataset:

``` r
data(dataset.ex)
MyData <- dataset.ex[!dataset.ex$tax %in% c("species_1", "species_2"), ]
MyData$tax <- as.character(MyData$tax)
str(MyData)
#> 'data.frame':    316 obs. of  3 variables:
#>  $ ddlat: num  0.75 3.57 1.18 3.24 4.09 ...
#>  $ ddlon: num  29.75 16.12 9.87 10.58 9.05 ...
#>  $ tax  : chr  "Psychotria minuta" "Psychotria minuta" "Psychotria minuta" "Psychotria minuta" ...
```

The extent of occurrence (EOO, in km<sup>2</sup>) is computed for every
taxon at once:

``` r
EOO.computing(MyData, show_progress = FALSE)
#>                      tax     eoo issue_eoo
#> 1     Berlinia bruneelii 2646614        NA
#> 2     Oncocalamus mannii  660638        NA
#> 3 Platycoryne guingangae 3437198        NA
#> 4      Psychotria minuta  763732        NA
```

So is the area of occupancy (AOO, in km<sup>2</sup>), here on the 2 km
grid recommended by IUCN:

``` r
AOO.computing(MyData, show_progress = FALSE)
#>                      tax aoo issue_aoo
#> 1     Berlinia bruneelii 404        NA
#> 2     Oncocalamus mannii 172        NA
#> 3 Platycoryne guingangae 148        NA
#> 4      Psychotria minuta  40        NA
```

And the number of subpopulations, using a 5 km circular buffer around
each occurrence:

``` r
subpop.comp(MyData, resol_sub_pop = 5, show_progress = FALSE)
#>                      tax subpop
#> 1     Berlinia bruneelii     88
#> 2     Oncocalamus mannii     29
#> 3 Platycoryne guingangae     34
#> 4      Psychotria minuta     10
```

These metrics, along with the number of locations, feed the assessment
functions `criterion_A()`, `criterion_B()`, `criterion_C()` and
`criterion_D()`. See the manual linked above for the full workflow.

## Funding

The development of this package was supported by:

- the European Union’s Horizon 2020 research and innovation program
  under the Marie Skłodowska-Curie grant agreement No 795114 ([THREAT
  project](https://cordis.europa.eu/project/id/795114))

- CESAB (Centre for the Synthesis and Analysis of Biodiversity).
  Research program of the FRB (Foundation for Research on Biodiversity)
  under the [RAINBIO
  project](https://gdauby.github.io/rainbio/index.html).

<!-- ## Acknowledgements -->

## Citation

G. Dauby & R. A. F. de Lima (2023). ConR: Computation of Parameters Used
in Preliminary Assessment of Species Conservation Status. R package
(version 2.1).

## See Also

Other R packages related to IUCN assessments:

- [`rredlist`](https://cran.r-project.org/web/packages/rredlist/rredlist.pdf)

- [`red`](https://cran.r-project.org/web/packages/red/red.pdf)

- [`redlistr`](https://cran.r-project.org/web/packages/redlistr/redlistr.pdf)
