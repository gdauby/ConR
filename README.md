
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

    install.packages("devtools")
    devtools::install_github("gdauby/ConR")

## Install R, ConR and dependent packages

**First step**. Install [R](https://cran.r-project.org/).

<!-- **Second step**. A proper way to work with R is to define a [working directory](https://bookdown.org/ndphillips/YaRrr/the-working-directory.html). If you are working with Rstudio, you can create the first time a [project](https://bookdown.org/ndphillips/YaRrr/projects-in-rstudio.html), which much simplify handling of scripts and data. -->

## Loading ConR

**Attach ConR package** This should be done everytime you open an R
session.

    library(ConR)

**Help files** Any function in R is documented by a help file which can
be obtained by the following code:

    ?EOO.computing
    ?AOO.computing
    ?subpop.comp
    ?locations.comp
    ?criterion_A
    ?EOO.sensitivity

A manual for using the package is available
[here](https://raw.githubusercontent.com/gdauby/ConR/devel/vignettes/articles/ConR.html).

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
(version 2.1.0).

## See Also

Other R packages related to IUCN assessments:

- [`rredlist`](https://cran.r-project.org/web/packages/rredlist/rredlist.pdf)

- [`red`](https://cran.r-project.org/web/packages/red/red.pdf)

- [`redlistr`](https://cran.r-project.org/web/packages/redlistr/redlistr.pdf)
