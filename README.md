
<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/gdauby/ConR.svg?branch=master)](https://travis-ci.com/gdauby/ConR)
<!-- badges: end -->

The ideas behind this package, and its testing involved many people and
institutes. See the original paper published in [Ecology and
Evolution](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3704).

## Install R, ConR and dependent packages

**First step**. Install [R](https://cran.r-project.org/).

**Second step**. A proper way to work with R is to define a [working
directory](https://bookdown.org/ndphillips/YaRrr/the-working-directory.html).
If you are working with Rstudio, you can create the first time a
[project](https://bookdown.org/ndphillips/YaRrr/projects-in-rstudio.html),
which much simplify handling of scripts and data.

**Third step**. Install ConR package (one of the two following) :

To install the version on CRAN

    install.packages("ConR")

To install the version under development.

    install.packages("devtools")
    devtools::install_github("gdauby/ConR")

# Loading ConR

**Attach ConR package** This should be done everytime you open an R
session.

    library(ConR)

**Help files** Any function in R is documented by a help file which can
be obtained by the following code:

    ?EOO.computing
    ?IUCN.eval
    ?map.res
    ?subpop.comp
    
## Authors

  Gilles Dauby <gildauby@gmail.com>
  
  Renato A. F. de Lima <renato.lima@naturalis.nl>

## References

Dauby, G, Stévart, T, Droissart, V, et al. ConR : An R package to assist large‐scale     multispecies preliminary conservation assessments using distribution data. Ecol     Evol. 2017; 7: 11292– 11303. https://doi.org/10.1002/ece3.3704

IUCN 2019. Guidelines for Using the IUCN Red List Categories and
    Criteria. Version 14. Standards and Petitions Committee. Downloadable from:
    http://www.iucnredlist.org/documents/RedListGuidelines.pdf.

## Funding

The development of this package was supported by: 

  - the European Union’s Horizon 2020 research and innovation program under the Marie Skłodowska-Curie grant agreement No 795114.

## Acknowledgements


## See Also

Other R packages related to IUCN assessments:

  - `rredlist` (https://cran.r-project.org/web/packages/rredlist/rredlist.pdf)
  
  - `red` (https://cran.r-project.org/web/packages/red/red.pdf)
  
  - `redlistr` (https://cran.r-project.org/web/packages/redlistr/redlistr.pdf)
