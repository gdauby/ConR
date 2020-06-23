
<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/gdauby/ConR.svg?branch=master)](https://travis-ci.com/gdauby/ConR)
<!-- badges: end --> [![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/ConR)](https://www.r-pkg.org/pkg/ConR)

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

## Funding

The development of this package was supported by:

  - the European Union’s Horizon 2020 research and innovation program
    under the Marie Skłodowska-Curie grant agreement No 795114.

## Acknowledgements

## See Also

Other R packages related to IUCN assessments:

  - `rredlist`
    (<https://cran.r-project.org/web/packages/rredlist/rredlist.pdf>)

  - `red` (<https://cran.r-project.org/web/packages/red/red.pdf>)

  - `redlistr`
    (<https://cran.r-project.org/web/packages/redlistr/redlistr.pdf>)
