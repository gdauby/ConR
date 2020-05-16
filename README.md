
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
