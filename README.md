
# How to use ConR (for beginners in R)

## Install R, ConR and dependent packages

**First step**. Install [R](https://cran.r-project.org/).

**Second step**. A proper way to work with R is to define a working
directory. I advise to create a directory near the root (easier to
handle), for example, here, named *Your\_R\_directory* Now you can set
your working directory by the following code:

    setwd("C:/Your_R_directory/"")

**Third step**. Install ConR package (one of the two following) :

To install the version on CRAN

    install.packages("ConR")

To install the version under development.

    install.packages("devtools")
    devtools::install_github("gdauby/ConR")

The package also comes with
[vignette](https://cran.r-project.org/web/packages/ConR/vignettes/my-vignette.html)
providing more examples.

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
