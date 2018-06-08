# ConR
ConR


# Install ConR

Install [R](https://cran.r-project.org/).
You may want to install the latest version from github by using the following code:

```
install.packages("devtools")
devtools::install_github("gdauby/ConR")
```


See the associated [vignette](https://CRAN.R-project.org/package=ConR/vignettes/my-vignette.html) for a overview of the different functions.


*__UPDATE June 2018__*

It is now possible to run the analysis in parallel which can fasten considerably the computation for large datasets (it does not a significant difference for small datasets or when mapping is not requested).

The example below run the evaluation in 4 cores in parallel.

```
data(Malagasy.amphibian)
MyResults <- IUCN.eval(Malagasy.amphibian, parallel=T, NbeCores=4, DrawMap=T)
```





