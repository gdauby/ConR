\name{subpop.comp}
\alias{subpop.comp}
\title{Subpopulations computation following Circular buffer method}
\description{
Compute subpopulations follwing the Circular buffer method (see reference)
}
\usage{
subpop.comp(XY,Resol_sub_pop)
}
\arguments{
\item{XY}{Dataframe, coordinates in decimal degrees of one species, latitude in first colum,, longitude in second column}

\item{Resol_sub_pop}{Numeric, radius of circles in kilometres}

}

\details{
Circular buffer method of Rivers et al. 2010
Each unique occurrence is buffered with a circle of a radius of \code{Resol_sub_pop km}. 
Overlapping circles are merged to form a single subpopulation, while non-overlapping circles are considered as
separate subpopulations.
}

\value{

A \code{list}, two elements:

\enumerate{
\item{first element} {The number of subpopulation}
\item{second element} {A shapefile of the subpopulations}
}

}
\author{
Gilles Dauby
}

\references{
Rivers CM, Bachman SP & Meagher TR, 2010, Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data, Biodiversity Conservation 19:2071-2085
}

\examples{
\dontrun{
data(dataset.ex)
data(land)

SUB <- subpop.comp(dataset.ex, Resol_sub_pop=25)
plot(SUB[["Platycoryne guingangae"]][["subpop.poly"]], col="red")
plot(land, add=TRUE)
}
}
