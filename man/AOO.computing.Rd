\name{AOO.computing}
\alias{AOO.computing}
\title{Area of occupancy}
\description{
Compute areas of occupancy (AOO) for multiple taxa in square kilometres
}
\usage{
AOO.computing(XY, Cell_size_AOO = 2, nbe.rep.rast.AOO = NULL, parallel = FALSE, NbeCores = 2, show_progress=FALSE)
}
\arguments{
\item{XY}{\code{dataframe} see Details}

\item{Cell_size_AOO}{a numeric, value indicating the grid size in decimal degrees used for estimating Area of Occupancy.  By default, equal to 2}

\item{nbe.rep.rast.AOO}{a numeric, indicate the number of raster with random starting position for estimating the AOO. By default, it is NULL but some minimal translation of the raster are still done}

\item{parallel}{a logical. Wether running in parallel. By default, it is FALSE}

\item{NbeCores}{an integer. Register the number of cores for parallel execution. By default, it is 2}

\item{show_progress}{a logical. Whether a bar showing progress in computation should be shown. By default, it is FALSE}


	}
	
	
\details{
\strong{Input} as a \code{dataframe} should have the following structure:

\strong{It is mandatory to respect field positions, but field names do not matter}

\tabular{ccc}{
[,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
[,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
[,3] \tab tax \tab character or factor, optinal field with taxa names\cr
}

}
\value{
A vector

}
\author{
Gilles Dauby


\email{gildauby@gmail.com}
}


\references{
Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9

}
\examples{
data(dataset.ex)
\dontrun{
AOO <- AOO.computing(dataset.ex)

}

}
