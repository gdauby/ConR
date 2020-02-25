\name{locations.comp}
\alias{locations.comp}
\title{Number of locations}
\description{
  Compute number of locations for multiple taxa
}
\usage{
  locations.comp(XY, method = "fixed_grid",
                           nbe_rep = 0,
                           protec.areas = NULL, 
                           Cell_size_locations = 10,
                           method_protected_area="no_more_than_one", 
                           ID_shape_PA="WDPA_PID", 
                           Rel_cell_size=0.05,
                           parallel=FALSE,
                           NbeCores=2)
}
\arguments{
\item{XY}{\code{dataframe} see Details}
  
\item{method}{a character string,  indicating the method used for estimating the number of locations. 
"fixed_grid" or "sliding scale". See details. By default, it is "fixed_grid"}
  
\item{nbe_rep}{a numeric, indicate the number of raster with random starting position for estimating the number of location. By default, it is 0 but some minimal translation of the raster are still done}

\item{protec.areas}{a \code{SpatialPolygonsDataFrame}, shapefile with protected areas. 
If provided, this will be taken into account for calculating number of location (see Details and \code{method_protected_area}). 
By default, no shapefile is provided}

\item{Cell_size_locations}{a numeric, value indicating the grid size in kilometres used for estimating the number of location. By default, equal to 10}

\item{method_protected_area}{a character string. By default is "no_more_than_one"", which means occurrences within protected areas (if provided)
will not be taken into account for estimating the number of locations following the grid system, see Details. By default, it is "no_more_than_one"}

\item{ID_shape_PA}{a character string, indicating the field name of \code{protec.areas} with ID of the \code{SpatialPolygonsDataFrame} of protected areas}

\item{Rel_cell_size}{a numeric, if \code{method_locations="sliding scale"}, \code{Cell_size_locations} is ignored and the resolution is 
given by the maximum distance separating two occurrences multiplied by \code{Rel_cell_size}. By default, it is 0.05}

\item{parallel}{a logical. Wether running in parallel. By default, it is FALSE}

\item{NbeCores}{an integer. Register the number of cores for parallel execution. By default, it is 2}
  
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
  A list with one list for each species containing [[1]]SpatialPolygonDataframe and [[2]]vector of the number of location.
  
}
\author{
  Gilles Dauby
  
  
  \email{gildauby@gmail.com}
}

\examples{
data(dataset.ex)
\dontrun{
locations <- locations.comp(dataset.ex)

}

}
