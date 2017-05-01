\name{EOO.computing}
\alias{EOO.computing}
\title{Extent of Occurrences}
\description{
Compute extent of occurrences (EOO) for multiple taxa in square kilometres using \code{\link[geosphere]{geosphere}} package
and provide \code{SpatialPolygons} used for EOO computation 
}
\usage{
EOO.computing(XY, exclude.area=FALSE, buff_width=0.1, country_map=NULL, 
	Name_Sp="Species1", export_shp=FALSE,
                          alpha=1, buff.alpha=0.1, method.range="convex.hull", 
						  write_shp=FALSE, 
						  method.less.than3="not comp", 
						  write.results=TRUE, 
                          file.name="EOO.results", verbose=TRUE)
}
\arguments{
\item{XY}{\code{dataframe} see Details}
	
\item{exclude.area}{a logical, if TRUE, areas outside of \code{country_map} are cropped of \code{SpatialPolygons} used for calculating EOO. By default is FALSE}

\item{country_map}{a \code{SpatialPolygonsDataFrame} or \code{SpatialPolygons} showing for example countries or continent borders. This shapefile will be used for cropping the \code{SpatialPolygons}l if exclude.area is TRUE}

\item{export_shp}{a logical, whether shapefiles should be exported or not, see Value. By default is FALSE}

\item{write_shp}{a logical, if TRUE, export \code{SpatialPolygons} used for EOO computation as ESRI shapefiles in the working directory. By default is FALSE}

\item{alpha}{a numeric, if \code{method.range} is "alpha.hull", value of alpha of the alpha hull, see \code{\link[alphahull]{ahull}}. By default is 1}

\item{buff.alpha}{a numeric, if \code{method.range} is "alpha.hull", define the buffer in decimal degree added to alpha hull. By default is 0.1}

\item{method.range}{a character string, "convex.hull" or "alpha.hull". By default is "convex.hull"}

\item{Name_Sp}{a character string, if \code{XY} is for one taxon and field containing taxon names is not provided, 
	this item provide taxon name. By default is "Species1"}

\item{buff_width}{a numeric. For a specific case where all points of a taxa are on a straight line, 
	see Details. By default is 0.1}

\item{method.less.than3}{a character string. If equal to "arbitrary", will give a value to species with two unique occurrences, see Details. By default is "not comp"}

\item{write.results}{a logical. If TRUE, results will be exported in the working environment as a csv file. By default is TRUE}

\item{file.name}{a character string. Name file for exported results in csv file. By default is "EOO.results"}

\item{verbose}{a logical. Wether R should report information on progress}
	}
	
	
\details{
\strong{Input} as a \code{dataframe} should have the following structure:

\strong{It is mandatory to respect field positions, but field names do not matter}

\tabular{ccc}{
[,1] \tab ddlat \tab numeric, latitude (in decimal degrees)\cr
[,2] \tab ddlon \tab numeric, longitude (in decimal degrees)\cr
[,3] \tab tax \tab character or factor, optinal field with taxa names\cr
}

\strong{Important notes:}

EOO will only be computed if there is at least three unique occurrences unless \code{method.less.than3} is put to "arbitrary".
In that specific case, EOO for species with two unique occurrences will be equal to 
Dist*Dist*0.1 
where Dist is the distance in kilometres separating the two points.

For the very specific (and infrequent) case where all occurrences 
are localized on a straight line (in which case EOO would be null),
EOO is estimated by the area of polygon surrounding this straight line with a buffer of \code{buff.alpha} decimal degree.
There is a warning when this happen.
}
\value{
If \code{export_shp} is FALSE, a \code{dataframe} with one field containing EOO in square kilometres.
 \code{NA} is given when EOO could not be computed because there is less than three unique occurrences 
 (or two if \code{method.less.than3} is put to "arbitrary").

If \code{export_shp} is TRUE, a \code{list} with:
\enumerate{
	\item EOO in square kilometres
	\item \code{SpatialPolygons} used for EOO computation
}

}
\author{
Gilles Dauby


\email{gildauby@gmail.com}
}


\seealso{
  \code{\link[alphahull]{ahull}}
  
  
  \code{\link[speciesgeocodeR]{CalcRange}}
}

\references{
Gaston & Fuller 2009 The sizes of species'geographic ranges, Journal of Applied Ecology, 49 1-9

}
\examples{
data(dataset.ex)
data(land)
\dontrun{
EOO <- EOO.computing(dataset.ex)

## This exclude areas outside of land (i.e. ocean) for EOO computation
EOO <- EOO.computing(dataset.ex, 
  exclude.area=TRUE, country_map=land)
}

}
