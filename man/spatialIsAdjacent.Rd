\name{spatialIsAdjacent}
\alias{spatialIsAdjacent}
\title{spatialIsAdjacent}
\description{
Determines if coordinates are adjacent or overlapping a set of query coordinates.
}
\usage{
spatialIsAdjacent(
    spotcalldf,
    queryLocations,
    gridDistance = 1
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. Coordinates extracted from \code{WX} and \code{WY} columns. Alternatively, also accepts a complex vector representing coordinates (e.g. X + 1i * Y). }
\item{ queryLocations }{ If a single character or numeric, and a dataframe supplied as \code{spotcalldf}, will assume \code{queryLocations} to refer to a given column name or column position in said dataframe: \code{ spotcalldf[,queryLocations]} must then return a boolean vector. Alternatively, also accepts a complex vector representing coordinates (e.g. X + 1i * Y). }
\item{ gridDistance }{ Numeric (>0). Size of search region in pixels. If non-integer, rounded up to nearest whole number. }
}
\value{
  A boolean vector for each pixel \code{spotcalldf}. Vector indicates if given pixel is at most \code{gridDistance} away from coordinates specified by \code{queryLocations}.
}
\details{
}
\seealso{
}
\examples{
isAdjacentToBlanks <- spatialHNNDist(spotcalldf, 'BLANK')
spotcalldf$NEXT_TO_BLANKS <- isAdjacentToBlanks
}
