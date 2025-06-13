\name{getPNNMatrix}
\alias{getPNNMatrix}
\title{getPNNMatrix}
\description{
Creates a Point's Nearest Neighbour matrix.
}
\usage{
getPNNMatrix(
  x,
  y,
  label,

  delaunayTriangulation = TRUE,
  delaunayDistanceThreshold = 20,
  delaunayNNDegrees = c(1:3),

  euclideanDistances = NULL,

  seed = 12345,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{ x }{ Numeric vector. X coordinates of point. }
\item{ y }{ Numeric vector. Y coordinates of point. }
\item{ label }{ String, integer, or factor vector. Indicates the identity of each point (e.g. gene name, cell type). }

\item{ delaunayTriangulation }{ Boolean (default \code{TRUE}). Determines whether to use the Delaunay Triangulation approach to finding neighbours (see details). If \code{FALSE}, then \code{euclideanDistances} must be specified. }
\item{ delaunayDistanceThreshold }{ Numeric. Determines the maximum allowable distance for a connection in the Delaunay Triangulation graph. Edges longer than this distance are dropped. Units depend on the units of \code{x} and \code{y}.}
\item{ delaunayNNDegrees }{ Integer vector. Determines the degree bins to find nearest neighbours. Bins are inclusive: i.e. "1" contains just 1st order neighbours, while "2" contains 1st and 2nd order neighbours. }

\item{ euclideanDistances }{ Numeric vector. Determines the distance bins to find nearest neighbours. Bins are inclusive: i.e. a 10-distance bin contains a 5-distance bin. Units depend on the units of \code{x} and \code{y}.  }

\item{ seed }{ Integer. Seed for pseudo-random number generation. }
\item{ verbose }{ Boolean (default \code{TRUE}). If \code{TRUE}, returns messages from each function. Otherwise, messages not returned (but warnings and stop messages still reported). }

\item{ ... }{ Other parameters passed to the \code{tripack::tri.mesh} function (for Delaunay Triangulation). }
}
\value{
  A PNNMatrix: with points in rows, and neighbourhoods in columns. Tabulates the frequency that a point of label X is a neighbour of a given reference point.
}
\details{
  The default approach is to perform a Delaunay Triangulation (DT) to connect all points into a simplicial mesh. Nearest neighbours of each point -- for each bin in \code{delaunayNNDegrees} -- are then found and tabulated according to their label.

  First order neighbours (\code{delaunayNNDegrees=1}) are those directly connected to the reference spot. Second order neighbours are those connected to first order neighbours. When \code{delaunayNNDegrees=2}, both first and second order neighbours are considered. Same logic applies for other degrees.

  Alternatively, bins can be determined according to physical distances with \code{euclideanDistances}. However, this approach is slower.

  When \code{delaunayTriangulation} is \code{TRUE} and  \code{euclideanDistances} is specified, then both the DT and Euclidean methods are employed. Columns are labelled with a "DT" prefix for DT, and "EUD" for Euclidean. Columns will also contain the distance / degree of the bin used.
}
\seealso{
\link{ tripack::tri.mesh }
}
\examples{
nearestNeighbourMatrix <- getPNNMatrix(

  x = spotcalldf$Xm, # The x coordinate
  y = spotcalldf$Ym, # The y coordinate
  label = spotcalldf$g, # The category that a given spot belongs to

  delaunayTriangulation = T, # Whether to use Delaunay Triangulation approach
  delaunayDistanceThreshold = 20, # Ignore links longer than this distance (Xm is in microns, so this is 20 microns)
  delaunayNNDegrees = c(1:3), # Consider 1st, 1st + 2nd, and 1st + 2nd + 3rd order neighbours

  euclideanDistances = c(5, 10) # Consider spots <5, and <10 microns away
)
}
