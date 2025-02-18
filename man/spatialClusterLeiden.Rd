\name{spatialClusterLeiden}
\alias{spatialClusterLeiden}
\title{spatialClusterLeiden}
\description{
Clusters proximal pixels (using Leiden clustering), and assigns a single pixel as cluster centroid.
}
\usage{
spatialClusterLeiden(
    spotcalldf,
    leidenResolution = 0.01,
    minNeighbours = 3,
    maxInterSpotDistancePixels = 5,
    distanceMetric = 'COS',
    seed = 12345,
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ leidenResolution }{ Numeric (>0). Resolution for Leiden clustering: see \code{igraph::cluster_leiden()}. }
\item{ minNeighbours }{ Numeric (>0). Minimum number of pixels per cluster for cluster to be considered valid. }
\item{ maxInterSpotDistancePixels }{ Numeric (>0). Maximum allowable distance (in pixels) between pixels in a given cluster. }
\item{ distanceMetric }{ Character. References a column in the supplied dataframe that contains the desired distance metric. The pixel with the lowest specified distance metric for its cluster will be labelled the 'centroid' of that cluster. }
\item{ seed }{ Integer. Seed for pseudo-random number generation. As default, read from \code{params$seed}. If unavailable, taken from \code{seed} value in this function. }
\item{ params }{  The parameters environment. }
}
\value{
  A dataframe with the same the number of rows as \code{spotcalldf}. Columns returned: (1) \code{CLUSTER}, the cluster identifier; (2) \code{CLUSTER_SIZE}, number of pixels in the cluster; (3) \code{CENTROID}, if the pixel is considered the cluster's 'centroid'. When clustering is not possible (e.g. too few pixels within acceptable distance), \code{CLUSTER} returns 0, and \code{CLUSTER_SIZE} returns 1.
}
\details{
}
\seealso{
}
\examples{
clusterInfo <- spatialClusterLeiden(spotcalldf)
spotcalldf <- cbind(spotcalldf, clusterInfo)
}
