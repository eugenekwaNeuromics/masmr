\name{spatialHNNDist}
\alias{spatialHNNDist}
\title{spatialHNNDist}
\description{
Measures shortest homotypic neighbour distances.
}
\usage{
spatialHNNDist(
    spotcalldf,
    seed = 12345,
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ seed }{ Integer. Seed for pseudo-random number generation. As default, read from \code{params$seed}. If unavailable, taken from \code{seed} value in this function. }
\item{ params }{  The parameters environment. }
}
\value{
  A numeric vector of nearest neighbour distances for each pixel (row) in supplied dataframe. For each given pixel, nearest neighbour must be the same gene (\code{g} column) - i.e. be a 'homotypic' neighbour.
}
\details{
}
\seealso{
}
\examples{
hnndist <- spatialHNNDist(spotcalldf)
spotcalldf$HNNDIST <- hnndist
}
