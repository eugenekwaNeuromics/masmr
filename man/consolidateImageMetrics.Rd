\name{consolidateImageMetrics}
\alias{consolidateImageMetrics}
\title{consolidateImageMetrics}
\description{
Consolidates image metrics into a single dataframe.
}
\usage{
consolidateImageMetrics(
    maskList,
    bitFloor,
    bitCeil,
    imMetrics = get('imMetrics', envir = globalenv()),
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ maskList }{ A list of masks. If unspecified, read from \code{imMetrics$MASK}. }
\item{ bitFloor }{ Numeric (>=0). Minimum number of images a pixel has been flagged as ON for it to be saved in the dataframe. If unspecified, 1. }
\item{ bitCeil }{ Numeric (>0). Maximum number of images a pixel has been flagged as ON for it to be saved in the dataframe. If unspecified, \code{params$nbits * 2}. }
\item{ params }{  The parameters environment. }
\item{ imMetrics }{ The image metrics environment (see \code{getImageMetrics}). }
}
\value{
  A dataframe, with pixels as rows and image metrics as columns.
}
\details{
  After performing this step, all subsequent processing will be formed on this dataframe. It is thus recommended to clean the environment to remove unneeded objects (while keeping the new dataframe): see Examples below.

  For memory conservation purposes, a \code{maskList} (or \code{imMetrics$MASK}) should be provided. Otherwise, metrics for all pixels will be returned.

  Column names of the dataframe:

  \itemize{
  \item{WX: X coordinate in current FOV (pixels).}
  \item{WY: Y coordinate in current FOV (pixels).}
  \item{IDX: An integer identifier unique to each pixel in a given FOV.}
    \item{METRIC_XX: Several columns with the format of metric name as prefix (e.g. ORIG) and a number indicating bit index.}
  }

  Coordinates have been corrected with reference to \code{params$shifts} (see \code{registerImages}).
}
\seealso{
}
\examples{
spotcalldf <- consolidateImageMetrics()
cleanUp( c('spotcalldf', cleanSlate) )
}
