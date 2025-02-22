\name{decode}
\alias{decode}
\title{decode}
\description{
Performs decoding of image list.
}
\usage{
decode(
    spotcalldf,
    nCores,
    currentFovName,
    decodeMetric = 'DECODE',
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ nCores }{ Integer. Number of threads. Defaults to -1 from total number of detected cores.}
\item{ currentFovName }{ Character. Name of current FOV. If not specified, will first attempt to read \code{params$current_fov} for same information. }
\item{ decodeMetric }{ Character. A previously calculated image metric that decoding will be performed against (checked against \code{params$imageMetrics}). }
\item{ params }{  The parameters environment. }
}
\value{
  An updated dataframe where pixels are now annotated with gene information.
}
\details{

Function removes rows where the decoding metric is always 0 across all images, and calculates various distance metrics to the code vector in \code{params$ordered_codebook}. These distance metrics are appended as columns to \code{spotcalldf}:

\itemize{
\item{  COS: Cosine distance. }
\item{  COSLAB: Decoded gene according to nearest \code{COS}.}
\item{  EUC: Euclidean distance. }
\item{  EUCLAB: Decoded gene according to nearest \code{EUC}. }
\item{  L2E: L2-norm euclidean distance. }
\item{  L2ELAB: Decoded gene according to nearest \code{L2E}. }
}

The corresponding values for the \strong{next} closest gene are similarly provided. These columns are similar to above but with an additional \code{B} prefix (e.g. \code{BCOS} and \code{BCOSLAB}).

Additionally, global coordinate position information is also appended:

\itemize{
\item{  Xm: X coordinate in microns. }
\item{  Ym: Y coordinate in microns. }
\item{  fov: Current FOV. }
\item{  g: Gene name (determined from \code{COSLAB}).  }
}

}
\seealso{
\link{ subsetDF }
\link{ scoreDeltaF }
}
\examples{}
