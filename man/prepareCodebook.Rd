\name{prepareCodebook}
\alias{prepareCodebook}
\title{prepareCodebook}
\description{
Prepares codebook for downstream use.
}
\usage{
prepareCodebook(
    params = get('params', envir = globalenv()),
    exhaustiveBlanks = T,
    hammingDistanceThreshold = 2
  )
}
\arguments{
\item{ params }{ The parameters environment. }
\item{ exhaustiveBlanks }{ Boolean. Determines if more blanks are to be appended to existing codebook.}
\item{ hammingDistanceThreshold }{ Numeric. Minimum Hamming distance a proposed blank code should have to existing codes. }
}
\value{
  An ordered codebook, saved as a \code{ORDEREDCODEBOOK.csv.gz} file. Also updates the \code{params} environment.
}
\details{
}
\seealso{
}
\examples{}
