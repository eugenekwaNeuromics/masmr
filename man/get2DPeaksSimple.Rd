\name{get2DPeaksSimple}
\alias{get2DPeaksSimple}
\title{get2DPeaksSimple}
\description{
Finds local peaks in a 2D matrix.
}
\usage{
get2DPeaksSimple(
  imRaster,
  searchspace = c('queen', 'bishop', 'rook'),
  verbose = FALSE )
}
\arguments{
\item{ imRaster }{ 2D numeric matrix / image. }
\item{ searchspace }{ Character. Accepts \code{queen} (default), \code{bishop}, \code{rook}, or their first letter.}
\item{ verbose }{ Boolean. If \code{TRUE}, returns details. }
}
\value{ A matrix with local peaks set to 1, and 0 otherwise. }
\details{
For \code{searchSpace}: If \code{queen} or \code{q} (default), identifies peaks by examining all eight pixels surrounding the reference pixel. If \code{bishop} or \code{b}, searches only four pixels diagonal to the reference pixel. If \code{rook} or \code{r}, searches only four pixels non-diagonal to reference pixel (i.e. left-right, top-down).
}
\examples{}
