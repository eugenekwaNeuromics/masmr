\name{crossCorrelate2D}
\alias{crossCorrelate2D}
\title{crossCorrelate2D}
\description{
Calculates cross correlation.
}
\usage{
crossCorrelate2D( referenceImageMatrix, queryImageMatrix, normalized = FALSE, pad = TRUE )
}
\arguments{
\item{ referenceImageMatrix }{ 2D numeric matrix / image. }
\item{ queryImageMatrix }{ 2D numeric matrix / image. }
\item{ normalized }{ Boolean. If \code{TRUE}, returns normalised cross correlation values. }
\item{ pad }{ Boolean. If \code{TRUE}, pads images to have the same dimensions. If \code{FALSE}, images need to have the same dimensions. }
}
\value{
  A matrix of normalised/unnormalised cross-correlation coefficients between two images at all pixel locations.
}
\examples{}
