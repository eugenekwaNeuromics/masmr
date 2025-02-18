\name{imHessianDeterminant}
\alias{imHessianDeterminant}
\title{imHessianDeterminant}
\description{
Calculates the determinant of the Hessian.
}
\usage{
imHessianDeterminant(
    im,
    smallBlur
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ smallBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}.  }
}
\value{
  Processed image matrix.
}
\details{
  Useful for blob detection.
}
\seealso{
\link{getImageMetrics}
}
\examples{}
