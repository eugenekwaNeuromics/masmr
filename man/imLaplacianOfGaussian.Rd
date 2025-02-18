\name{imLaplacianOfGaussian}
\alias{imLaplacianOfGaussian}
\title{imLaplacianOfGaussian}
\description{
Calculates the difference between Gaussian blurs.
}
\usage{
imLaplacianOfGaussian(
    im,
    smallBlur,
    bigBlur
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ smallBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}.  }
\item{ bigBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}. Should be bigger than \code{smallBlur}.  }
}
\value{
  Processed image matrix.
}
\details{
  Useful for correcting uneven lighting and blob detection.
}
\seealso{
\link{getImageMetrics}
}
\examples{}
