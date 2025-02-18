\name{imLowPass}
\alias{imLowPass}
\title{imLowPass}
\description{
Low pass filtering of image.
}
\usage{
imLowPass(
    im,
    smallBlur,
    currentIteration,
    ...
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ smallBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}.  }
\item{ currentIteration }{ Character. For keeping track of progress. }
\item{ ... }{ Other parameters passed to \code{imNormalise}. }
}
\value{
  Processed image matrix.
}
\details{
  Image processing involves the following steps:

  \enumerate{
  \item{Identification of background from original image.}
  \item{Gaussian blurring of original image.}
  \item{Autobrightening of image to maximise separation between median intensity values in foreground from background in blurred image.}
  }

  Useful for separating foreground from background.
}
\seealso{
\link{getImageMetrics}
}
\examples{}
