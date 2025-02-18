\name{imsBrightnessMin}
\alias{imsBrightnessMin}
\title{imsBrightnessMin}
\description{
Summarises image: finds background brightness.
}
\usage{
imsBrightnessMin(
    im,
    smallBlur,
    bigBlur,
    ...
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ smallBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}.  }
\item{ bigBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}. Should be bigger than \code{smallBlur}.  }
\item{ ... }{ Other functions passed from \code{getAnchorParams}.  }
}
\value{
  A numeric: minimum acceptable brightness.
}
\details{
  Image processing involves the following steps:

  \enumerate{
  \item{ Gaussian blurring of image. }
  \item{ Calculation of Laplacian of Gaussians (LoG). }
  \item{ Quantile thresholding to separate blurred intensities between LoG>0 and LoG<0.}
  }
}
\seealso{
\link{getAnchorParams}
}
\examples{
}
