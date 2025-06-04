\name{imsBrightnessMax}
\alias{imsBrightnessMax}
\title{imsBrightnessMax}
\description{
Summarises image: finds saturation brightness.
}
\usage{
imsBrightnessMax(
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
  A numeric: maximum acceptable brightness.
}
\details{
  Image processing involves the following steps:

  \enumerate{
  \item{ Gaussian blurring of image. }
  \item{ Background substraction of blurred image. }
  \item{ Calculation of determinant of Hessian (detHess). }
  \item{ Masking of pixels with high detHess and intensity values.  }
  \item{ Quantile thresholding to separate blurred intensities between best spots and all other pixels.}
  }
  
  For image stacks, use \code{imsBrightnessMax_MIP} instead.
}
\seealso{
\link{getAnchorParams}
\link{imsBrightnessMax_MIP}
}
\examples{
}
