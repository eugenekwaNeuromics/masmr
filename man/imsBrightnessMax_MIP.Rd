\name{imsBrightnessMax_MIP}
\alias{imsBrightnessMax_MIP}
\title{imsBrightnessMax_MIP}
\description{
Summarises image: finds saturation brightness. Designed to work on image stacks by first collapsing stacks into a single Z plane using maximum intensity projection.
}
\usage{
imsBrightnessMax_MIP(
    im,
    zDim = 3,
    ...
)
}
\arguments{
\item{ im }{ A 2-4D numeric array. }
\item{ zDim }{ Integer. Specifies the dimension that Z planes are found in. Passed to \code{maxIntensityProject}.  }
\item{ ... }{ Other functions passed from \code{getAnchorParams} to \code{imsBrightnessMax}.  }
}
\value{
  A numeric: maximum acceptable brightness.
}
\details{
  Image processing involves the following steps:

  \enumerate{
  \item{ If image is a stack, collapse Z planes using maximum intensity projection.}
  \item{ Gaussian blurring of image. }
  \item{ Background substraction of blurred image. }
  \item{ Calculation of determinant of Hessian (detHess). }
  \item{ Masking of pixels with high detHess and intensity values.  }
  \item{ Quantile thresholding to separate blurred intensities between best spots and all other pixels.}
  }
}
\seealso{
\link{getAnchorParams}
\link{imsBrightnessMax}
\link{maxIntensityProject}
}
\examples{
}
