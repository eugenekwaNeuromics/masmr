\name{imsBrightnessMin_MIP}
\alias{imsBrightnessMin_MIP}
\title{imsBrightnessMin_MIP}
\description{
Summarises image: finds background brightness. Designed to work on image stacks by first collapsing stacks into a single Z plane using maximum intensity projection.
}
\usage{
imsBrightnessMin_MIP(
    im,
    zDim = 3,
    ...
)
}
\arguments{
\item{ im }{ A 2-4D numeric array. }
\item{ zDim }{ Integer. Specifies the dimension that Z planes are found in. Passed to \code{maxIntensityProject}.  }
\item{ ... }{ Other functions passed from \code{getAnchorParams} to \code{imsBrightnessMin}.  }
}
\value{
  A numeric: minimum acceptable brightness.
}
\details{
  Image processing involves the following steps:

  \enumerate{
  \item{ If image is a stack, collapse Z planes using maximum intensity projection.}
  \item{ Gaussian blurring of image. }
  \item{ Calculation of Laplacian of Gaussians (LoG). }
  \item{ Quantile thresholding to separate blurred intensities between LoG>0 and LoG<0.}
  }
}
\seealso{
\link{getAnchorParams}
\link{imsBrightnessMin}
\link{maxIntensityProject}
}
\examples{
}
