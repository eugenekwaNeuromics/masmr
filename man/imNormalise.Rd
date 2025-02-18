\name{imNormalise}
\alias{imNormalise}
\title{imNormalise}
\description{
Normalises the intensity values of an image.
}
\usage{
imNormalise(
    im,
    floorVal,
    ceilVal,
    floorQuant,
    ceilQuant
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. Intensities normalised between 0 and 1. }
\item{ floorVal }{ Numeric (>0). Floor of dynamic range. Values below this floor set to 0.  }
\item{ ceilVal }{  Numeric (<1). Ceiling of dynamic range. Values above this ceiling set to 1. }
\item{ floorQuant }{ Numeric (>0). Floor of dynamic range, expressed as a quantile.  }
\item{ ceilQuant }{  Numeric (<1). Ceiling of dynamic range, expressed as a quantile. }
}
\details{
  Image is min-max normalised between \code{floorVal} and \code{ceilVal}. When \code{floorQuant} specified, quantile for intensities above \code{floorVal} serves as floor instead. Likewise, if \code{ceilQuant} specified, quantile for intensities below \code{ceilVal} serves as ceiling.

  When unspecified, \code{floorVal} is automatically assigned the minimum intensity in the image; and \code{ceilVal}, the maximum.
}
\value{
  Normalised image matrix with values between 0 and 1.
}
\seealso{
}
\examples{}
