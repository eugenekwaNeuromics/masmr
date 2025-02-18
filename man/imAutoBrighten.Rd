\name{imAutoBrighten}
\alias{imAutoBrighten}
\title{imAutoBrighten}
\description{
Automatically brightens/darkens an image.
}
\usage{
imAutoBrighten(
    im,
    floorVal,
    ceilVal
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. Intensities normalised between 0 and 1. }
\item{ floorVal }{ Numeric (>0). Floor of dynamic range.  }
\item{ ceilVal }{  Numeric (<1). Ceiling of dynamic range. }
}
\details{
  Image matrix raised to power \code{x}, where \code{x} maximises the value of \code{ceilVal^x - floorVal^x}.

  Useful for separating two values (e.g. foreground-background), and thus has uses for image registration purposes. Image decoding may or may not benefit from this transformation - careful application recommended.
}
\value{
  Brightened/darkened image matrix.
}
\seealso{
}
\examples{}
