\name{crossCorrelate2D}
\alias{crossCorrelate2D}
\title{crossCorrelate2D}
\description{
Calculates cross correlation (without any normalisation).
}
\usage{
crossCorrelate2D( referenceImageMatrix, queryImageMatrix, pad = FALSE )
}
\arguments{
\item{ referenceImageMatrix }{ 2D numeric matrix / image. }
\item{ queryImageMatrix }{ 2D numeric matrix / image. }
\item{ pad }{ Boolean. If \code{TRUE}, pads images with 0 values to have the same dimensions (dimension of image 1 + dimension of image 2 - 1). }
}
\value{
  A matrix of unnormalised cross-correlation coefficients between two images at all pixel locations.
}
\details{
  The default implementation is meant to be compatible with \code{registerImages}. It expects two numeric matrices of the same dimensions, and calculates cross correlation by calculating the product of \code{fft(Mean-centered Reference image)} and \code{fft(Flipped mean-centred Query Image)}, before calculating the inverse FFT and obtaining the real component. This is a modified version of the default approach to getting cross correlation: which would be the product of \code{fft(Mean-centered Reference image)} and \code{Conj(fft(Mean-centred Query Image))} -- note the use of the conjugate of the FFT and the absence of flipping. In our experience, the flipped version of this algorithm performs better for registration and stitching (i.e. better return of spot numbers, and higher FPKM correlation). 
  
  When \code{pad} is \code{TRUE}, we work with two images (that need not be of the same dimensions) and perform some padding with \code{0} values. Padding is often necessary for \code{stitchImages}.
}
\seealso{
\link{registerImages}
\link{stitchImages}
}
\examples{}
