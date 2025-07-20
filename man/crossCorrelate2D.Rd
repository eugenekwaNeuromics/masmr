\name{crossCorrelate2D}
\alias{crossCorrelate2D}
\title{crossCorrelate2D}
\description{
Calculates cross correlation (without any normalisation).
}
\usage{
crossCorrelate2D( referenceImageMatrix, queryImageMatrix, forStitch = FALSE )
}
\arguments{
\item{ referenceImageMatrix }{ 2D numeric matrix / image. }
\item{ queryImageMatrix }{ 2D numeric matrix / image. }
\item{ forStitch }{ Boolean. If \code{TRUE}, implements the version of this algorithm designed to work with \code{stitchImages} (see details). }

}
\value{
  A matrix of unnormalised cross-correlation coefficients between two images at all pixel locations.
}
\details{
  The default implementation is meant to be compatible with \code{registerImages}. It expects two numeric matrices of the same dimensions, and calculates cross correlation by calculating the product of \code{fft(Mean-centered Reference image)} and \code{Conj(fft(Mean-centred Query Image))}, before calculating the inverse FFT and obtaining the real component. This is the typical approach for obtaining cross-correlation.
  
  When \code{forStitch} is \code{TRUE}, we work with two images (that need not be of the same dimensions) and perform some padding with \code{0} values. We do not obtain the conjugate of the mean-centred query image (i.e. \code{Conj()}) as before. Instead, we first flip the mean-centred query image along both axes, then pad, then calculate the discrete Fourier transform. We multiply this by the \code{fft(Mean-centered Padded Reference image)}, before calculating the inverse FFT and obtaining the real component. 
  
  Flipping of the query image typically returns comparable results to obtaining the conjugate of the unflipped image; but works better with our implementation of \code{stitchImages} -- which performs some image cropping for better memory management under the hood.
}
\seealso{
\link{registerImages}
\link{stitchImages}
}
\examples{}
