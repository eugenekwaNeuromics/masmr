\name{imForMask}
\alias{imForMask}
\title{imForMask}
\description{
Masks regions of interest in image.
}
\usage{
imForMask(
    im,
    smallBlur,
    minBlobSize,
    currentIteration,
    ...)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ smallBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}.  }
\item{ minBlobSize }{ Numeric (>0). Smallest acceptable blob size in pixels. }
\item{ currentIteration }{ Character. For keeping track of progress. }
\item{ ... }{ Other parameters passed to \code{imLowPass}, \code{imHessianDeterminant}, \code{findThreshold}, and \code{imLaplacianOfGaussian}. }
}
\value{
  Processed image matrix.
}
\details{

  Image processing involves the following steps:

  \enumerate{
  \item{  Clipping of image intensities.  }
  \item{  Low pass filtering of clipped image.  }
  \item{  Calculation of the determinant of the Hessian (detHess).  }
  \item{ Blurring of detHess (to account for large puncta: where edges have high detHess but centres have low.). }
  \item{ Elbow thresholding to find optimal detHess value that separates puncta from background. Threshold used to mask spots (i.e. detHess > threshold). }
  \item{ Identification of pixels brighter than their immediate background with Laplacian of Gaussian (LoG). }
  \item{ Summing of low-pass and high-pass filtered images, and normalisation of intensities between 0 and 1: i.e. creation of band-pass images. }
  \item{ Creating a mask where pixels of interest have LoG > 0, detHess > threshold, or are saturated. }
  \item{ Quantile thresholding that separates the median of masked and unmasked band-pass image values. Creates a mask of regions of interest (ROI). }
  \item{ Blobs that are too small (<\code{minBlobSize}) in ROI mask are filtered out. }
  }

 Masks pixels worth saving for subsequent decoding. Masks should comprise pixels nearby to puncta.

}
\seealso{
\link{getImageMetrics}
}
\examples{}
