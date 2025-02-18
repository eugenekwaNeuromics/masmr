\name{imForDecode}
\alias{imForDecode}
\title{imForDecode}
\description{
Default image processing to return images optimised for decoding.
}
\usage{
imForDecode(
    im,
    smallBlur,
    currentIteration,
    ...)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ smallBlur }{ Numeric (>0). Standard deviation for isotropic blurring: see \link{imager::isoblur}.  }
\item{ currentIteration }{ Character. For keeping track of progress. }
\item{ ... }{ Other parameters passed to \code{imHessianDeterminant}, \code{findThreshold}, and \code{imLaplacianOfGaussian}. }
}
\value{
  Processed image matrix.
}
\details{

  Image processing involves the following steps:

  \enumerate{
  \item{  Clipping of image intensities.  }
  \item{  Calculation of the determinant of the Hessian (detHess).  }
  \item{ Blurring of detHess (to account for large puncta: where edges have high detHess but centres have low.). }
  \item{ Elbow thresholding to find optimal detHess value that separates puncta from background. Threshold used to mask spots (i.e. detHess > threshold). }
  \item{ Correction of uneven lighting with Laplacian of Gaussian (LoG). }
  \item{ Normalisation of LoG. Floor set to mean LoG values where detHess < threshold. Ceiling set to 99.9\% quantile of non-negative LoG values.}
  \item{ Ensure previously saturated pixels remain saturated (in large regions of saturation, LoG processing returns 0).}
  }

  Effectively returns high-pass features (i.e. puncta).

}
\seealso{
\link{getImageMetrics}
}
\examples{}
