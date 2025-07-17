\name{imTVDenoise}
\alias{imTVDenoise}
\title{imTVDenoise}
\description{
Performs total variation (TV) denoising of images.
}
\usage{
imTVDenoise(
    im,
    denoisingWeight = 0.001,
    stopThreshold = 1,
    maxIterations = 200,
    returnIntermediateImages = FALSE,
    verbose = NULL
)
}
\arguments{
\item{ im }{ A 2D/3D numeric image array. If provided a 3D numeric array, TV denoising is performed on each Z slice separately. }
\item{ denoisingWeight }{ Numeric (>0). Higher values reduce image fidelity (i.e. more denoising). }
\item{ stopThreshold }{ Numeric (>0). Algorithm converges when the ratio of previous to current TV exceeds or equals to this number. }
\item{ maxIterations }{ Integer (>0). Maximum number of iterations. }
\item{ returnIntermediateImages }{ Boolean (default FALSE). If true, a list of images for each iteration of denoising prior to convergence is returned. }
\item{ verbose }{ Boolean (default NULL). When unspecified, defaults to \code{params$verbose}. }
}
\value{
  A denoised image matrix, or a list of denoised images.
}
\details{
  Implementation here is of the Chambolle (2004) algorithm [1]. The lambda parameter in [1] is specified by the \code{denoisingWeight} parameter. Do note that the lambda parameter here is the reciprocal of the lambda parameter defined in later texts [2]; we are adhering to the definitions from the original paper.

  Convergence occurs when TV no longer decreases significantly with new iterations. The default threshold is that the new TV must be <100\% of the previous iteration's TV (indicated by 1 for \code{stopThreshold}).

  [1] Chambolle, A. (2004). An algorithm for total variation minimization and applications. \emph{Journal of Mathematical imaging and vision, 20(1)}, 89-97.
  [2] Duran, J., Coll, B., & Sbert, C. (2013). Chambolle's projection algorithm for total variation denoising. \emph{Image processing on Line, 2013}, 311-331.
}
\seealso{
}
\examples{}
