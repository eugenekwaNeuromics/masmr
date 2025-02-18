\name{imWinsorIntensities}
\alias{imWinsorIntensities}
\title{imWinsorIntensities}
\description{
Returns image with winsorised intensities.
}
\usage{
imWinsorIntensities(
    im,
    currentIteration,
    ...
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ currentIteration }{ Character. For keeping track of progress. }
\item{ ... }{ Other parameters passed to \code{imNormalise}. }
}
\value{
  Processed image matrix.
}
\details{
  Functionally identical to \code{imNormalise}, but with the additional \code{currentIteration} parameter.
}
\seealso{
\link{getImageMetrics}
}
\examples{}
