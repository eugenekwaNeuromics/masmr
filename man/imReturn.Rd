\name{imReturn}
\alias{imReturn}
\title{imReturn}
\description{
Returns image without processing.
}
\usage{
imReturn(
    im,
    currentIteration
)
}
\arguments{
\item{ im }{ 2D numeric image matrix. }
\item{ currentIteration }{ Character. For keeping track of progress. }
}
\value{
  Normalised image matrix with values between 0 and 1.
}
\details{
  Other than image intensitites being normalised between 0 (minimum intensity) and 1 (maximum intensity), no other processing is performed.
}
\seealso{
\link{getImageMetrics}
}
\examples{}
