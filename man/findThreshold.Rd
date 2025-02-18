\name{findThreshold}
\alias{findThreshold}
\title{findThreshold}
\description{
Automatically calculates threshold given a numeric vector.
}
\usage{
findThreshold(
    values,
    thresholdMethod,
    ...
)
}
\arguments{
\item{ values }{ Numeric vector. Values to calculate thresholds for. }
\item{ thresholdMethod }{ Function to perform thresholding with. }
\item{ ... }{ Other parameters passed to \code{thresholdMethod}. }
}
\value{
  A numeric threshold.
}
\details{
  Several options exist for \code{thresholdMethod}:

  \itemize{
  \item{\code{thresh_Otsu} (default)}
  \item{\code{thresh_Elbow}}
  \item{\code{thresh_Quantile}}
  }

  See respective functions for details.

}
\seealso{
\link{thresh_Otsu}
\link{thresh_Elbow}
\link{thresh_Quantile}
}
\examples{}
