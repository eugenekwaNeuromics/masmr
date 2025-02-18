\name{thresh_Otsu}
\alias{thresh_Otsu}
\title{thresh_Otsu}
\description{
Otsu thresholding.
}
\usage{
thresh_Otsu(
    values,
    thresholdsToTry = NULL,
    nbreaks = NULL,
    ...
)
}
\arguments{
\item{ values }{ Numeric vector. Values to calculate thresholds for. }
\item{ thresholdsToTry }{ Numeric vector. Specifies which threshold values to test. If unspecified, automatically selected from range of \code{values}. }
\item{ nbreaks }{ Integer. Number of histogram breaks: see \link{hist}.}
\item{ ... }{ Other parameters passed from \code{findThreshold} to \link{hist}. }
}
\value{
  A numeric: the Otsu threshold.
}
\details{
}
\seealso{
\link{findThreshold}
}
\examples{}
