\name{thresh_Elbow}
\alias{thresh_Elbow}
\title{thresh_Elbow}
\description{
Finds the elbow in a histogram.
}
\usage{
thresh_Elbow(
    values,
    rightElbow = T,
    nbreaks = NULL,
    span = 0.5,
    ...
)
}
\arguments{
\item{ values }{ Numeric vector. Values to calculate thresholds for. }
\item{ rightElbow }{ Boolean. Determines if elbow should be searched to the left or right of the highest peak. }
\item{ nbreaks }{ Integer. Number of histogram breaks: see \link{loess}.}
\item{ span }{ Numeric. Span parameter for LOESS regression: see \link{loess}. }
\item{ ... }{ Other parameters passed from \code{findThreshold} to \link{loess} and \link{hist}. }
}
\value{
  A numeric: the elbow of a histogram.
}
\details{
}
\seealso{
\link{findThreshold}
}
\examples{}
