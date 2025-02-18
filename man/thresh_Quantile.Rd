\name{thresh_Quantile}
\alias{thresh_Quantile}
\title{thresh_Quantile}
\description{
Finds the mean between a quantile value in group A, and a quantile value in group B.
}
\usage{
thresh_Quantile(
    values,
    labels,
    quantileFalse = NULL,
    quantileTrue = NULL,
    ...
)
}
\arguments{
\item{ values }{ Numeric vector. Values to calculate thresholds for. }
\item{ labels }{ Boolean (or integer) vector. Determines which elements in \code{values} belong to which group. }
\item{ quantileFalse }{ Numeric. Quantile of \code{values} when \code{labels} is \code{FALSE}. If unspecified, defaults to median. }
\item{ quantileTrue }{ Numeric. Quantile of \code{values} when \code{labels} is \code{TRUE}. If unspecified, defaults to median. }
\item{ ... }{ Other parameters passed from \code{findThreshold} to \link{quantile}. }
}
\value{
  A numeric: the mean of quantiles.
}
\details{
}
\seealso{
\link{findThreshold}
}
\examples{}
