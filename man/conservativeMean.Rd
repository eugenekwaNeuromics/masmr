\name{conservativeMean}
\alias{conservativeMean}
\title{conservativeMean}
\description{
Calculates mean for a subset of values.
}
\usage{
conservativeMean(
    values,
    nSDS = 0.9,
    ...
)
}
\arguments{
\item{ values }{ Numeric vector. }
\item{ nSDs }{ Numeric (>0). Number of standard deviations away from mean, within which acceptable \code{values} fall.  }
\item{ ... }{ Other functions passed from \code{getAnchorParams}.  }
}
\value{
  A numeric: conservative estimation of the mean of \code{values}.
}
\details{
}
\seealso{
\link{getAnchorParams}
}
\examples{
}
