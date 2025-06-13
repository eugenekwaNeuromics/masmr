\name{getPNNLatentFactors}
\alias{getPNNLatentFactors}
\title{getPNNLatentFactors}
\description{
Performs non-negative matrix factorization (NMF) on PNNMatrix (output from \code{getPNNMatrix}).
}
\usage{
getPNNLatentFactors(

  PNNMatrix,

  nFactors = 4,

  mergeSimilarFactors = TRUE,
  correlationThreshold = 0.1,

  seed = 12345,
  verbose = T,
  ...
)
}
\arguments{
\item{ PNNMatrix }{ Matrix or dataframe. E.g. output from \code{getPNNMatrix}.  }
\item{ nFactors }{ Integer. Number of factors (aka rank).  }

\item{ mergeSimilarFactors }{ Boolean (default \code{TRUE}).  Determines if similar factors are to be merged. }
\item{ correlationThreshold }{ Numeric. Maxmimum allowable Pearson correlation of any two factors. If exceeded, factors are merged. }

\item{ seed }{ Integer. Seed for pseudo-random number generation. }
\item{ verbose }{ Boolean (default \code{TRUE}). If \code{TRUE}, returns messages from each function. Otherwise, messages not returned (but warnings and stop messages still reported). }

\item{ ... }{ Other parameters passed to the \code{RcppML::nmf} function. }
}
\value{
  A list containing the eventual number of factors after merging (\code{n_factors}), the feature matrix (W) with values min-max scaled per factor (\code{point_scores}), and the coefficient matrix (H) with values min-max scaled per factor (\code{coefficients}).
}
\details{
  Designed to work with the output of \code{getPNNMatrix}, but generaliseable use to any non-negative matrix.
}
\seealso{
\link{ RcppML::nmf }
}
\examples{

}
