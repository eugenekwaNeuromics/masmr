\name{eigenMapMatMult2}
\alias{eigenMapMatMult2}
\title{eigenMapMatMult2}
\description{
Fast matrix multiplication in C++.
}
\usage{
eigenMapMatMult2( A, B, n_cores )
}
\arguments{
\item{ A }{ First matrix. }
\item{ B }{ Second matrix. }
\item{ n_cores }{ Integer. Number of threads for parallel computation. }
}
\value{
  Matrix multiplication between matrices \code{A} and \code{B}.
}
\details{
  Depends on \code{Rcpp}, \code{RcppEigen}, and \code{RcppArmadillo}.
}
\seealso{
}
\examples{}
