\name{filterDF}
\alias{filterDF}
\title{filterDF}
\description{
Filters dataframe while reporting some QC metrics.
}
\usage{
filterDF(
    spotcalldf,
    filterOut,
    logBase = exp(1),
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ filterOut }{ Boolean vector. Same length as nrows of \code{spotcalldf}. If \code{TRUE}, row in \code{spotcalldf} is filtered out.  }
\item{ logBase }{ Numeric (>0). Base for log-transformation: see \code{log}. }
\item{ params }{ The parameters environment. }
}
\value{
  A filtered dataframe.
}
\details{
  Functionally identical to \code{spotcalldf[!filterOut,]}.

  If \code{params$n_probes_file} and/or \code{params$fpkm_file} detected, will attempt to parse. If successful, correlations of N puncta per gene to respective metrics will be reported. Correlation is calculated after log-transformation with 1 pseudocount.
}
\seealso{
\link{log}
}
\examples{}
