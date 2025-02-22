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
    returnTroubleShootPlots = FALSE,
    troubleShootCoordinates = NULL,
    logBase = exp(1),
    params = get('params', envir = globalenv()),
    ...
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ filterOut }{ Boolean vector. Same length as nrows of \code{spotcalldf}. If \code{TRUE}, row in \code{spotcalldf} is filtered out.  }
\item{ returnTroubleShootPlots }{ Boolean. Determines if troubleshooting plots are to be returned. }
\item{ troubleShootCoordinates }{ Complex or complex vector. Indicates which central coordinates around which users would like to plot (see \code{spotcall_troubleshootPlots}). If \code{NULL}, defaults to the pixel with the lowest \code{COS} among those to be kept.}
\item{ logBase }{ Numeric (>0). Base for log-transformation: see \code{log}. }
\item{ params }{ The parameters environment. }
\item{ ... }{ Other parameters passed to \code{spotcall_troubleshootPlots}. }
}
\value{
  A filtered dataframe. If \code{returnTroubleShootPlots}, a \code{troubleshootPlots} environment is also created.
}
\details{
  Functionally identical to \code{spotcalldf[!filterOut,]}. If \code{filterOut} is instead an expression to derive a boolean vector, e.g. \code{spotcalldf$HD > 1}, that is also acceptable.

  If \code{params$n_probes_file} and/or \code{params$fpkm_file} detected, will attempt to parse. If successful, correlations of N puncta per gene to respective metrics will be reported. Correlation is calculated after log-transformation with 1 pseudocount.

  If \code{returnTroubleShootPlots}, a \code{troubleshootPlots} environment is created containing a an image before and after filtering. See \code{spotcall_troubleshootPlots} for additional parameters users can specify.
}
\seealso{
\link{log}
\link{spotcall_troubleshootPlots}
}
\examples{}
