\name{scoreDeltaF}
\alias{scoreDeltaF}
\title{scoreDeltaF}
\description{
Calculates delta F from image metrics.
}
\usage{
scoreDeltaF(
    spotcalldf,
    prefixes,
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ prefixes }{  Character vector. Column prefixes to subset by: see \code{subsetDF}. If unspecified, defaults to \code{params$imageMetrics}. }
\item{ params }{ The parameter environment. }

}
\value{
  Updated dataframe with additional delta F metric columns.
}
\details{
  For each prefix indiciating an image metric (e.g. \code{DECODE}), the following are calculated:

  \itemize{
  \item{ F1: Mean image metric value when bit is ON (1 in codebook).}
  \item{ F0: Mean image metric value when bit is OFF (0 in codebook).}
  \item{ DF: The difference between F1 and F0.}
  \item{ DFF0: The difference between F1 and F0, divided by F0. When F0 is 0, the maximum non-infinite DFF0 value is returned instead. }
  \item{ F1VAR: Variance of image metric when bit is ON (1 in codebook). }
  \item{ F0VAR: Variance of image metric when bit is OFF (0 in codebook). }
  \item{ F1CV: Coefficient of variation of image metric when bit is ON (1 in codebook). }
  \item{ F0CV: Coefficient of variation of image metric when bit is OFF (0 in codebook). }
  }

  New column names are formatted with delta F as prefix, and image metric name as suffix: e.g. \code{F1_DECODE}.

  \code{spotcalldf} must have a column titled \code{g}, which indicates rownames in \code{params$ordered_codebook}.
}
\seealso{
\link{subsetDF}
}
\examples{
spotcalldf <- scoreDeltaF(spotcalldf)
}
