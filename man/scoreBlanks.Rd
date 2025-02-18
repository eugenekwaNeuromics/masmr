\name{scoreBlanks}
\alias{scoreBlanks}
\title{scoreBlanks}
\description{
Determines if pixel is to be labelled as a 'blank' according to different criteria.
}
\usage{
scoreBlanks(
    spotcalldf,
    criteria = c(
      'gBlank',
      'gbPossible',
      'consistentLabels'
      ),
    distanceMetrics = c('COS', 'EUC', 'L2E'),
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ criteria }{  Character vector. Determines which tests to perform. See Details below for currently acceptable tests. }
\item{ distanceMetrics }{  Character vector. Specifies which distance metrics to consider for various tests: see \code{decode}. }
\item{ params }{ The parameter environment. }
}
\value{
  An integer vector specifying number of tests that a pixel has failed. E.g. a pixel that has failed no tests (0) is likely \strong{not} a blank.
}
\details{
  Currently acceptable tests:
   \itemize{
   \item{ gBlank: Whether the \code{g} column in the dataframe refers to a blank. Cross-references \code{params$isblank}. }
   \item{ bPossible: Whether the next closest gene label for each distance metric in \code{distanceMetrics} is possible according to the codebook. [Not included as default.] }
    \item{ gbPossible: Similar to \code{bPossible} but only for the distance metric that was used to derive the \code{g} column (typically \code{COSLAB}). }
    \item{ consistentLabels: Determines if the closest labels agree across distance metrics (e.g. if \code{COSLAB} agrees with \code{EUCLAB} and \code{L2ELAB}). Accordingly, more than one distance metric needs to be specified.}
   }

   These tests are not meant to be exhaustive or definitive. Users are strongly recommended to also perform distance-based thresholding and any other relevant downstream filtering. [Filtering performance can be monitored by considering enrichment/depletion of blanks.]
}
\seealso{
\link{decode}
}
\examples{
blankTestFails <- scoreBlanks(spotcalldf)
spotcalldf$BLANK <- blankTestFails>0
}
