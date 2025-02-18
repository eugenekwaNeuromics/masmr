\name{saveDF}
\alias{saveDF}
\title{saveDF}
\description{
Save dataframe.
}
\usage{
saveDF(
    spotcalldf,
    params = get('params', envir = globalenv()),
    currentFovName = NULL
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ params }{ The parameters environment. }
\item{ currentFovName }{ Character. Name of current FOV. If not specified, will first attempt to read \code{params$current_fov} for same information. }
}
\value{
  A \code{.csv.gz} file in \code{params$out_dir}. Named in the format of \code{SPOTCALL_{FOV}.csv.gz}.
}
\details{
}
\seealso{
\link{data.table::fwrite}
}
\examples{}
