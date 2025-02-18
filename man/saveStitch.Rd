\name{saveStitch}
\alias{saveStitch}
\title{saveStitch}
\description{
Saves stitch dataframe.
}
\usage{
saveStitch(
    stitchDF,
    stitchParams = get('stitchParams', envir=globalenv())
)
}
\arguments{
\item{ stitchDF }{ Dataframe returned from \code{stitchImages}, containing a \code{fov} column and at least one other column with shift vectors as complex numbers. }
\item{ stitchParams }{ The stitching parameters environment (generated from \code{readImagesForStitch}). }
}
\value{
  Saves dataframe in \code{stitchParams$out_dir} as a \code{.csv}.
}
\details{
}
\seealso{
\link{stitchImages}
}
\examples{
}
