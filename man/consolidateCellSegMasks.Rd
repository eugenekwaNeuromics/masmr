\name{consolidateCellSegMasks}
\alias{consolidateCellSegMasks}
\title{consolidateCellSegMasks}
\description{
Consolidates a list of masks into a single mask matrix.
}
\usage{
consolidateCellSegMasks(
    maskList,
    minFlags = 1,
    cellSeg = get('cellSeg', envir = globalenv()),
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{  maskList  }{ List of 2D mask matrices. If not a list, \code{maskList} is returned as is. }
\item{  minFlags  }{ Numeric. Minimum number of matrices that flag a given pixel, for that pixel to be saved as a mask. As default, a pixel must be masked at least 1 time for it to remain masked. Stricter criteria possible by raising the \code{minFlags} value, up to and including the total number of images in \code{maskList}. }
\item{ cellSeg }{ The cell segmentation environment. }
\item{ params }{ The parameters environment. }
}
\value{
  A single summarised mask matrix.
}
\details{
}
\seealso{
}
\examples{
}
