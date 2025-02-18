\name{runCellSegModel}
\alias{runCellSegModel}
\title{runCellSegModel}
\description{
Runs cell segmentation models.
}
\usage{
runCellSegModel(
    im,
    masksOnly = T,
    recordLastImage = T,
    cellSeg = get('cellSeg', envir = globalenv()),
    params = get('params', envir = globalenv()),
    ...
)
}
\arguments{
\item{ im }{ A 2D numeric image matrix, or a list of image matrices. If list is named, names are preserved. }
\item{  masksOnly  }{ Boolean. Determines if only the mask matrix is returned. }
\item{  recordLastImage  }{ Boolean. Determines if last processed image is to be saved as \code{params$cell_seg_image} (which can be later converted into an image for downstream stitching: see \code{saveCellSegMasks}). }
\item{ cellSeg }{ The cell segmentation environment. }
\item{ params }{ The parameters environment. }
\item{ ... }{ Other parameters passed to respective evaluation functions for each type of cell segmentation model: see \code{eval} for Cellpose, and \code{predict_instances} for Stardist. }
}
\value{
  A 2D integer matrix, or a list of matrices. Each matrix contains masks of cells / nuclei.
}
\details{
}
\seealso{
\link{establishCellSegModel}
\link{saveCellSegMasks}
}
\examples{
}
