\name{saveCellSegMasks}
\alias{saveCellSegMasks}
\title{saveCellSegMasks}
\description{
Saves cell segmentation masks as a file.
}
\usage{
saveCellSegMasks(
    masks,
    im = NULL,
    currentFOVName = NULL,
    cellSeg = get('cellSeg', envir = globalenv()),
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ masks }{ The 2D integer mask matrix. Saved as a dataframe (background 0 values filtered out). File labelled \code{CELLSEG_{Cell segmentation model}_{currentFOVName}.csv.gz} in \code{params$out_dir}. }
\item{ im }{ A 2D numeric image matrix to save. Can be used for stitching downstream. File labelled as a \code{REGIM_{currentFOVName}.png} in \code{params$out_dir}. If unspecified, will attempt to read from \code{params$cell_seg_image}: see \code{runCellSegModel}. }
\item{ currentFovName }{ Character. Name of current FOV. If unspecified, will first attempt to read from \code{params$current_fov}. }
\item{ cellSeg }{ The cell segmentation environment. }
\item{ params }{ The parameters environment. }
}
\value{
  Files in \code{params$out_dir}: consistently, a \code{.csv.gz} containing coordination information per cell mask. Also, a \code{.png}, if \code{im} is specified.
}
\details{
  Columns in the \code{CELLSEG_XXX.csv.gz}:
  \itemize{
  \item{Xm, Ym: Global X and Y coordinates in microns.}
  \item{WX, WY: X and Y coordinates in pixels, for the current FOV.}
  \item{fov: Current FOV. }
  \item{cell: Unique cell label for the given FOV.}
  }
  When multiple FOVs are consolidated later, a combination of \code{cell} and \code{fov} returns a unique identifier for any given cell mask in any FOV.
}
\seealso{
\link{runCellSegModel}
}
\examples{
}
