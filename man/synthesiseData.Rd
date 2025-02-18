\name{synthesiseData}
\alias{synthesiseData}
\title{synthesiseData}
\description{
Synthesises spotcall and cell segmentation data across multiple FOVs.
}
\usage{
synthesiseData(
    referenceFOV,
    stitchChosenColumn = 1,
    nLeadingZeroes = 5,
    cellOverlapFraction = 0.25,
    removeRedundancy = TRUE,
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{  referenceFOV  }{ Character. Choice of FOV that other FOVs are stitched with respect to. If unspecified, automatically chosen from FOV with the most spotcalls. }
\item{  stitchChosenColumn  }{ A single or a vector of integer(s) or character(s. Specifies the column in stitch dataframes to use for stitching. If an integer, the first column after \code{fov} is chosen. If character, the specific stitch column is used. If a vector, assumes multiple columns specified; in which case, \strong{mean} stitch vector calculated. }
\item{ nLeadingZeroes }{ Numeric. Number of leading zeroes to append before cell number. Ideally, select a value bigger than the maximum number of cells per FOV. E.g. if there are at most 9999 cells in any FOV, then a safe \code{nLeadingZeroes} value is >4. }
\item{ cellOverlapFraction }{ Numeric between 0 and 1. Fraction of cell overlap, above which one of the overlapping cells has to be dropped (as opposed to removing overlapping pixels). }
\item{ removeRedundancy }{ Boolean. Determines if redundant cell segmentation and spotcalls are filtered out. See Details. }
\item{ params }{ The parameters environment. }
}
\value{
  Outputs files to the \code{OUT} subdirectory in \code{params$parent_out_dir}.
}
\details{
  Outputs files are:
  \itemize{
   \item{ \code{OUT_GLOBALCOORD.csv}: Updated global coordinates for each FOV. }
   \item{ \code{OUT_SPOTCALL_PIXELS.csv.gz}: Concatenated spotcall dataframes across all FOVs. Dataframe has updated coordinates according to stitch vectors (\code{Xm}, \code{Ym}), and a \code{CELLNAME} column indicating which cell a given spot belongs to (NA if spot does not overlap a cell). }
   \item{ \code{OUT_CELLSEG_PIXELS.csv.gz}: Concatenated cell segmentation masks across all FOVs. Dataframe has updated coordinates according to stitch vectors (\code{Xm}, \code{Ym}), and a \code{CELLNAME} column (which matches the same column in the spotcall dataframe). }
   \item{ \code{OUT_CELLEXPRESSION.csv.gz}: Gene expression profiles summarised by cell. }
   \item{ \code{OUT_CELLS.csv.gz}: Cellular meta data such as centroid coordinates, cell sizes in pixels, and number of counts for that cell. }
  }

  When \code{removeRedundancy} is \code{TRUE}:
  \itemize{
   \item{ Removes redundant spots: Overlapping FOV regions are likely to return more spots, because of spotcalling occuring more than once. Hence, spots belonging to one of the FOVs are ignored. The FOV chosen is the one with fewer spots in the overlapping region.  }
  \item{ Removes redundant cells: If two cells overlap <\code{X}\%, overlapping pixels between cell masks from different FOVs are ignored. If two cells overlap >\code{X}\%, the smaller cell is dropped. If a given cell in the current FOV overlaps with more than two cells significantly (>\code{X}\%), that cell is also dropped. \code{X} is determined from \code{cellOverlapFraction}. }
  }

}
\seealso{
}
\examples{
}
