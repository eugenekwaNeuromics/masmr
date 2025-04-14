\name{readImagesForStitch}
\alias{readImagesForStitch}
\title{readImagesForStitch}
\description{
Reads FOVs around a chosen FOV into memory to prepare for stitching.
}
\usage{
readImagesForStitch(
    currentFOVName,
    subDirectory = 'IM',
    loadProcessedImages = T,
    registerTo = 1,
    zSlice = 1,
    imageFunction = function(im){ return(im) },
    params = get('params', envir=globalenv()),
    ...
)
}
\arguments{
\item{ currentFovName }{ Character. Name of current FOV. If not specified, will first attempt to read \code{params$current_fov} for same information. }
\item{ subDirectory }{ Character. Name of subdirectory containing required \code{GLOBALCOORD} and \code{META_RESOLUTION} files, generated from \code{readImageMetaData}. This directory should be in \code{params$parent_out_dir}. Typically \code{IM} or \code{DAPI}. }
\item{ loadProcessedImages }{ Boolean. Determines if processed images -- \code{.png} image files prefixed with \code{REGIM_} in \code{subDirectory} --- are loaded. Processed images are saved optionally as part of \code{registerImages} and \code{saveCellSegMasks} functions. If \code{TRUE}, attempts to find processed images first. If processed images do not exist, automatically switched to \code{FALSE} and raw images are loaded. }
\item{ registerTo }{ Integer. Specifies bit position of image to be taken as reference image (for MERFISH images) -- should be the same as \code{registerImages}. As default, first image taken as reference. }
\item{ zSlice }{ Integer. Specifies z slice of image to be taken as reference image, if more than one z slice exists in loaded image. As default, first z slice taken as reference. }
\item{ imageFunction }{ Function. Image processing to be performed on each loaded image. }
\item{ params }{ The parameters environment. }
\code{ ... }{ Other parameters passed to \code{readImage}, \code{readImageList}, or \code{imFunction}. }
}
\details{
  Function is able to load processed images that are generated as a normal part of image processing, provided these have been saved earlier (see \code{registerImages}, \code{saveCellSegMasks}). It is thus the default recommendation that stitching be performed \strong{after} spotcalling and cell segmentation.

  However, more advanced users can save time by running spotcalling and stitching scripts in parallel. This can be achieved by setting \code{loadProcessedImages} to \code{FALSE}; and ensuring that at least \code{GLOBALCOORD} and \code{META_RESOLUTION} files have been generated (as part of running \code{readImageMetaData} function early in the spotcalling pipeline). This function then loads raw images instead of processed images.

  If processed images are unavailable, the raw images loaded can still undergo image processing using the \code{imFunction} field. For \code{imFunction}, the default is to return the image as is. However, some image processing (e.g. brightening, band-pass filtering) is generally recommended before attempting to stitch. See \code{stitchImages} for details.

  Images loaded are typically those of the reference bit (for \code{subDirectory = 'IM'}) or the DAPI channel (for \code{subDirectory = 'DAPI'}). Images for the current FOV and spatially adjacent FOVs are loaded, with the current FOV's image being first in the returned list.
}
\value{
  A named list of images to perform stitching on (as input to \code{stitchImages}), and a new environment \code{stitchParams} in the global environment that houses parameters for image stitching.
}
\seealso{
\link{readImageMetaData}
\link{registerImages}
\link{saveCellSegMasks}
\link{readImage}
\link{readImageList}
\link{stitchImages}
}
\examples{}
