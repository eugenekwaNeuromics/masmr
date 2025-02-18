\name{stitchImages}
\alias{stitchImages}
\title{stitchImages}
\description{
Stitches adjacent FOVs together.
}
\usage{
stitchImages(
    imList,
    registerTo = 1,
    stitchParams = get('stitchParams', envir=globalenv())
)
}
\arguments{
\item{ imList }{ Named list of 2D numeric image matrices to register. Names should be \code{fov_names}, that indicate which FOV an image comes from. }
\item{  registerTo  }{ Numeric (>0). Indicates which image in \code{imList} houses the reference image (of the current FOV). If \code{imList} from \code{readImagesForStitch}, then should reference the first image (i.e. \code{registerTo = 1}). }
\item{ stitchParams }{ The stitching parameters environment (generated from \code{readImagesForStitch}). }
}
\details{
  Function relies on maximising cross-correlation (see \code{crossCorrelate2D}) to stitch adjacent FOVs to the current FOV. The offset in coordinates for neighbouring FOVs (listed under column \code{fov}) is provided under column (\code{shift_pixel}) as a complex number (i.e. offset in pixels for x-coordinate is real; and y-coordinate is imaginary).

  Users may wish to attempt stitching of \strong{both} MERFISH and DAPI images. An example of how to do so is provided below.
}
\value{
  A dataframe with \code{fov} and \code{shift_pixel} columns.
}
\seealso{
\link{crossCorrelate2D}
\link{readImagesForStitch}
}
\examples{
## Assume establishParams() has been run
## Assume readImageMetaData() has been run for both MERFISH and DAPI images - creating respective GLOBALCOORD and META_RESOLUTION files
## Assume current FOV is listed under params$current_fov

## Stitch MERFISH images
imFunc <- function(im){ imNormalise(im)^0.1 } #Brighten image by raising to power 0.1
imList <- readImagesForStitch( subDirectory = 'IM', imFunction = imFunc )
merfish_stitch <- stitchImages( imList )
rm( stitchParams )

## Stitch DAPI images
imFunc <- function(im){ imNormalise(im)^1 } #DAPI typically bright, so additional brightening may not be needed
imList <- readImagesForStitch( subDirectory = 'DAPI', imFunction = imFunc )
dapi_stitch <- stitchImages( imList )
rm( stitchParams )

colnames( merfish_stitch ) <- c('fov', 'merfish_shift_pixel')
colnames( dapi_stitch )    <- c('fov', 'dapi_shift_pixel'   )
combined_stitch <- data.frame(
  merfish_stitch,
  'dapi_shift_pixel' = dapi_stitch[ match(merfish_stitch$fov, dapi_stitch$fov), 'dapi_shift_pixel']
)
head(combined_stitch)
}
