\name{readImageMetaData}
\alias{readImageMetaData}
\title{readImageMetaData}
\description{
Attempts to parse image metadata.
}
\usage{
readImageMetaData(
    dir = 'im_dir',
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ dir }{ Character. Directory to read metadata from. Can be \code{im_dir} or \code{dapi_dir}, or some object in \code{params} that directs to a folder. }
\item{ params }{ The parameters environment. }
}
\value{
  Updates \code{params$out_dir}: appends an additional subfolder named after \code{dir}.

  Outputs the \code{GLOBALCOORD.csv} and \code{META_RESOLUTION.txt} files there, and updates the parameters environment \code{params} to contain resolution and global coordinate information.
}
\details{
  Function first checks if \code{resumeMode} is TRUE in \code{params}. If yes, and \code{GLOBALCOORD.csv} and \code{META_RESOLUTION.txt} exist, function will only update the \code{params} object.

  \code{GLOBALCOORD.csv} provides the x, y, z coordinates (in microns) of each image, location of image files, and assigns a name to every bit.

  \code{META_RESOLUTION.txt} contains resolution information: number of microns per pixel, dimensions of image in pixels and microns, order of fluorescence channels, and FOV names sorted by acquisition order.

  If unable to parse metadata, template \code{GLOBALCOORD.csv} and \code{META_RESOLUTION.txt} files are generated, for users to manually edit.

}
\seealso{}
\examples{}
