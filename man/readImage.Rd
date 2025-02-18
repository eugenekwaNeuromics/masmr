\name{readImage}
\alias{readImage}
\title{readImage}
\description{
Reads image data using \code{RBioFormats}, or directly for \code{.dax} files.
}
\usage{
readImage(
    fileName,
    nrows = NULL,
    ncols = NULL,
    endian = c('little', 'big'),
    nbytes = 2,
    signed = T,
    normalise = T,
    ... )
}
\arguments{
\item{ fileName }{ Character. Image file name. }
\item{ nrows }{ Integer. Number of rows the image should have. }
\item{ ncols }{ Integer. Number of columns the image should have.}
\item{ endian }{ Character. Endianness. Used for binary file reading. Accepts \code{little} (default), \code{big}, or their initial letters. See \code{readBin}. }
\item{ nbytes }{ Integer. Number of bytes (8 bits) encoding each element. Used for binary file reading. See \code{readBin}. }
\item{ signed }{ Boolean. Determines if file should return signed or unsigned integers. Used for binary file reading. See \code{readBin}.}
\item{ normalise }{ Boolean. Determines if image intensity should be normalised to between 0 to 1. }
\item{ ... }{  Other arguments passed to \code{read.image} from \code{RBioFormats}. }
}
\value{
  A 2D matrix or 3D array containing image data.
}
\details{
  If \code{nrows} or \code{ncols} unspecified, will attempt to find appropriate values by assuming a square image.
}
\seealso{
\link{readBin}
\link{RBioFormats::read.image}
\link{readImageList}
}
\examples{}
