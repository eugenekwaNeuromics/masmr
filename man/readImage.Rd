\name{readImage}
\alias{readImage}
\title{readImage}
\description{
Reads image data using \code{RBioFormats}, or directly for \code{.dax} files.
}
\usage{
readImage(
    fileName,

    imageDimensions = NULL,

    nchannels = NULL,
    nrows = NULL,
    ncols = NULL,
    nzs = NULL,

    channelIndex = NULL,
    zIndex = NULL,

    endian = 'little',
    nbytes = 2,
    signed = T,
    normalise = T,
    ... )
}
\arguments{
\item{ fileName }{ Character. Image file name. }
\item{ imageDimensions }{ Integer vector. Indicates the dimensions of the image being loaded. If unspecified, defaults to \code{c(Number of Channels, Number of rows, Number of Columns, Number of Z slices)} (see Details).}
\item{ nchannels }{ Integer. Number of channels the image should have. }
\item{ nrows }{ Integer. Number of rows the image should have. }
\item{ ncols }{ Integer. Number of columns the image should have.}
\item{ nzs }{ Integer. Number of Z slices the image should have.}
\item{ channelIndex }{ Integer. When specified, only the indicated channel is loaded.}
\item{ zIndex }{ Integer. When specified, only the indicated Z slice is loaded.}
\item{ endian }{ Character. Endianness. Used for binary file reading. Accepts \code{little} (default), \code{big}, or their initial letters. See \code{readBin}. }
\item{ nbytes }{ Integer. Number of bytes (8 bits) encoding each element. Used for binary file reading. See \code{readBin}. }
\item{ signed }{ Boolean. Determines if file should return signed or unsigned integers. Used for binary file reading. See \code{readBin}.}
\item{ normalise }{ Boolean. Determines if image intensity should be normalised to between 0 to 1. }
\item{ ... }{  Other arguments passed to \code{read.image} from \code{RBioFormats}. }
}
\value{
  A 2-4D array containing image data, or a list of 2-3D arrays containing image data.
}
\details{
  If dimensions of the image is not specified, the dimensions are assumed to be: [Number of Channels, Number of rows, Number of Columns, Number of Z slices]. This information is obtained from the \code{nchannels}, \code{nrows}, \code{ncols}, and \code{nzs} fields. Providing \code{imageDimensions} is optional for \code{ome.tif} formats but mandatory for \code{.dax} binary formats (otherwise a single 1D vector will be returned).

  \code{channelIndex} is meant for a specific use case in \code{getAnchorParams}, where specific channels are loaded to save memory. The equivalent when using \code{RBioFormats::read.image} is to include the parameter \code{subset =  list( c = channelIndex )}. Note: because \code{readImage} function loads images one channel and one Z-slice at a time (to avoid Java heap space error), users \strong{cannot} simultaneously provide a \code{channelIndex} and \code{ subset = list( c = channelIndex ) } due to clash in parameter names. The same logic applies to \code{zIndex}, which replaces \code{ subset = list(z = zIndex) } for \code{RBioFormats::read.image}.
}
\seealso{
\link{readBin}
\link{RBioFormats::read.image}
\link{readImageList}
\link{getAnchorParams}
}
\examples{}
