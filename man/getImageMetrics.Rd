\name{getImageMetrics}
\alias{getImageMetrics}
\title{getImageMetrics}
\description{
Recommended image processing.
}
\usage{
getImageMetrics(
    imList,
    params = get('params', envir = globalenv()),
    imageFunctions = list(
      'ORIG' = imReturn,
      'NORM' = imWinsorIntensities,
      'COARSE' = imLowPass,
      'DECODE' = imForDecode,
      'MASK' = imForMask
    )
)
}
\arguments{
\item{ imList }{ Named list of 2D numeric image matrices to register. Names should be \code{bit_names}, that indicate image and channel. }
\item{ params }{  The parameters environment. }
\item{ imageFunctions }{  Named list of image processing functions. }
}
\value{
  The \code{imMetrics} environment, and updated \code{params} environment (to contain \code{params$imageMetrics}).
}
\details{
  Provides recommended image metrics for downstream decoding and registration purposes. As default, five types of processing are done, and will each return a corresponding named list in \code{imMetrics} according to names provided in \code{imFunctions}: e.g. \code{imMetrics$ORIG}.

  \code{imFunctions} is customisable: users can provide their own functions. As default, will return:

  \itemize{
  \item{ORIG: Unprocessed image.}
  \item{NORM: Image with intensities clipped between \code{params$brightness_min} and \code{params$brightness_max} (if provided, otherwise between image min and max values). }
  \item{COARSE: Low-pass filtered image (optimised to separate foreground from background).}
  \item{DECODE: Decoding-optimised image (preserving high-pass features).}
  \item{MASK: Masks of pixels of interest to bring forwards for decoding.}
  }

  Details on image processing functions can be found in their corresponding documentations.

  Notes on usage:
  \itemize{
  \item{ For image registration, a combination of low-pass (COARSE) and high-pass (DECODE) features are recommended. }
  \item{ For decoding, currently returning continuous values. Users may or may not wish to binarise continuous values to 0/1 to avoid ambiguity. }
    }

}
\seealso{
}
\examples{}
