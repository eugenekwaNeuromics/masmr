\name{readImageList}
\alias{readImageList}
\title{readImageList}
\description{
Reads a list of images for a given FOV.
}
\usage{
readImageList(
    chosenFOV = NULL,
    fileNames = NULL,
    safeLoad = T,
    waitSeconds = 1,
    params = get('params', envir = globalenv()),
    ...
  )
}
\arguments{
\item{ chosenFOV }{ Character or character vector. Chosen field of view (FOV). All images belonging to a given FOV will be loaded. More than one FOV can be specified. }
\item{ fileNames }{ Character vector. If \code{chosenFOV} not provided, alternative input is to provide a list of image file names to load. }
\item{ safeLoad }{ Boolean. If \code{TRUE}, images are loaded channel by channel. Otherwise, no subsetting performed. See Details for more information. }
\item{ waitSeconds }{ Numeric. Number of seconds to wait between image loads when \code{safeLoad} is \code{FALSE}. See Details for more information. }
\item{ params }{ The parameters environment. }
\item{ ... }{  Other arguments passed to \code{readImage}. }
}
\value{
  A named list of images. If more than one FOV provided, a list of lists will be returned.
}
\details{
  As default \code{safeLoad} should be left as \code{TRUE} if the goal is to load all channels for a given FOV. This is to avoid running into the \code{java.lang.OutOfMemoryError: Java heap space} error, seen when too many images are loaded in close succession. [\code{unloadNamespace('RBioFormats')} is implemented in \code{readImage} to similarly avoid this issue.]

  However, if only a single channel is desired across multiple FOVs (as in \code{getAnchorParams}), \code{safeLoad} has to be set to \code{FALSE} to avoid multiple \code{subset} parameters being passed to \link{readImage}. To compensate, a short wait duration is implemented (\code{waitSeconds}) between image loads.
}
\seealso{
\link{readImage}
\link{getAnchorParams}
}
\examples{}
