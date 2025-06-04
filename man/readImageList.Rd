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
    params = get('params', envir = globalenv()),
    ...
  )
}
\arguments{
\item{ chosenFOV }{ Character or character vector. Chosen field of view (FOV). All images belonging to a given FOV will be loaded. More than one FOV can be specified. }
\item{ fileNames }{ Character vector. If \code{chosenFOV} not provided, alternative input is to provide a list of image file names to load. }
\item{ params }{ The parameters environment. }
\item{ ... }{  Other arguments passed to \code{readImage}. }
}
\value{
  A named list of images. If more than one FOV provided, a list of lists will be returned.
}
\details{
}
\seealso{
\link{readImage}
\link{getAnchorParams}
}
\examples{}
