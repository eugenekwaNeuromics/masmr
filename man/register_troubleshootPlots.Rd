\name{register_troubleshootPlots}
\alias{register_troubleshootPlots}
\title{register_troubleshootPlots}
\description{
Creates plots for evaluating registration.
}
\usage{
register_troubleshootPlots(
    imList,
    chosenCoordinate,
    alphaRange = c(0, 0.75),
    plotWindowRadius = 50,
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ imList }{ Named list of numeric image matrices / arrays to register. Names should be \code{bit_names}, that indicate image and channel. }
\item{ chosenCoordinate }{ Complex or complex vector. Indicates coordinate central pixel(s) of window to plot around.}
\item{ alphaRange }{ A vector of two numbers. First number is the transparency (alpha) when pixel intensity is minimum, second for maximum. }
\item{ plotWindowRadius }{ Integer. Number of pixels from central pixel to window boundary. }
\item{ params }{ The parameters environment. }
}
\value{
  A list of plots, showing the reference FOV in blue and the registered FOVs in red.
}
\details{
  If \code{chosenCoordinate} is not specified, defaults to brightest pixel in reference FOV.

  Designed to work with \code{registerImages}.
}
\seealso{
\link{registerImages}
}
\examples{
## Below returns two plots in the plotList: one centred on (x,y) of (200,100), and another at (100,120)
plotList <- register_troubleshootPlots(imList, chosenCoordinate = c(200+100i, 100+120i))
}
