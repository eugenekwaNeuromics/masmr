\name{spotcall_troubleshootPlots}
\alias{spotcall_troubleshootPlots}
\title{spotcall_troubleshootPlots}
\description{
Creates plots for troubleshooting spotcalling pipeline.
}
\usage{
spotcall_troubleshootPlots(
    spotcalldf,
    imList,
    chosenCoordinate,
    plotWindowRadius = 10,
    decodeMetric = 'DECODE',
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ imList }{ Named list of 2D numeric image matrices. Names should be \code{bit_names}, that indicate image and channel. }
\item{ chosenCoordinate }{ Complex or complex vector. Indicates coordinate central pixel(s) of window to plot around.}
\item{ plotWindowRadius }{ Integer. Number of pixels from central pixel to window boundary. }
\item{ decodeMetric }{ Character. References an object in the \code{imMetrics} environment (if \code{imMetrics} still exists). }
\item{ params }{ The parameters environment. }
}
\value{
  A list of plots. Additionally, adds a colour palette to \code{params$genePalette}, one colour for each gene.
}
\details{
  For each coordinate specified in \code{chosenCoordinate}, a plot is generated showing pixel intensity values across bits, with gene identities (determined by \code{spotcalldf}) overlaid atop the images. If  \code{chosenCoordinate} is not specified, defaults to the first pixel in \code{spotcalldf}.

  A list of images (\code{imList}) can be provided. Otherwise, will look for the object referenced by \code{decodeMetric} in \code{imMetrics} (generated by \code{getImageMetrics}). If \code{imMetrics} does not exist (and \code{imList} not provided), an error message is returned.

  Designed to work with \code{filterDF}.
}
\seealso{
\link{filterDF}
\link{getImageMetrics}
}
\examples{
## Below returns two plots in the plotList: one centred on (x,y) of (200,100), and another at (100,120)
plotList <- spotcall_troubleshootPlots(spotcalldf, chosenCoordinate = c(200+100i, 100+120i))
}
