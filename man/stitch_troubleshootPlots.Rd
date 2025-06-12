\name{stitch_troubleshootPlots}
\alias{stitch_troubleshootPlots}
\title{stitch_troubleshootPlots}
\description{
Creates plots for evaluating stitching.
}
\usage{
stitch_troubleshootPlots(
    imList,
    stitchResults,
    alphaRange = c(0, 0.75),
    shiftColumn = NULL,
    marginToPlot = NULL,
    stitchParams = get('stitchParams', envir = globalenv())
)
}
\arguments{
\item{ imList }{ Named list of numeric image matrices / arrays to register. Names should be \code{bit_names}, that indicate image and channel. }
\item{ stitchResults }{ A dataframe with at least one column containing stitch vector (such as the output of \code{stitchImages}), or a complex vector housing the same information. }
\item{ alphaRange }{ A vector of two numbers between 0 and 1. First number is the transparency (alpha) when pixel intensity is minimum, second for maximum. }
\item{ shiftColumn }{ Integer or string indicating which column in \code{stitchResults} contains the stitch vector. If unspecified, defaults to first column housing complex numbers. Ignored when \code{stitchResults} is not a dataframe. }
\item{ marginToPlot }{ Integer. Number of rows / columns of the reference and neighbour images to plot. When unspecified, automatically derived from how much shift has occurred. }
\item{ stitchParams }{ The stitching parameters environment (generated from \code{readImagesForStitch}). }
}
\value{
  A list of plots, showing the reference FOV in red and the neighbour FOVs in blue.
}
\details{
  Designed to work with \code{stitchImages}.
}
\seealso{
\link{stitchImages}
}
\examples{
## Below plots 100 rows / columns for each neighbour
plotList <- stitch_troubleshootPlots(imList, stitchResults, marginToPlot=100)
}
