\name{plotQC}
\alias{plotQC}
\title{plotQC}
\description{
Creates default quality control (QC) plots for evaluating MERFISH performance.
}
\usage{
plotQC(
    includeNewBlanks = T,
    synthesisDir = NULL,
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ includeNewBlanks }{ Boolean. Determines if any new blanks added during \code{prepareCodebook} are to be included in evaluation. }
\item{ synthesisDir }{ File path. If specified, files read from this file path. When unspecified, defaults to \code{params$out_dir}. }
\item{ params }{  The parameters environment. }
}
\value{
  Various quality control plots as \code{.png} in the \code{/OUT/QC/} subdirectory of \code{params$parent_out_dir}. Where there are multiple subsets, \code{/QC/} subdirectories are created within subset directories.
}
\details{
  Plots returned include:

  \itemize{
  \item{ \code{FPKMCorrelation}: Correlation of spots per gene to expected expression values (FPKM).}
  \item{ \code{CosineDistanceDistribution}: Distribution of cosine distances per gene (low cosine distance implies more confidence during decoding).}
  \item{ \code{SpotFOVDistribution}: Mean spatial distribution of spots per FOV. Note that this plot was made from all spots called, rather than just those that end up in the final \code{OUT_SPOTCALL_PIXELS.csv}. }
  \item{ \code{CosineSpatialDistribution}: Mean cosine distance per 2D spatial bin (across entire tissue). Note that this plot was made from all spots called, rather than just those that end up in the final \code{OUT_SPOTCALL_PIXELS.csv}. }
  \item{ \code{BitDetectionRateAll}: Boxplots showing if the number of spots per gene depends on bit status (ON vs OFF) for each bit. }
  \item{ \code{BitDetectionRateNonBlank}: Same as \code{BitDetectionRateAll} but only for spots that are non-blanks. }
  \item{ \code{BitwiseErrorRates}: Percentage of "ON to OFF" and "OFF to ON" errors per bit. }
  \item{ \code{FOVMetrics}: Heatmap of FOVs, with each tile showing average metric values. Several metrics reported. }
  \item{ \code{CellMetrics_CellSize}: Mean cell size per 2D spatial bin (across entire tissue). }
  \item{ \code{CellMetrics_nCounts}: Mean cellular counts per 2D spatial bin (across entire tissue). }
  \item{ \code{SpatialPatterns_SpatialPatterns}: Mean spatial pattern score per 2D spatial bin (across entire tissue). }
  \item{ \code{SpatialPatterns_GeneMembership}: Heatmap showing extent of gene's expression explained by a given spatial pattern. }
  }

  Unless otherwise stated above, QC plots are derived from output in the \code{OUT} directory produced by \code{synthesiseData}.

  For more information on \code{SpatialPatterns}, see documentation for \code{getPNNMatrix} and \code{getPNNLatentFactors}. Additionally, the output of \code{getPNNMatrix} is saved as \code{PNNMatrix.csv} in the QC folder.

  If missing certain files -- e.g. if no \code{fpkmFileName} is provided for \code{establishParams} -- some plots will be skipped.
}
\seealso{
\link{ prepareCodebook }
\link{ establishParams }
\link{ getPNNMatrix }
\link{ getPNNLatentFactors }
\link{ synthesiseData }
}
\examples{}
