\name{getAnchorParams}
\alias{getAnchorParams}
\title{getAnchorParams}
\description{
Obtains metrics across a subset of "anchor" images.
}
\usage{
getAnchorParams(
    anchors,
    out_dir,
    returnTroubleShootPlots = FALSE,
    imageFunctions = list(
      'brightness_min' = imsBrightnessMin,
      'brightness_max' = imsBrightnessMax
    ),
    summaryFunctions = list(
      'brightness_min' = conservativeMean,
      'brightness_max' = conservativeMean
    ),
    anchorMode = c('sample', 'grid'),
    nSamples = 10,
    anchorSpacing = 3,
    anchorOffset = 2,
    freshAnchors = FALSE,
    params = get('params', envir = globalenv())
    ...
)
}
\arguments{
\item{ anchors }{ Character vector. Subset of \code{params$fov_names} to load. If unspecified, automatically determined: see Details below. }
\item{ out_dir }{ Character. Output directory. If unspecified, read from \code{params$out_dir}. }
\item{ returnTroubleShootPlots }{ Boolean. Determines if troubleshooting plots are to be returned. }
\item{ imageFunctions }{ Named list of functions, with names also shared by \code{summaryFunctions}. Functions chosen must return a single numeric per image. Currently defaults to \code{imsBrightnessMin} and \code{imsBrightnessMax}. }
\item{ summaryFunctions }{ Named list of functions, with names also shared by \code{imageFunctions}. Functions chosen to summarise one or more values. Currently defaults to \code{conservativeMean}. }
\item{ anchorMode }{ Character. Determines if anchor FOVs are to be selected from regularly spaced FOVs (\code{'grid'}) or randomly sampled (\code{'sample', default}). }
\item{ nSamples }{ Integer. Number of FOVs to sample. Used when \code{anchors} is unspecified, and \code{anchorMode} is \code{sample}. }
\item{ anchorSpacing }{ Integer. Number of FOVs between anchor FOVs in same row. Used when \code{anchors} is unspecified and \code{anchorMode} is \code{grid}. }
\item{ anchorOffset }{ Integer. Number of FOVs to offset by for even numbered rows. Used when \code{anchors} is unspecified and \code{anchorMode} is \code{grid}. }
\item{ freshAnchors }{ Boolean. Determines if values are to be re-obtained. Functionally equivalent to setting \code{params$resumeMode} to \code{FALSE}.}
\item{ params }{ The parameters environment. }
\item{ ... }{  Other arguments passed to \code{readImageList}, \code{imageFunctions}, and \code{summaryFunctions}. }
}
\value{
  New objects in the parameters environment, and new \code{.csv} files. If \code{returnTroubleShootPlots}, a \code{troubleshootPlots} environment is also created.
}
\details{
  "Anchors" are a subset of FOVs for which global parameters are derived: for instance, background fluorescence intensity across images (\code{brightness_min}), or fluorescence intensity cap above which saturation is assumed (\code{brightness_max}).

  When unspecified, anchors are selected from evenly spaced FOVs interior to the edge (\code{anchorMode = 'grid'}), or sampled from all existing FOVs ( \code{anchorMode = 'sample'} ).

  Function returns several outputs:

  \itemize{
  \item{ Objects in the \code{params} environment, named with respect to \code{imageFunctions}: e.g. \code{params$brightness_min}. }
  \item{ \code{.csv} files in \code{params$out_dir} with the corresponding metrics per bit: e.g. \code{BRIGHTNESS_MIN.csv} for \code{params$brightness_min}.}
  \item{ A \code{ggplot2} plot (\code{params$anchors_plot}) indicating which FOVs were selected as anchors (red). Can be visualised with \code{print}. }
  \item{ Other records in \code{params}: \code{anchors}, \code{anchors_metrics}. }
  }

  Functions in \code{imageFunctions} and \code{summaryFunctions} are user-selected and customisable. For more details on default functions, see their relevant documentations.

  For each \code{imageFunction}, we also report the proportion of pixels below the obtained threshold per bit per anchor FOV. These can be found in the \code{ANCHOR_QUANTILEQC.csv}.

  If \code{returnTroubleShootPlots}, a \code{troubleshootPlots} environment created containing histograms of pixel intensities per FOV per bit, and thresholds (e.g. max or min brightnesses) as a vertical red dashed line. Plots are be specified in \code{troubleshootPlots} according to metric then bit: e.g. \code{troubleshootPlots[['brightness_max']][[12]]} specifies intensity histograms across all anchor FOVs for bit 12 with the max brightness threshold indicated in red.

}
\seealso{
\link{readImageList}
\link{imsBrightnessMin}
\link{imsBrightnessMax}
\link{conservativeMean}
}
\examples{
getAnchorParams()
print( params$anchors_plot )
}
