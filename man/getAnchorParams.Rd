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
    params = get('params', envir = globalenv())
    ...
)
}
\arguments{
\item{ anchors }{ Character vector. Subset of \code{params$fov_names} to load. If unspecified, automatically determined: see Details below. }
\item{ out_dir }{ Character. Output directory. If unspecified, read from \code{params$out_dir}. }
\item{ imageFunctions }{ Named list of functions, with names also shared by \code{summaryFunctions}. Functions chosen must return a single numeric per image. Currently defaults to \code{imsBrightnessMin} and \code{imsBrightnessMax}. }
\item{ summaryFunctions }{ Named list of functions, with names also shared by \code{imageFunctions}. Functions chosen to summarise one or more values. Currently defaults to \code{conservativeMean}. }
\item{ anchorMode }{ Character. Determines if anchor FOVs are to be selected from regularly spaced FOVs (\code{'grid'}) or randomly sampled (\code{'sample', default}). }
\item{ nSamples }{ Integer. Number of FOVs to sample. Used when \code{anchors} is unspecified, and \code{anchorMode} is \code{sample}. }
\item{ anchorSpacing }{ Integer. Number of FOVs between anchor FOVs in same row. Used when \code{anchors} is unspecified and \code{anchorMode} is \code{grid}. }
\item{ anchorOffset }{ Integer. Number of FOVs to offset by for even numbered rows. Used when \code{anchors} is unspecified and \code{anchorMode} is \code{grid}. }
\item{ params }{ The parameters environment. }
\item{ ... }{  Other arguments passed to \code{readImageList}, \code{imageFunctions}, and \code{summaryFunctions}. }
}
\value{
  New objects in the parameters environment, and new \code{.csv} files.
}
\details{
  "Anchors" are a subset of FOVs for which global parameters are derived: for instance, background fluorescence intensity across images (\code{brightness_min}), or fluorescence intensity cap above which saturation is assumed (\code{brightness_max}).

  When unspecified, anchors are selected from evenly spaced FOVs interior to the edge (\code{anchorMode = 'grid'}, default), or sampled from all existing FOVs ( \code{anchorMode = 'sample'} ).

  Function returns several outputs:

  \itemize{
  \item{ Objects in the \code{params} environment, named with respect to \code{imageFunctions}: e.g. \code{params$brightness_min}. }
  \item{ \code{.csv} files in \code{params$out_dir} with the corresponding metrics per bit: e.g. \code{BRIGHTNESS_MIN.csv} for \code{params$brightness_min}.}
  \item{ A \code{ggplot2} plot (\code{params$anchors_plot}) indicating which FOVs were selected as anchors (red). Can be visualised with \code{print}. }
  \item{ Other records in \code{params}: \code{anchors}, \code{anchors_metrics}. }
  }

  Functions in \code{imageFunctions} and \code{summaryFunctions} are user-selected and customisable. For more details on default functions, see their relevant documentations.

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
