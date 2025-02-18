\name{bitsFromDecode}
\alias{bitsFromDecode}
\title{bitsFromDecode}
\description{
Binarises/thresholds continuous metrics per bit into ON (1) and OFF (0).
}
\usage{
bitsFromDecode(
    spotcalldf,
    trainIndex,
    metricToMaximise = c('f1', 'fb', 'precision', 'recall'),
    fBeta = 2,
    nThresholdTests = 100,
    saveThresholds = T,
    decodeMetric = 'DECODE',
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ trainIndex }{ Boolean vector. Indicate which rows in \code{spotcalldf} to calculate thresholds from. If unspecified, defaults to all rows. }
\code{ metricToMaximise }{ Character. Indicates which metric to maximise in finding thresholds per bit. Currently accepts \code{precision}, \code{recall}, \code{f1}, and \code{fb}. See Details.}
\item{ fBeta }{ Numeric (>0). Beta value for \code{fb} calculation: indicates how much recall is weighted more than precision. Higher values (e.g. 2) weighs recall higher than precision and minimises false negatives; lower values (e.g. 0.5) weights recall lower and minimises false positives. A value of 1 returns the same as \code{f1}. }
\item{ nThresholdTests }{ Integer (>0). Number of thresholds to examine. Thresholds are evenly spaced between minimum and maximum values for a given bit. }
\item{ saveThresholds }{ Boolean. Determines if found thresholds per bit are saved as \code{THRESHBITCALL_{metricToMaximise}_{decodeMetric}_{Current FOV}.csv} in \code{params$out_dir}.}
\item{ decodeMetric }{ Character. A previously calculated image metric that decoding will be performed against (checked against \code{params$imageMetrics}). Should be the same as the metric specified in \code{decode}. }
\item{ params }{  The parameters environment. }
}
\value{
  An updated dataframe with additional \code{BIT_{Bit position}} columns, where 0 indicates OFF (subthreshold values), and 1 indicates ON (above threshold values). Also, a \code{HD} column in the dataframe: showing Hamming distance from expected gene (as indicated in the \code{g} column).
  If \code{saveThresholds}, a \code{.csv} containing threshold values is saved in \code{params$out_dir}.
}
\details{
  Function first finds per-bit threshold that maximises "accuracy" between metric > threshold, and expected codebook values for that bit. Expected codebook value is dependent on the \code{g} column in the dataframe (hence, \code{decode()} needs to have been run first).

  "Accuracy" here refers to a user-specified metric: either Precision, Recall, F1 or FB (F Beta) scores. See https://en.wikipedia.org/wiki/F-score for definitions.

  On usage: in the theoretical pipeline below, this function is called in step 4.

  \enumerate{
  \item{ Image metrics are consolidated per pixel. }
  \item{ Calculation of distance between image metric vector and codebook vector per pixel. Accomplished with \code{decode}. }
  \item{ Preliminary distance-based filtering to remove obvious erroneous calls. Accomplished with \code{filterDF} and \code{findThreshold}. }
  \item{ Hard thresholding and Hamming distance calculation - to serve as basis for further downstream filtering. Accomplished with this function. }
  \item{ Further filtering as needed. }
  }

  Hamming distance is saved as a new column \code{HD} in the dataframe, and can be the basis for further filtering.

  Downstream 0-to-1 and 1-to-0 error calculation rates should rely on hard thresholded values.
}
\seealso{
\link{ decode }
}
\examples{}
