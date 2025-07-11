\name{registerImages}
\alias{registerImages}
\title{registerImages}
\description{
Registers images by maximising weighted cross-correlation.
}
\usage{
registerImages(
    imList,
    registerTo = 1,
    chosenZslice = NULL,
    params = get('params', envir = globalenv()),
    currentFovName = NULL,
    maxAcceptableShiftDistance = Inf,
    gaussianWeightedAlignment = TRUE,
    gaussianWeightAmplitude = NULL,
    gaussianWeightSigmaX = NULL,
    gaussianWeightSigmaY = NULL,
    tolerableBestAlignmentPercentage = 50,
    standardiseShifts = FALSE,
    saveShifts = TRUE,
    returnTroubleShootPlots = FALSE,
    ...
)
}
\arguments{
\item{ imList }{ Named list of 2D numeric image matrices to register. Names should be \code{bit_names}, that indicate image and channel. }
\item{ registerTo }{ Integer. Specifies position of image in \code{imList} to be taken as reference image. As default, first image taken as reference. }
\item{ chosenZslice }{ Integer. Specifies which z slice to perform registration with respect to (only applied when an image stack with more than 1 z slice is detected). If unspecified, the default behaviour is to perform a maximum intensity projection to collapse 3D array into a 2D matrix. }
\item{ params }{  The parameters environment. }
\item{ currentFovName }{ Character. Name of current FOV (used for writing files if \code{saveShifts} is \code{TRUE}). If not specified, will first attempt to read \code{params$current_fov} for same information. }
\item{ maxAcceptableShiftDistance }{  Numeric. Maximum tolerable distance an image can be moved relative to the reference image. }
\item{ gaussianWeightedAlignment }{  Boolean. Determines if a Gaussian weight is applied to cross-correlation calculations. }
\item{ gaussianWeightAmplitude }{  Numeric. Maximum amplitude of weight. When unspecified, set to 1. }
\item{ gaussianWeightSigmaX }{  Numeric. Standard deviation in x direction. When unspecified, set to 20\% of x dimension of image. }
\item{ gaussianWeightSigmaY }{  Numeric. Standard deviation in y direction. When unspecified, set to 20\% of y dimension of image. }
\item{ tolerableBestAlignmentPercentage }{ Numeric between 0 and 100. If shift exceeds \code{maxAcceptableShiftDistance}, highest shift within \code{maxAcceptableShiftDistance} returned instead, provided cross-correlation of new shift is at least \code{tolerableBestAlignmentPercentage}\% of highest cross-correlation value.}
\item{standardiseShifts}{ Boolean. If \code{TRUE}, will attempt to standardise shifts (see Details below). }
\item{saveShifts}{ Boolean. If \code{TRUE}, will save shift vector in a \code{SHIFTS_XXX.csv.gz} file. If \code{standardiseShifts} too, unstandardised values saved in corresponding \code{RAWSHIFTS_XXX.csv.gz}. \code{XXX} indicates  \code{currentFovName}.  }
\item{ returnTroubleShootPlots }{ Boolean. Determines if troubleshooting plots are to be returned (see \code{register_troubleshootPlots}). }
\item{ ... }{ Other parameters passed to \code{register_troubleshootPlots}. }
}
\details{
  Images registered by maximising cross-correlation in Fourier transformed space (see \code{crossCorrelate2D}. A Gaussian filter can be applied to prioritise shifts nearer to 0 pixels moved in x and y directions (i.e. \code{ 0 + 0i } when expressed as a complex number).

  If \code{standardiseShifts}, will assume that main sources of image shifting comes from (1) chromatic aberration across channels, and (2) microscope movement across image acquisition (assumed to be otherwise negligible for the same image, across channels). Thus, mean shift calculated across channels and across images, and summed accordingly to return a standardised shift vector.

  When an image stack is detected, \code{maxIntensityProject} is first performed to collapse the Z slices into a single slice. To bypass this behaviour, registration can instead be performed with respect to a single chosen Z slice.
}
\value{
  A complex vector \code{params$shifts} that contains the calculated shifts, and a numeric vector \code{params$intersecting_window} that specifies the bounding box (min x, max x, min y, max y) of intersecting coordinate values post-transformation.

  If \code{saveShifts}, \code{SHIFTS_*.csv.gz} files containing the shifts for each image is generated. The reference image is also saved in same folder as \code{REGIM_*.png}.
}
\seealso{
\link{crossCorrelate2D}
\link{register_troubleshootPlots}
}
\examples{}
