## This script contains different image processing functions

require(imager)
require(EBImage)

imReturn_raw <- function(
    im,
    currentIteration,
    ...
){
  if(!missing(currentIteration)){
    message(currentIteration, appendLF = F)
  }
  return(im)
}

imReturn <- function(
    im,
    currentIteration,
    ...
){
  if(!missing(currentIteration)){
    message(currentIteration, appendLF = F)
  }
  im <- imNormalise(im)
  return(im)
}

##

imNormalise <- function(
    im,
    floorVal,
    ceilVal,
    floorQuant,
    ceilQuant,
    ...
){
  if(missing(floorVal)){ floorVal = NULL }
  if(missing(ceilVal)){ ceilVal = NULL }
  if(missing(floorQuant)){ floorQuant = NULL }
  if(missing(ceilQuant)){ ceilQuant = NULL }

  if(is.null(floorVal)){
    floorVal = min(im)
  }

  if(!is.null(floorQuant)){
    if( (floorQuant<0) | (floorQuant>1) ){
      stop('Invalid quantile value!')
    }
    floorVal <- quantile(im[im>=floorVal], floorQuant)
  }

  if(is.null(ceilVal)){
    ceilVal = max(im)
  }

  if(!is.null(ceilQuant)){
    if( (ceilQuant<0) | (ceilQuant>1) ){
      stop('Invalid quantile value!')
    }
    ceilVal <- quantile(im[im<=ceilVal], ceilQuant)
  }

  if(floorVal >= ceilVal){
    stop('floorVal cannot be bigger than ceilVal!')
  }

  norm <- (im - floorVal) / (ceilVal - floorVal)
  norm[norm<0] <- 0
  norm[norm>1] <- 1
  return(norm)
}

##

imWinsorIntensities <- function(
    im,
    currentIteration,
    ...
){
  if(!missing(currentIteration)){
    message(currentIteration, appendLF = F)
  }
  norm <- imNormalise(im, ...)
  return(norm)
}

##

imAutoBrighten <- function(
    im,
    floorVal,
    ceilVal,
    ...
){
  if(missing(floorVal)){ floorVal = 0.25 }
  if(missing(ceilVal)){ ceilVal = 0.75 }

  if(floorVal <=0){
    stop('Invalid floor value!')
  }
  if(ceilVal >= 1){
    stop('Invalid ceiling value!')
  }
  if(ceilVal <= floorVal){
    stop('Ceiling value cannot be smaller than floor value!')
  }

  # Maximise ceil^x - floor^x
  # => Solve ceil^x * log(ceil) - floor^x * log(floor) = 0
  x = log(log(ceilVal) / log(floorVal)) / log(floorVal / ceilVal)
  return( im^x )
}

##

imLowPass <- function(
    im,
    smallBlur,
    currentIteration,
    ...
){
  if(!missing(currentIteration)){
    message(currentIteration, appendLF = F)
  }
  if(missing(smallBlur)){ smallBlur = 1 }

  background <- imNormalise(im, ...)==0
  imblur <- array( imager::isoblur(suppressWarnings(imager::as.cimg(im)), smallBlur), dim=dim(im) )
  lowthresh <- median( imblur[background] )
  hithresh <- median( imblur[!background & (imblur > lowthresh)] )
  norm <- imAutoBrighten( imblur, floorVal = lowthresh, ceilVal = hithresh )
  norm <- imNormalise(norm)
  return(norm)
}

##

imHessianDeterminant <- function(
    im,
    smallBlur,
    ...
){
  if(missing(smallBlur)){ smallBlur = 1 }
  hess <- imager::imhessian(imager::isoblur(suppressWarnings(imager::as.cimg(im)), smallBlur))
  dethess <- hess$xx * hess$yy - (hess$xy^2)
  dethess <- array(dethess, dim=dim(im))
  return(dethess)
}


##

imLaplacianOfGaussian <- function(
    im,
    smallBlur,
    bigBlur,
    ...
){
  if(missing(smallBlur)){
    smallBlur = 1
  }
  if(missing(bigBlur)){
    bigBlur = 5
  }
  norm <- imager::isoblur(suppressWarnings(imager::as.cimg(im)), smallBlur) -
    imager::isoblur(suppressWarnings(imager::as.cimg(im)), bigBlur)
  norm <- array(norm, dim=dim(im))
  return(norm)
}

##

imForDecode <- function(
    im,
    smallBlur,
    currentIteration,
    ...){

  if(!missing(currentIteration)){
    message(currentIteration, appendLF = F)
  }
  if(missing(smallBlur)){ smallBlur = 1 }

  # Cap between low and high values
  norm <- imNormalise(im, ...)

  # Calculate determinant of hessian
  dethess <- imHessianDeterminant(norm, ...)
  dethess <- imNormalise(dethess)
  dethessb <- array(imager::isoblur(suppressWarnings(imager::as.cimg(dethess)), smallBlur),
                       dim = dim(im))
  vals <- dethessb[dethessb > quantile(dethessb, 0.95)]
  dethessbFloor <- findThreshold(vals, method = thresh_Elbow, ...)

  # LoG to correct for uneven lighting
  log <- imLaplacianOfGaussian(norm, ...)

  # Normalise
  logFloor <- mean( log[dethessb < dethessbFloor] )
  logCeil <- quantile(log[log>0], 0.999)
  final <- imNormalise(log, floorVal = logFloor, ceilVal = logCeil)

  # If image is saturated, LoG may not be helpful
  final[norm==1] <- 1
  return(final)
}

##

imForMask <- function(
    im,
    smallBlur,
    minBlobSize,
    currentIteration,
    ...){

  if(!missing(currentIteration)){
    message(currentIteration, appendLF = F)
  }
  if(missing(smallBlur)){ smallBlur = 1 }
  if(missing(minBlobSize)){ minBlobSize = 9 }

  # Cap between low and high values
  norm <- imNormalise(im, ...)
  lowpass <- imLowPass(im, ...)

  # Calculate determinant of hessian
  dethess <- imHessianDeterminant(norm, ...)
  dethess <- imNormalise(dethess)
  dethessb <- array(imager::isoblur(suppressWarnings(imager::as.cimg(dethess)), smallBlur),
                       dim=dim(im))
  vals <- dethessb[dethessb > quantile(dethessb, 0.95)]
  dethessbFloor <- findThreshold(vals, method = thresh_Elbow, ...)

  # LoG to correct for uneven lighting
  LoG <- imLaplacianOfGaussian(norm, ...)
  LoGthresh <- findThreshold(
    LoG[LoG>0], method = thresh_Elbow, ...
  )
  LoG <- LoG > LoGthresh

  vals <- imNormalise( (norm * lowpass) * (LoG * norm)  )
  lab <- ( LoG + (dethessb>dethessbFloor) + (im==1) ) > 0
  thresh <- findThreshold(
    vals, method = thresh_Quantile, labels = lab,
    ...)

  # Mask values
  mask <- vals > thresh
  labmask <- EBImage::bwlabel(mask)
  toosmall <- which(tabulate(labmask) < minBlobSize)
  mask[labmask %in% toosmall] <- F

  final <- matrix(as.integer(mask), nrow=nrow(mask), ncol=ncol(mask))
  return(final)
}





