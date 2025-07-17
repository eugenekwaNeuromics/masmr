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

  if( ( (lowthresh>0) & (hithresh<1) & (lowthresh<hithresh) & !is.na(lowthresh) & !is.na(hithresh) ) ){
    norm <- imAutoBrighten( imblur, floorVal = lowthresh, ceilVal = hithresh )
  }else{
    warning( paste0(
      'Skipping imAutoBrighten because of invalid floorVal (', lowthresh,
      ') and/or ceilVal (', hithresh, ')...'
    ))
    norm <- imblur
  }
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

imTVDenoise <- function(
    im,
    denoisingWeight = 0.001,
    stopThreshold = 1,
    maxIterations = 200,
    returnIntermediateImages = F,
    verbose = NULL
){
  if(is.null(verbose)){
    verbose <- get('params', globalenv())$verbose
  }
  if(is.null(verbose)){
    verbose <- TRUE
  }

  ## We will call images u, as in the original ROF paper
  if( length(dim(im))>3 | is.list(im) | length(dim(im)) < 2 ){
    stop('imTVDenoise can only accept a 2D matrix or 3D array')
  }
  if( length(dim(im))==3 ){

    resultArray <- array(0, dim=dim(im))
    resultList <- list()
    for(zi in dim(im)[3]){
      resultx <- imTVDenoise(
        im[,,zi],
        denoisingWeight = denoisingWeight,
        stopThreshold = stopThreshold,
        maxIterations = maxIterations,
        returnIntermediateImages = returnIntermediateImages,
        verbose = verbose
      )
      if(returnIntermediateImages){
        resultList[[zi]] <- resultx
      }else{
        resultArray[,,zi] <- resultx
      }
    }
    if(returnIntermediateImages){
      return( resultList )
    }else{
      return( resultArray )
    }

  }

  if( length(dim(im))==2 ){

    g <- im #The observed image is 'g' and the cleaned image is 'u'
    px <- py <- array(0, dim=dim(g))
    i = 1
    converged = F
    if(returnIntermediateImages){
      imList <- list()
      imList[[i]] <- g
    }

    ## Now let's define some functions as per Chambolle
    ## The nabla function: getting gradients
    nabla <- function( u ){
      gradx = rbind( u[-1,], u[nrow(u),] ) - u #when i=N, gx = 0
      grady = cbind( u[,-1], u[,ncol(u)] ) - u
      return( list('x' = gradx, 'y' = grady) )
    }

    ## The modulus function: getting norms
    modulus <- function( xylist ){
      sqrt( xylist[['x']]^2 + xylist[['y']]^2 )
    }

    ## The div function: getting divergence
    div <- function( xylist ){
      # p is a list / vector
      x <- xylist[['x']]
      y <- xylist[['y']]

      trim <- rbind(x[-nrow(x),], 0)
      shift <- rbind(0, x[-nrow(x),])
      divpx <- trim - shift

      trim <- cbind(y[,-ncol(y)], 0)
      shift <- cbind(0, y[,-ncol(y)])
      divpy <- trim - shift

      divp = divpx + divpy

      return(divp)
    }

    ## Now we can start looping
    TV_init = TV_previous = sum(modulus(nabla(g)))
    ubest <- u <- g
    while( (i < maxIterations) & !converged ){

      if(verbose){ message(paste0('\nEpoch: ', i, '...'), appendLF = F) }

      ## We are updating p every cycle (p starts off as 0)
      divp <- div( list('x' = px, 'y' = py) )

      ## Now calculate the next p(n+1)
      ## Image fidelity is our lambda term
      tau = 1/4 # Recommended in Chambolle to ensure convergence
      forUpdating <- nabla(divp - g / denoisingWeight)
      denominator <- 1 + tau * modulus(forUpdating)

      ## Update px and py
      px <- (px + tau * forUpdating$x) / denominator
      py <- (py + tau * forUpdating$y) / denominator

      ## Estimate nonlinear project plk(g)
      ## Which converges to lambda * div( pn )
      plkg <- div( list('x' = px, 'y' = py) ) * denoisingWeight
      u = g - plkg
      TV = sum( modulus(nabla(u)) )

      ## Check for convergence
      if(verbose){
        message(paste0(
          'TV: ', TV, ' ( ', round( 100 * TV / TV_previous, digits=3), '% previous TV )...'
        ), appendLF = F)
      }
      if( (TV / TV_previous) >= stopThreshold ){
        converged = T
        if(verbose){message('\nConverged!', appendLF = T)}
      }
      TV_previous = TV

      if(!converged){
        ubest <- u
        if( returnIntermediateImages ){
          imList[[i+1]] <- ubest
        }
      }

      i = i + 1
    }

    if(returnIntermediateImages){
      return(imList)
    }
    return(ubest)
  }
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
  lowpass <- imLowPass(norm, ...)

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

  vals <- imNormalise( lowpass * imNormalise(LoG, floorVal = LoGthresh, ceilVal = max(LoG)) )
  vals[( (im>=1) + (norm>=1) )>0] <- 1
  lab <- ( (LoG>LoGthresh) + (dethessb>dethessbFloor) + (im==1) ) > 0
  thresh <- findThreshold(
    vals, method = thresh_Quantile, labels = lab,
    ...)

  # Mask values
  mask <- vals > thresh
  labmask <- EBImage::bwlabel(mask)
  toosmall <- which(tabulate(labmask) < minBlobSize)
  mask[labmask %in% toosmall] <- F

  return(mask)
}
###
