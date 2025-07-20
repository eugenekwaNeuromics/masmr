## Generic functions for handling 2D data

## Get coordinates for an image.
getRasterCoords <- function( raster, realY = FALSE ){
  rowLength = nrow(raster)
  colLength = ncol(raster)
  squareLength = pmax(rowLength, colLength)
  coordX <- matrix( rep( seq(0, squareLength-1, length.out=squareLength) ), nrow=squareLength, ncol=squareLength)
  coordY <- t(coordX)
  if(realY){
    coords <- coordY + coordX * sqrt(as.complex(-1))
  }else{
    coords <- coordX + coordY * sqrt(as.complex(-1))
  }
  coords <- coords[c(1:rowLength), c(1:colLength)]
  return(coords)
}

## Calculate cross correlation.
crossCorrelate2D <- function(referenceImageMatrix, queryImageMatrix, pad = FALSE){
  x = referenceImageMatrix
  h = queryImageMatrix

  h = h - mean(h)
  x = x - mean(x)

  hdim = dim(as.matrix(h))
  xdim = dim(as.matrix(x))
  cdim = hdim + xdim - 1 # Cross-correlation dimension

  if(pad){
    ## When images are not the same dimension
    hpad = xpad = matrix(0, nrow=cdim[1], ncol=cdim[2])
    xpad[1:xdim[1], 1:xdim[2]] = x
    hpad[1:hdim[1], 1:hdim[2]] = h[hdim[1]:1, hdim[2]:1]
  }else{
    if(any(hdim != xdim)){
      stop('Expecting query and reference images to have the same dimensions!')
    }
    xpad = x
    hpad = h[hdim[1]:1, hdim[2]:1]
  }
  fftx = fft(xpad)
  ffth = fft(hpad) # The more accurate version is Conj(fft(unflipped h)), but this works better for some reason
  res = fft(fftx * ffth, inverse = TRUE)
  return(Re(res))
}

## Identify local peaks in 2D data.
get2DPeaksSimple <- function(imRaster, searchspace = 'queen', verbose = F){
  valid_ss <- c('queen', 'bishop', 'rook', 'q', 'b', 'r')
  searchspace <- tolower(searchspace)
  if(!all(searchspace %in% valid_ss)){
    stop("Invalid searchspace specified!")
  }
  if(length(searchspace) > 1){
    warning('Only first searchspace considered')
    searchspace <- searchspace[1]
  }
  if(searchspace %in% c('queen', 'q')){
    searchspace = 'queen'
  }
  if(searchspace %in% c('rook', 'r')){
    searchspace = 'rook'
  }
  if(searchspace %in% c('bishop', 'b')){
    searchspace = 'bishop'
  }

  ispeak <- matrix(T, nrow=nrow(imRaster), ncol=ncol(imRaster))

  if(searchspace != 'bishop'){
    ## Check top
    if(verbose){cat('Top...')}
    subim <- imRaster[-1,]
    subim <- rbind(subim, rep(0, ncol(subim) ))
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)

    ## Check left
    if(verbose){cat('Left...')}
    subim <- imRaster[,-1]
    subim <- cbind(subim, rep(0, nrow(subim) ))
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)

    ## Check right
    if(verbose){cat('Right...')}
    subim <- imRaster[,-ncol(imRaster)]
    subim <- cbind(rep(0, nrow(subim)), subim)
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)

    ## Check bottom
    if(verbose){cat('Bottom...')}
    subim <- imRaster[-nrow(imRaster),]
    subim <- rbind(rep(0, ncol(subim) ), subim)
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)
  }

  if(searchspace != 'rook'){
    ## Check top left
    if(verbose){cat('Top left...')}
    subim <- imRaster[-1,-1]
    subim <- rbind(subim, rep(0, ncol(subim) ))
    subim <- cbind(subim, rep(0, nrow(subim) ))
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)

    ## Check top right
    if(verbose){cat('Top right...')}
    subim <- imRaster[-1,-ncol(imRaster)]
    subim <- rbind(subim, rep(0, ncol(subim) ))
    subim <- cbind(rep(0, nrow(subim)), subim )
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)

    ## Check bottom left
    if(verbose){cat('Bottom left...')}
    subim <- imRaster[-nrow(imRaster),-1]
    subim <- rbind(rep(0, ncol(subim) ), subim)
    subim <- cbind(subim, rep(0, nrow(subim) ))
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)

    ## Check bottom right
    if(verbose){cat('Bottom right...')}
    subim <- imRaster[-nrow(imRaster),-ncol(imRaster)]
    subim <- rbind(rep(0, ncol(subim) ), subim)
    subim <- cbind(rep(0, nrow(subim)), subim )
    isgreater <- imRaster - subim
    ispeak <- ispeak & (isgreater > 0)
  }

  ispeak <- matrix(as.integer(ispeak), nrow=nrow(ispeak), ncol=ncol(ispeak))
  return(ispeak)
}
