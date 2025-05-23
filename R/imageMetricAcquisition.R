## Calculates image metrics

require(imager)

getImageMetrics <- function(
    imList,
    params = get('params', envir = globalenv()),
    imageFunctions = list(
      'ORIG' = imReturn,
      'NORM' = imWinsorIntensities,
      'COARSE' = imLowPass,
      'DECODE' = imForDecode,
      'MASK' = imForMask
    )
){

  if(!exists('imMetrics', envir = globalenv())){
    imMetrics <<- new.env()
  }

  message('\nPreparing imMetrics environment...')
  for(i in 1:length(imageFunctions)){
    message(paste0('\nLooping ', names(imageFunctions)[i], ' function...'))

    if(exists( names(imageFunctions)[i], envir = imMetrics)){
      warning( paste0(names(imageFunctions)[i]), ' will be overridden' )
    }
    imFunc <- imageFunctions[[i]]

    processedList <- lapply(
      1:length(imList),
      function(j){
        imFunc(
          imList[[j]],
          floorVal = params$brightness_min[j],
          ceilVal = params$brightness_max[j],
          currentIteration = paste0(j, ' of ', length(imList), '...')
          )
      })
    names(processedList) <- names(imList)
    imMetrics[[names(imageFunctions)[i]]] <<- processedList
  }

  params$imageMetrics <<- names(imageFunctions)

}


consolidateImageMetrics <- function(
    maskList,
    bitFloor,
    bitCeil,
    imMetrics = get('imMetrics', envir = globalenv()),
    params = get('params', envir = globalenv())
    ){

  shifts <- params$shifts
  window <- params$intersecting_window
  if(is.null(shifts) | is.null(window)){
    stop('It appears image registration has not been performed yet!')
  }

  maskFilt = T
  if(missing(maskList)){
    maskList <- imMetrics$MASK
    if(is.null(maskList)){
      warning('maskList not provided and imMetrics$MASK does not exist...Assuming want to save all information...')
      maskFilt = F
    }else{
      #warning('Using imMetrics$MASK to filter pixels...')
    }
  }
  if(missing(bitFloor)){ bitFloor = 1 }
  if(missing(bitCeil)){ bitCeil = params$nbits * 2 }

  acceptable_idx <- 0:prod(as.numeric(as.numeric(params$resolutions$xydimensions_pixels)))
  if(maskFilt){
    message('\nGetting acceptable pixel locations...')
    acceptable_idx <- lapply(1:length(maskList), function(idx){
      message(paste0(idx, ' of ', length(maskList), '...'), appendLF = F)
      imx <- matrix(as.numeric(maskList[[idx]]), nrow=nrow(maskList[[idx]]), ncol=ncol(maskList[[idx]]))
      df <- as.data.frame(imager::as.cimg(imx))
      df$WX <- df$x + Re(shifts[idx])
      df$WY <- df$y + Im(shifts[idx])
      df <- df[(df$WX >= window[1]) & (df$WX < window[2]) & (df$WY >= window[3]) & (df$WY < window[4]), ]
      df$IDX <- as.integer(df$WX + (window[2] * (df$WY - 1) ))
      df <- df[df$value>0,]
      return(df$IDX)
    })
    quant <- tabulate(unlist(acceptable_idx))
    acceptable_idx <- intersect( which( (quant >= bitFloor) ), which( (quant <= bitCeil) ) )
    acceptable_idx <- unique(acceptable_idx)
    message(paste0(
      round( 100 * length(acceptable_idx) / (window[2] * window[4]), digits=2), '% of pixels to keep...'
    ))
  }

  message('\nConsolidating metrics per pixel...')
  available_metrics <- ls(envir=imMetrics)
  raw_images <- imMetrics[[available_metrics[1]]]
  spotcalldf <- list()
  for(idx in 1:length(shifts)){
    message(paste0(idx, ' of ', length(shifts), '...'), appendLF = F)
    base <- as.data.frame(imager::as.cimg(raw_images[[idx]]))
    if(sum( (as.vector(raw_images[[idx]]) - base$value)^2 ) != 0 ){
      print("ERROR: WRONG INDEXING!")
      break
    }
    base$value <- NULL
    base$WX <- base$x + Re(shifts[idx])
    base$WY <- base$y + Im(shifts[idx])
    base$IDX <- as.numeric(base$WX + (window[2] * (base$WY-1)))
    base <- base[,c("WX", "WY", 'IDX')]

    bool <- (base$WX >= window[1]) & (base$WX < window[2]) & (base$WY >= window[3]) & (base$WY < window[4])
    df <- base[bool,]
    for(mi in 1:length(available_metrics)){
      val <- imMetrics[[ available_metrics[mi] ]]
      val <- val[[idx]]
      df$new <- as.vector( val[bool] )
      colnames(df)[colnames(df)=='new'] <- available_metrics[mi]
    }
    df <- df[,]

    df <- df[match(acceptable_idx, df$IDX),]
    if(idx==1){
      coorddf <- df[,which(colnames(df) %in% colnames(base))]
    }
    df <- df[,-which(colnames(df) %in% colnames(coorddf))]
    colnames(df) <- paste0(colnames(df), '_',  sprintf(paste0('%0', nchar(length(shifts)), 'd'), idx))
    spotcalldf[[idx]] <- df
  }

  nrow_check <- nrow(coorddf)
  if(!all(sapply(spotcalldf, function(x) nrow(x)==nrow_check))){
    stop("nrows not consistent across images!")
  }
  spotcalldf <- data.frame(do.call(cbind, spotcalldf), check.names = F)
  spotcalldf <- data.frame(cbind(coorddf, spotcalldf), check.names = F)

  message('\nRemoving pixels where we have incomplete information...')
  filtout <- rowSums(is.na(spotcalldf)) > 0
  spotcalldf <- spotcalldf[!filtout,]

  return(spotcalldf)
}
