## Function for getting max and min brightness

require(imager)
require(ggplot2)
require(data.table)

getAnchorParams <- function(
    anchors,
    out_dir,
    returnTroubleShootPlots = FALSE,
    imageFunctions = list(
      'brightness_min' = imsBrightnessMin_MIP,
      'brightness_max' = imsBrightnessMax_MIP
    ),
    summaryFunctions = list(
      'brightness_min' = conservativeMean,
      'brightness_max' = conservativeMean
    ),
    anchorMode = 'sample',
    nSamples = 10,
    anchorSpacing = 3,
    anchorOffset = 2,
    freshAnchors = FALSE,
    params = get('params', envir = globalenv()),
    ...
){
  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose
  }else{
    verbose = T
  }

  if( exists('seed', envir = params) ){
    if( is.numeric( params$seed ) ){
      set.seed(params$seed)
    }
  }else{
    set.seed(12345)
  }
  anchorMode <- tolower(anchorMode)[1]
  valid_modes <- c('sample', 'grid')
  if(!(anchorMode %in% valid_modes)){
    stop('Invalid anchorMode specified!')
  }

  ## Get necessary parameters
  if( !all(names(imageFunctions) %in% names(summaryFunctions)) ){
    stop('imageFunctions and summaryFunctions must have the same names!')
  }
  if(!exists('global_coords', envir=params)){
    stop('This function requires readImageMetaData() to have been run first!')
  }
  global_coords <- params$global_coords
  if(!exists('ordered_codebook', envir=params)){
    f <- list.files(params$parent_out_dir,
                    pattern='ORDEREDCODEBOOK.csv.gz',
                    full.names = T, recursive = T)[1]
    if( length(f)!=1 ){
      stop('This function requires prepareCodebook() to have been run first!')
    }else{
      if(verbose){ message('Loading ORDEREDCODEBOOK.csv.gz...') }
    }
    codebook <- data.table::fread(f, data.table=FALSE)
    if(colnames(codebook)[1] == 'V1'){
      rownames(codebook) <- codebook$V1
      codebook$V1 <- NULL
    }
    params$ordered_codebook <<- codebook
    params$capture_order <<- colnames(codebook)
    params$isblank <<- grepl('^blank-', tolower(rownames(codebook)))
    if(verbose){ message(paste0('Codebook has ', sum(params$isblank), ' blanks and ', sum(!params$isblank), ' genes...' )) }
  }
  codebook <- params$ordered_codebook
  if(missing(out_dir)){
    out_dir <- params$out_dir
    if(is.null(out_dir)){
      stop('out_dir not specified and params$out_dir missing!')
    }
  }


  ## If anchors not specified, create
  # freshAnchors = F
  if(missing(anchors)){
    anchors <- params$anchors

    ## If fresh anchors, then re-run
    # freshAnchors = T
    if(!params$resumeMode | freshAnchors ){
      anchors <- NULL
    }

    if(is.null(anchors)){
      if(verbose){ message('\nGetting new anchors...') }
      tileCx <- global_coords[!duplicated(global_coords$fov),]
      tileCx$X <- as.integer(factor(round(tileCx$x_microns)))
      tileCx$Y <- as.integer(factor(round(tileCx$y_microns)))
      tileCx$IDX <- tileCx$X + (tileCx$Y - 1) * (max(tileCx$X))
      tileCx$FOV <- as.integer(factor(tileCx$fov, levels=params$fov_names))

      if( anchorMode=='grid' ){
        tileMatrix <- try( suppressWarnings( matrix(
          tileCx[match(1:max(tileCx$IDX), tileCx$IDX), 'FOV'],
          nrow=max(tileCx$X), ncol=max(tileCx$Y), byrow=T) ) )
        if(inherits(t, 'try-error')){
          stop('Unable to perform grid based selection of FOVs: consider random sampling or user input instead!')
        }
        m <- tileMatrix
        m <- m[-c(1, nrow(m)),]
        m <- m[,-c(1, ncol(m))]
        chosen_fovs <- as.vector( unlist(apply(m, 1, function(x){
          x[is.na(x)] <- -1
          bool <- x %in% na.omit(diag(m))
          if(all(bool==F)){ return(-1) }
          if(which(bool) %% 2 == 0){
            x <- x[-c(1:anchorOffset)]
          }
          y <- x[c(T, rep(F, anchorSpacing))]
          y <- y[y>0]
          return(y)
        })))
        chosen_fovs <- chosen_fovs[ chosen_fovs > 0 ]
        anchors <- sort( tileCx$fov[tileCx$FOV %in% chosen_fovs] )
      }

      if(anchorMode=='sample'){
        if(nSamples >= length(unique(tileCx$fov))){
          nSamples <- ceiling(length(unique(tileCx$fov)) / 10)
          warning('nSamples too big...Defaulting to 10% of FOVs (', nSamples, ')...')
        }
        anchors <- sample( unique(tileCx$fov), nSamples, replace = F)
      }

    }
  }
  if( !all(anchors %in% params$fov_names) ){
    stop('Anchors are not present in params$fov_names!')
  }
  params$anchors <<- anchors
  if(verbose){ message(paste0('\nUsing ', length(anchors), ' anchor FOVs...')) }


  ## Create anchors plot for visualisation
  tmp <- global_coords[!duplicated(global_coords$fov),]
  tmp$fovi <- as.integer(factor(tmp$fov, levels=params$fov_names))
  tmp$anchor <- tmp$fov %in% anchors
  p <-
    ggplot2::ggplot(tmp, ggplot2::aes(x=x_microns, y=y_microns, label=fovi, colour=anchor)) +
    ggplot2::scale_colour_manual(values=c('TRUE' = 'red', 'FALSE' = 'black')) +
    ggplot2::geom_text(size=8) +
    ggplot2::theme_minimal(base_size=20) +
    ggplot2::coord_fixed() + ggplot2::scale_y_reverse() +
    ggplot2::xlab('X (microns)') + ggplot2::ylab('Y (microns)') +
    ggplot2::theme(legend.position = 'none')
  params$anchors_plot <<- p
  params$anchors_metrics <<- names(imageFunctions)


  if(returnTroubleShootPlots){
    troubleshootPlots <<- new.env()
  }

  ## Check if desired functions have been run
  for( i in 1:length(imageFunctions) ){
    current_name <- names(imageFunctions)[i]
    out_file <- paste0(out_dir, toupper(current_name), '.csv')
    if( file.exists(out_file) & params$resumeMode & !freshAnchors ){
      info <- data.table::fread(out_file, data.table = F)[,1]
      if(length(info) == ncol(codebook)){
        params[[current_name]] <- info
      }
    }
  }
  if(freshAnchors){
    skip <- rep(F, length(imageFunctions))
  }else{
    skip <-  (names(imageFunctions)) %in% ls(params)
  }

  ## If any functions left, perform
  # STACKWARNED = F
  if( (length(imageFunctions) > 0) & (sum(skip) < length(imageFunctions)) ){

    imageFunctions <- imageFunctions[!skip]
    ## Initialise
    results <- list()
    for(fi in 1:length(imageFunctions)){
      results[[names(imageFunctions)[fi]]] <- rep(0, ncol(codebook))
    }

    ## Loop through every bit and get metrics
    if(verbose){ message('\nObtaining anchor metrics per bit...') }

    anchorqc_file = paste0(params$out_dir, 'ANCHOR_QUANTILEQC.csv')
    if(file.exists(anchorqc_file)){
      file.remove(anchorqc_file)
    }
    # if(!file.exists(anchorqc_file)){
    #   file.create(anchorqc_file)
    # }

    for( i in 1:ncol(codebook) ){
      if(verbose){ message(paste0('\nLoading bit ', i, ' of ', ncol(codebook), '...')) }
      current_bit <- colnames(codebook)[i]
      gcx <- global_coords[
        global_coords$fov %in% anchors
        & global_coords$bit_name==current_bit,]
      channel_index <- which(
        sapply(params$resolutions$channel_order, function(x){
          grepl(tolower(x), tolower(current_bit))
        }))

      # Load images
      mc <- match.call()
      if( !(any(names(mc) %in% 'subset')) ){
        imList <- readImageList(
          fileNames = unique(gcx$image_file),
          channelIndex = channel_index,
          ...
        )
      }else{
        if( 'c' %in% names(mc$subset) ){
          stop('Do not specify a channel for getAnchorParams! This is handled under the hood!')
        }
        imList <- readImageList(
          fileNames = unique(gcx$image_file),
          channelIndex = channel_index,
          ...
        )
      }
      nchannels <- length(params$resolutions$channel_order)
      imList <- setNames( lapply(imList, function(im){
        if(is.list(im)){ im <- im[[1]] }
        im <- array(as.vector(im), dim = dim(im)[dim(im)>1])
        if( (length(dim(im))==3) & (dim(im)[3]==nchannels) ){
          im <- im[,,channel_index]
        }
        return(im)
      }), names(imList))

      ## Check images
      imCheck <- imList[[1]]
      if(is.list(imCheck)){ imCheck <- imCheck[[1]] }
      if(!is.array(imCheck)){
        stop('Unable to parse image file!')
      }
      # if( length(dim(imCheck)) == 3 ){
      #   if(!STACKWARNED){
      #     STACKWARNED = T
      #     warning(
      #       'Image stack detected: will first perform a maximum intensity projection!
      #
      #       If desiring e.g. min and max brightnesses to be obtained across the whole stack,
      #       without MIP, use imsBrightnessMin and imsBrightnessMax instead (i.e. no _MIP suffix).
      #
      #       If desiring per-slice anchor parameters, loop this function across Z slices:
      #       Add "subset = list( z = {User chosen Z slice} )" to the parameters of this function!'
      #       )
      #   }
      #
      # }
      rm(imCheck)

      for( fi in 1:length(imageFunctions) ){
        current_func <- names(imageFunctions)[fi]
        if(verbose){ message(paste0('Getting ', current_func, '...')) }
        imFunc <- imageFunctions[[fi]]
        sumFunc <- summaryFunctions[[current_func]]
        intermediateResults <- lapply( imList, function(im){
          if(is.list(im)){ im <- im[[1]] }
          imFunc(im, ...)
        })
        intermediateResults <- as.vector(unlist(intermediateResults))
        intermediateResults <- intermediateResults[!is.na(intermediateResults)]
        identifiedThreshold <- sumFunc(intermediateResults, ...)
        results[[current_func]][i] <- identifiedThreshold

        ## Summary table
        summary_report <- lapply(imList, function(imx){
          sum(imx < identifiedThreshold) / length(imx)
        })
        summary_report <- c( 'metric' = current_func, 'bit' = i, summary_report )
        data.table::fwrite(summary_report, file = anchorqc_file, append = T)

        if(returnTroubleShootPlots){

          ## SLOW
          # dfp <- do.call(rbind, lapply( 1:length(imList), function(imidx){
          #   imxi <- imList[[imidx]]
          #   dfx <- data.frame(
          #     'intensity' = as.vector(imxi),
          #     'fov' = names(imList)[imidx])
          #   return(dfx)
          # }) )
          #
          # plot_title <- paste0( current_func, ': Bit ', i )
          # p <-
          #   ggplot2::ggplot( data=dfp, ggplot2::aes(x=intensity) ) +
          #   ggplot2::geom_histogram(fill='black', bins=100) +
          #   ggplot2::facet_wrap( ~factor(fov), scales = 'free_y') +
          #   ggplot2::theme_minimal(base_size=14) +
          #   ggplot2::xlab('Intensity') + ggplot2::ylab('Count') +
          #   ggplot2::geom_vline( xintercept = identifiedThreshold, colour='red', linetype='dashed') +
          #   ggplot2::ggtitle( plot_title )

          ## Default geom_histogram is too slow, so swap out with this:
          ## From https://stackoverflow.com/questions/56607124/plotting-histogram-of-a-big-matrix-in-ggplot2-is-20x-slower-than-base-hist
          dfp <- do.call(rbind, lapply( 1:length(imList), function(imidx){
            imxi <- imList[[imidx]]
            nbreaks = 50
            imName <- names(imList)[imidx]
            if(is.null(imName)){
              imName <- sprintf(paste0('%0', nchar(length(imList)), 'd'), imidx)
            }
            his <- hist(as.vector(imxi), plot=F, breaks=nbreaks)
            dfx <- data.frame(
              'fov' = imName,
              'xmin' = his$breaks[-length(his$breaks)],
              'xmax' = his$breaks[-1],
              'ymin' = 0,
              'ymax' = his$counts
              )
            return(dfx)
          }) )
          plot_title <- paste0( current_func, ': Bit ', i )
          p <-
            ggplot2::ggplot( data=dfp, ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) ) +
            ggplot2::geom_rect(fill='black', colour='black', size=0.5) +
            ggplot2::facet_wrap( ~factor(fov), scales = 'free_y') +
            ggplot2::theme_minimal(base_size=14) +
            ggplot2::xlab('Intensity') + ggplot2::ylab('Count') +
            ggplot2::geom_vline( xintercept = identifiedThreshold, colour='red', linetype='dashed') +
            ggplot2::ggtitle( plot_title )
          troubleshootPlots[[ current_func ]][[ i ]] <- p
        }

      }
    }

    ## Once results ready, save
    for( i in 1:length(results) ){
      current_name <- names(results)[i]
      out_file <- paste0(out_dir, toupper(current_name), '.csv')
      write.csv(results[[i]], out_file, row.names = F)
      params[[current_name]] <<- results[[i]]
    }
  }

  if(verbose){
  message(paste0(
    '\nAnchor metrics obtained: ', paste(params$anchors_metrics, collapse=', ')
    ))
  }
}

##

maxIntensityProject <- function(im, zDim = NULL){

  if(!is.null(zDim)){
    zDim = as.integer(zDim)
    if(zDim <1 | zDim>3){
      stop('zDim needs to be between 1 and 3! When im has 4 dimensions / is a list of 3D arrays, the first dimension is ignored!')
    }
  }

  if(is.list(im)){
    imageDimension <- c(length(im), dim(im[[1]]))
  }else{
    imageDimension <- dim(im)
  }

  ## If zDim not specified, guess which is zDim
  if(length(imageDimension) > 4){
    stop('Can only except a max 4 dimensional object for im (a list of 3D arrays is considered 4D)!')
  }

  if( !(length(imageDimension) %in% c(2, 4, 3) ) ){
    warning(
      paste0('Images loaded of unusual dimension size: ',
             paste(imageDimension, collapse=' x '),
             '...Skipping MIP...')
    )
    return( im )
  }

  if(length(imageDimension)==2){
    warning('2D matrix provided...Skipping MIP...')
    return(im)
  }

  if(length(imageDimension)==4){
    result <- list()
    for( ci in 1:imageDimension[1] ){
      if(is.list(im)){
        imx <- im[[ci]]
      }else{
        imx <- im[ci,,,]
      }
      result[[ci]] <- maxIntensityProject(imx, zDim=zDim)
    }
    names(result) <- names(im) #If there are names for im
  }

  if(length(imageDimension)==3){

    if(is.list(im) | imageDimension[1]==1){
      result <- list()
      for( ci in 1:imageDimension[1] ){
        if(is.list(im)){
          imx <- im[[ci]]
        }else{
          imx <- im[ci,,]
        }
        result[[ci]] <- maxIntensityProject(imx, zDim=zDim)
      }
    }else{
      if(is.null(zDim)){
        zDim = 3
      }
      newDim = imageDimension[-zDim]
      if(length(newDim)!=2){
        stop('zDim needs to be between 1 and 3! When im has 4 dimensions / is a list of 3D arrays, the first dimension is ignored!')
      }
      result <- array(0, dim = newDim)
      for( i in 1:imageDimension[zDim] ){
        if(zDim==1){
          ref <- im[i,,]
        }
        if(zDim==2){
          ref <- im[,i,]
        }
        if(zDim==3){
          ref <- im[,,i]
        }
        bool <- (ref - result) > 0
        result[bool] <- ref[bool]
      }
    }

  }

  return(result)
}

##


imsBrightnessMin <- function(
    im,
    smallBlur,
    bigBlur,
    ...
){
  if(missing(smallBlur)){ smallBlur = 2 }
  if(missing(bigBlur)){ bigBlur = 5 }
  imBlur <- as.array(imager::isoblur(suppressWarnings(imager::as.cimg(im)), smallBlur))
  LoG <- imLaplacianOfGaussian(im, smallBlur, bigBlur)
  thresh <- findThreshold(
    imBlur, thresh_Quantile,
    labels=LoG<0,
    quantileFalse = 0.5,
    quantileTrue = 0.01)
  return(thresh)
}

imsBrightnessMin_MIP <- function(
    im,
    zDim = 3,
    ...
){

  ## Check if image is a stack
  if( length(dim(im)) != 2  ){
    im <- suppressWarnings( maxIntensityProject(im, zDim=zDim) )
  }

  result <- imsBrightnessMin(im, ...)
  return(result)
}


##

imsBrightnessMax <- function(
    im,
    smallBlur,
    bigBlur,
    ...
){
  if(missing(smallBlur)){ smallBlur = 2 }
  if(missing(bigBlur)){ bigBlur = 5 }
  imBlur <- as.array(imager::isoblur(suppressWarnings(imager::as.cimg(im)), smallBlur))
  lowthresh <- imsBrightnessMin(im, smallBlur, bigBlur)
  imBlur <- imNormalise(imBlur, lowthresh)
  detHess <- imHessianDeterminant(imBlur, smallBlur)
  bestSpots <-
    (detHess>quantile(detHess[detHess>0], 0.99)) &
    (imBlur>quantile(imBlur[imBlur>0], 0.99))
  thresh <- findThreshold(
    imBlur, thresh_Quantile,
    labels = bestSpots,
    quantileFalse = 0.95,
    quantileTrue = 0.95)
  return(thresh)
}

imsBrightnessMax_MIP <- function(
    im,
    zDim=3,
    ...
){

  ## Check if image is a stack
  if( length(dim(im)) != 2  ){
    im <- suppressWarnings( maxIntensityProject(im, zDim=zDim) )
  }

  result <- imsBrightnessMax(im, ...)
  return(result)
}

##

conservativeMean <- function(
    values,
    nSDs = 0.9,
    ...
){
  values <- values[!is.na(values) & !is.infinite(values)]
  lowbound <- mean(values) - ( sd(values) * nSDs )
  hibound <- mean(values) + ( sd(values) * nSDs )
  final <- mean( values[(values>lowbound) & (values<hibound)] )
  return(final)
}
