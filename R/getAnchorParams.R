## Function for getting max and min brightness

require(imager)
require(ggplot2)
require(data.table)

getAnchorParams <- function(
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
    anchorMode = 'grid',
    anchorSpacing = 3,
    anchorOffset = 2,
    nSamples = 10,
    params = get('params', envir = globalenv()),
    ...
){

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
      message('Loading ORDEREDCODEBOOK.csv.gz...')
    }
    codebook <- data.table::fread(f, data.table=FALSE)
    if(colnames(codebook)[1] == 'V1'){
      rownames(codebook) <- codebook$V1
      codebook$V1 <- NULL
    }
    params$ordered_codebook <<- codebook
    params$capture_order <<- colnames(codebook)
    params$isblank <<- grepl('^blank-', tolower(rownames(codebook)))
    message(paste0('Codebook has ', sum(params$isblank), ' blanks and ', sum(!params$isblank), ' genes...' ))
  }
  codebook <- params$ordered_codebook
  if(missing(out_dir)){
    out_dir <- params$out_dir
    if(is.null(out_dir)){
      stop('out_dir not specified and params$out_dir missing!')
    }
  }


  ## If anchors not specified, create
  if(missing(anchors)){
    anchors <- params$anchors

    ## Deprecated check
    freshAnchors = T
    if(!params$resumeMode | freshAnchors ){
      anchors <- NULL
    }

    if(is.null(anchors)){
      tileCx <- global_coords[!duplicated(global_coords$fov),]
      tileCx$X <- as.integer(factor(round(tileCx$x_microns)))
      tileCx$Y <- as.integer(factor(round(tileCx$y_microns)))
      tileCx$IDX <- tileCx$X + (tileCx$Y - 1) * (max(tileCx$X))
      tileCx$FOV <- as.integer(factor(tileCx$fov, levels=params$fov_names))
      tileMatrix <- matrix(
        tileCx[match(1:max(tileCx$IDX), tileCx$IDX), 'FOV'],
        nrow=max(tileCx$X), ncol=max(tileCx$Y), byrow=T)

      if( anchorMode=='grid' ){
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
  message(paste0('Using ', length(anchors), ' anchor FOVs...'))


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


  ## Check if desired functions have been run
  for( i in 1:length(imageFunctions) ){
    current_name <- names(imageFunctions)[i]
    out_file <- paste0(out_dir, toupper(current_name), '.csv')
    if( file.exists(out_file) & params$resumeMode ){
      info <- data.table::fread(out_file, data.table = F)[,1]
      if(length(info) == ncol(codebook)){
        params[[current_name]] <- info
      }
    }
  }
  skip <-  (names(imageFunctions)) %in% ls(params)

  ## If any functions left, perform
  STACKWARNED = F
  if( (length(imageFunctions) > 0) & (sum(skip) < length(imageFunctions)) ){

    imageFunctions <- imageFunctions[!skip]
    ## Initialise
    results <- list()
    for(fi in 1:length(imageFunctions)){
      results[[names(imageFunctions)[fi]]] <- rep(0, ncol(codebook))
    }

    ## Loop through every bit and get metrics
    for( i in 1:ncol(codebook) ){
      message(paste0('Loading bit ', i, ' of ', ncol(codebook), '...'))
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
          safeLoad = F,
          subset = list( c = channel_index ),
          ...
        )
      }else{
        if( 'c' %in% names(mc$subset) ){
          stop('Cannot specify a channel for getAnchorParams!')
        }
        imList <- readImageList(
          fileNames = unique(gcx$image_file),
          safeLoad = F,
          ...
        )
        nchannels <- length(params$resolutions$channel_order)
        imList <- setNames( lapply(imList, function(im){
          if(is.list(im)){ im <- im[[1]] }
          im <- array(as.vector(im), dim = dim(im)[dim(im)>1])
          if( (length(dim(im))==3) & (dim(im)[3]==nchannels) ){
            im <- im[,,channel_index]
          }
          return(im)
        }), names(imList))
      }

      ## Check images
      imCheck <- imList[[1]]
      if(is.list(imCheck)){ imCheck <- imCheck[[1]] }
      if(!is.array(imCheck)){
        stop('Unable to parse image file!')
      }
      if( length(dim(imCheck)) == 3 ){
        if(!STACKWARNED){
          STACKWARNED = T
          warning(
            'Image stack detected: will summarise values across entire stack!
          If desiring per-slice anchor parameters, loop this function across Z slices:
          Add "subset = list( z = {User chosen Z slice} )" to the parameters of this function!')
        }

      }
      rm(imCheck)

      message('')
      for( fi in 1:length(imageFunctions) ){
        current_func <- names(imageFunctions)[fi]
        message(paste0('Getting ', current_func, '...'))
        imFunc <- imageFunctions[[fi]]
        sumFunc <- summaryFunctions[[current_func]]
        intermediateResults <- lapply( imList, function(im){
          if(is.list(im)){ im <- im[[1]] }
          imFunc(im, ...)
        })
        intermediateResults <- as.vector(unlist(intermediateResults))
        intermediateResults <- intermediateResults[!is.na(intermediateResults)]
        results[[current_func]][i] <- sumFunc(intermediateResults, ...)
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

  message(paste0(
    'Anchor metrics: ', paste(params$anchors_metrics, collapse=', ')
    ))
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
