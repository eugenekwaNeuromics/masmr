## This contains functions for stitching

require(ggplot2)
require(data.table)

readImagesForStitch <- function(
    currentFOVName,
    subDirectory = 'IM',
    loadProcessedImages = F,
    registerTo = 1,
    zSlice = NULL,
    imageFunction = function(im){ return(im) },
    params = get('params', envir=globalenv()),
    ...
){

  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose
  }else{
    verbose = T
  }

  ## Prepare params
  if(!is.function(imageFunction)){
    stop('imageFunction needs to be a function!')
  }
  if(missing(currentFOVName)){
    currentFOVName <- params$current_fov
    if(is.null(currentFOVName)){
      stop('currentFOVName not specified and params$current_fov not found!')
    }
  }
  if(length(subDirectory) > 1){
    warning( 'More than one subDirectory detected...Taking first entry...' )
    subDirectory <- subDirectory[1]
  }
  if(!dir.exists( paste0(params$parent_out_dir, '/', subDirectory))){
    stop('Invalid subDirectory: needs to be a directory in params$parent_out_dir!')
  }
  f <- list.files( paste0(params$parent_out_dir, '/', subDirectory),
                   pattern = 'GLOBALCOORD', full.names = T )[1]
  if(is.na(f)){
    stop('Unable to find GLOBALCOORD file in chosen subDirectory!')
  }
  globalcoords <- data.table::fread(f, data.table = F)
  f <- list.files( paste0(params$parent_out_dir, '/', subDirectory),
                   pattern = 'META_RESOLUTION', full.names = T )[1]
  if(is.na(f)){
    stop('Unable to find META_RESOLUTION file in chosen subDirectory!')
  }
  current_res_info <- readLines(f)
  n <- current_res_info[c(TRUE, FALSE)]
  val <- current_res_info[c(FALSE, TRUE)]
  val <- lapply(val, function(v) strsplit(v, ' ')[[1]])
  resolutions <- setNames(val, n)

  ## Find neighbours of image
  gcxi <- globalcoords
  gcxi$ref <- bool <- gcxi$fov==currentFOVName
  gcxi$x <- as.integer( factor(rank( round(abs( gcxi[,'x_microns'] - gcxi[bool,'x_microns'] )), ties.method = 'average')) )-1
  gcxi$y <- as.integer( factor(rank( round(abs( gcxi[,'y_microns'] - gcxi[bool,'y_microns'] )), ties.method = 'average')) )-1
  # gcxi$block_dist <- rowSums(gcxi[,c('x', 'y')])
  # gcxi$nn <- gcxi$block_dist==1

  ## Calculate physical distance
  cx <- (gcxi$x_microns + 1i * gcxi$y_microns)
  dist <- round( Mod(cx - cx[bool]), digits=2 ) #Within 0.01 microns
  dist_thresh <- sort(unique(dist[dist>0]))[4]

  miny_whenx0 <- sort(gcxi[gcxi$x==0 & dist<=dist_thresh,'y'])
  miny_whenx0 <- min(miny_whenx0[miny_whenx0>0])
  minx_wheny0 <- sort(gcxi[gcxi$y==0 & dist<=dist_thresh,'x'])
  minx_wheny0 <- min(minx_wheny0[minx_wheny0>0])
  gcxi$nn <- (gcxi$x==0 & gcxi$y==miny_whenx0) | (gcxi$y==0 & gcxi$x==minx_wheny0)

  ## Parse image read info
  gcx <- gcxi[gcxi$fov==currentFOVName,]
  gcx <- gcx[order(gcx$image_file),]
  if(length(gcx$bit_name) < registerTo){
    stop( paste0('Invalid registerTo: expecting values between 1 and ', length(gcx$bit_name), '!') )
  }
  chosen_bit <- gcx$bit_name[registerTo]
  gcxi <- gcxi[gcxi$bit_name==chosen_bit,]
  channel_index <- which(
    sapply(resolutions$channel_order, function(x){
      grepl(tolower(x), tolower(chosen_bit))
    }))
  gcxi <- gcxi[!duplicated(gcxi$fov),] #To account for scenarios with multiple Z depths
  chosen_fovs <- c( gcxi[gcxi$ref, 'fov'], gcxi[gcxi$nn, 'fov']  )
  # fov_idx <- match(chosen_fovs, resolutions$fov_names)
  fs <- c(gcxi[gcxi$ref, 'image_file'], gcxi[gcxi$nn, 'image_file'] )

  ## Create a temporary environment
  tmp_param <- new.env()
  tmp_param$global_coords <- gcxi
  tmp_param$im_format <- params$im_format
  tmp_param$fov_names <- resolutions$fov_names
  tmp_param$resolutions <- resolutions
  tmp_param$out_dir <- paste0(params$parent_out_dir, '/STITCH/')
  tmp_param$sub_directory <- subDirectory
  tmp_param$current_fov <- currentFOVName
  tmp_param$verbose <- verbose

  ## Read images
  if(loadProcessedImages){
    regims <- list.files( paste0(params$parent_out_dir, '/', subDirectory), pattern='^REGIM_', full.names = T)
    if(length(regims) > 0){
      names(regims) <- gsub('^REGIM_|[.]png', '', basename(regims))
      regims <- regims[match(chosen_fovs, names(regims))]
      names(regims) <- chosen_fovs
      if( all(is.na(regims)) ){
        warning('Unable to find any processed images: loaded raw images instead...')
        loadProcessedImages = F
      }
      if(loadProcessedImages){
        fs <- regims
        imList <- list()
        for( i  in 1:length(regims) ){
          if(is.na(regims[i])){
            imList[[names(regims)[i]]] <- NA
            next
          }
          imList[[names(regims)[i]]] <-
            readImage( regims[i],
                       nrows = as.numeric(resolutions$xydimensions_pixels)[1],
                       ncols = as.numeric(resolutions$xydimensions_pixels)[2],
                       # nzs = 1, nchannels = 1 #Shouldn't need to specify
                       ...)
        }
      }
    }else{
      warning('Unable to load processed images! Loading raw images...')
      loadProcessedImages = F
    }
  }

  if(!loadProcessedImages){

    if(as.numeric(resolutions$zslices)>1){

      if(!is.null(zSlice)){
        warning( paste0("Detected multiple z-slices...Defaulting to zSlice = ", zSlice, '...'))
        imList <- readImageList(
          fileNames = fs,
          params = tmp_param,
          channelIndex = channel_index,
          zIndex = zSlice,
          ...
        )
      }else{
        warning("Detected multiple z-slices...Defaulting to MIP...")
        imList <- readImageList(
          fileNames = fs,
          params = tmp_param,
          channelIndex = channel_index,
          ...
        )
        imList <- maxIntensityProject(imList)
      }

    }else{
      imList <- readImageList(
        fileNames = fs,
        params = tmp_param,
        channelIndex = channel_index,
        ...
      )
    }

  }

  tmp_param$fovs_loaded <- setNames(fs, chosen_fovs)
  tmp_param$processed_images_loaded <- loadProcessedImages
  stitchParams <<- tmp_param

  ## Process each image
  if(verbose){ message('\nProcessing images...') }

  imList <- setNames( lapply(imList, function(im){
    if( is.list(im) ){ im <- im[[1]] }
    if( any(is.na(im)) ){ return(NA) }
    imageFunction(im, ...)
    }), names(imList) )
  imList <- imList[!is.na(imList)]
  imList <- imList[c(currentFOVName, names(imList)[names(imList) != currentFOVName])]

  return(imList)
}

###

stitch_troubleshootPlots <- function(
    imList,
    stitchResults,
    alphaRange = c(0, 0.75),
    shiftColumn = NULL,
    marginToPlot = NULL,
    stitchParams = get('stitchParams', envir = globalenv())
){
  if (is.logical(stitchParams$verbose)) {
    verbose = stitchParams$verbose
  }else{
    verbose = T
  }

  if(is.data.frame(stitchResults) | is.matrix(stitchResults)){

    if(!is.null(shiftColumn)){
      stitchResults <- as.complex(stitchResults[,shiftColumn])
    }else{
      ## Attempt to infer which column is the shift column
      is_complex <- rep(F, ncol(stitchResults))
      for(i in 1:ncol(stitchResults)){
        is_complex[i] <- suppressWarnings( !all(is.na( as.complex(stitchResults[,i])) ))
      }
      if(sum(is_complex)<1){
        stop('Unable to parse stitchResults to find the column housing the shift vector: kindly specify shiftColumn!')
      }
      if(sum(is_complex)>1){
        warning('More than one possible shiftColumn...Picking the first...')
        stitchResults <- as.complex( stitchResults[,which(is_complex)[1]] )
      }
      if(sum(is_complex)==1){
        stitchResults <- as.complex( stitchResults[,is_complex] )
      }
    }
  }

  if(!is.complex(stitchResults)){
    stop('Provide stitchResults as a dataframe (with shifts in one column) or a complex vector!')
  }
  shiftVector <- stitchResults

  ## ImList is always assumed to have the reference FOV first, followed by neighbours
  if( length(shiftVector) != (length(imList) - 1) ){
    stop('Expecting a list of images, where first image is the reference, followed by neighbours (in order)')
  }
  refim <- imList[[1]]
  nnList <- imList[-1]

  if(verbose){ message('\nGenerating QC plots for each neighbour in stitching...')}
  plotList <- list()
  for( i in 1:length(shiftVector) ){
    if(verbose){ message( paste0(i, ' of ', length(shiftVector), '...') )}
    shiftvi <- shiftVector[i]
    nnim <- nnList[[i]]
    imName <- names(imList)[i]
    imDirection <- which.max(abs(c(Re(shiftvi), Im(shiftvi))))
    imSize <- dim(refim)[imDirection]
    if(is.null(marginToPlot)){
      marginToPlot <- pmin( round( 2 * abs( abs( c(Re(shiftvi), Im(shiftvi))[imDirection] ) - imSize ) ), imSize )
    }
    imDirection <- sign( c(Re(shiftvi), Im(shiftvi))[imDirection] ) * imDirection
    imDirection <- which(c(-2,-1,1,2)==imDirection)

    imDirectionName <- switch(imDirection, 'NORTH', 'WEST', 'EAST', 'SOUTH')
    if(is.null(imName)){
      imName <- imDirectionName
    }else{
      imName <- as.character(imName)
    }

    refImage <- suppressWarnings( as.data.frame(imager::as.cimg(refim)) )
    querImage <- suppressWarnings( as.data.frame(imager::as.cimg(nnim)) )

    refImage <- switch(
      imDirection,
      refImage[refImage$y<(marginToPlot),],
      refImage[refImage$x<(marginToPlot),],
      refImage[refImage$x>(imSize-marginToPlot),],
      refImage[refImage$y>(imSize-marginToPlot),]
      )

    querImage <- switch(
      imDirection,
      querImage[querImage$y>(imSize-marginToPlot),],
      querImage[querImage$x>(imSize-marginToPlot),],
      querImage[querImage$x<(marginToPlot),],
      querImage[querImage$y<(marginToPlot),]
    )

    refImage$image <- 'Reference'
    querImage$image <- 'Query'

    refImage$value <- refImage$value - min(refImage$value)
    refImage$value <- refImage$value / max(refImage$value)
    querImage$value <- querImage$value - min(querImage$value)
    querImage$value <- querImage$value / max(querImage$value)

    p_unstitched <-
      ggplot2::ggplot() +
      ggplot2::geom_raster(
        data=refImage,
        ggplot2::aes(x=x, y=y, alpha=value, fill = image),
      ) +
      ggplot2::geom_raster(
        data=querImage,
        ggplot2::aes(
          x=x + switch(imDirection, 0, -imSize, imSize, 0),
          y=y + switch(imDirection, -imSize, 0, 0, imSize),
          alpha=value, fill = image)) +
      ggplot2::scale_fill_manual(name='', values=c('Reference' = 'red', 'Query' = 'blue')) +
      ggplot2::scale_alpha_continuous( range = alphaRange ) +
      ggplot2::guides( alpha = 'none', fill = ggplot2::guide_legend(ncol=1) ) +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_void( base_size= 14 ) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::coord_fixed()

    p_stitched <-
      ggplot2::ggplot() +
      ggplot2::geom_raster(
        data=refImage,
        ggplot2::aes(x=x, y=y, alpha=value, fill = image),
      ) +
      ggplot2::geom_raster(
        data=querImage,
        ggplot2::aes(
          x=x + Re(shiftvi), y=y + Im(shiftvi),
          alpha=value, fill = image)) +
      ggplot2::scale_fill_manual(name='', values=c('Reference' = 'red', 'Query' = 'blue')) +
      ggplot2::scale_alpha_continuous( range = alphaRange ) +
      ggplot2::guides( alpha = 'none', fill = ggplot2::guide_legend(ncol=1) ) +
      ggplot2::scale_y_reverse() +
      ggplot2::theme_void( base_size= 14 ) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::coord_fixed()

    plotList[[ imName ]] <- list('UNSTITCHED' = p_unstitched, 'STITCHED' = p_stitched)
  }
  return(plotList)
}

###

stitchImages <- function(
    imList,
    registerTo = 1,
    returnTroubleShootPlots = FALSE,
    stitchParams = get('stitchParams', envir=globalenv())
){

  ## Get verbosity
  if(is.logical(stitchParams$verbose)){
    verbose = stitchParams$verbose
  }else{
    verbose = T
  }

  if(length(registerTo) != 1 | is.na(as.numeric(registerTo))){
    stop('Invalid registerTo parameter!')
  }
  registerTo <- as.numeric(registerTo)
  if( (registerTo > length(imList)) | (registerTo <= 0) ){
    stop( paste0('registerTo needs to between 1 and ', length(imList), '!'))
  }
  refx <- imList[[registerTo]]
  ZDETECTED = F
  if(length(dim(refx))>2){
    ZDETECTED = T
  }
  nnlist <- imList[-registerTo]
  if(ZDETECTED){
    warning('Z stack detected...Defaulting to MIP...')
    refx <- maxIntensityProject(refx)
    nnlist <- setNames( maxIntensityProject(nnlist), names(nnlist) )
  }

  if(length(nnlist)==0){
    if(verbose){ message('No neighbours found...Returning empty dataframe...') }
    return(data.frame())
  }
  if( all(is.na(refx)) | all(is.null(refx)) ){
    if(verbose){ message('No reference image found...Returning empty dataframe...') }
    return(data.frame())
  }

  gcxi <- stitchParams$global_coords
  if(!all(names(nnlist) %in% gcxi$fov)){
    stop('Names of imList needs to be in stitchParams$global_coords$fov!')
  }
  per_pixel_microns <- as.numeric(stitchParams$resolutions$per_pixel_microns)

  if(verbose){ message('\nAligning adjacent FOVs...') }
  shifts <- c()
  # coord <- getRasterCoords(refx) - sum(dim(refx) * c(1, 1i)) + (1+1i)
  for(imidxi in 1:length(nnlist)){
    if(verbose){ message(paste0(imidxi, ' of ', length(nnlist), '...'), appendLF = F) }
    epcx <- gcxi[gcxi$fov==names(nnlist)[imidxi],c('x_microns', 'y_microns')] - gcxi[gcxi$ref,c('x_microns', 'y_microns')]
    epcx <- (epcx['x_microns']/per_pixel_microns[1]) + 1i * (epcx['y_microns']/per_pixel_microns[2])
    epcx <- as.complex(epcx)

    dist <- sum(dim(refx) * c(1, 1i)) - (abs(Re(epcx)) + abs(Im(epcx))*1i)
    dist <- pmin(Re(dist), Im(dist))
    dist <- dist

    refim <- refx
    querim <- nnlist[[imidxi]]

    if( abs(Re(epcx)) > abs(Im(epcx)) ){
      if( Re(epcx) > 0 ){
        refim <- refim[(nrow(refim) - dist + 1):nrow(refim),]
        querim <- querim[1:dist,]
        offset <- sum(dim(refx) * c(1, 0i))
      }else{
        refim <- refim[1:dist,]
        querim <- querim[(nrow(querim) - dist + 1):nrow(querim),]
      }

    }else{
      # Left-right
      if( Im(epcx) > 0 ){
        refim <- refim[,(ncol(refim) - dist + 1):ncol(refim)]
        querim <- querim[,1:dist]
      }else{
        refim <- refim[,1:dist]
        querim <- querim[,(ncol(querim) - dist + 1):ncol(querim)]
      }
    }

    corr <- crossCorrelate2D(refim, querim, normalized=FALSE, pad = T)
    coord <- getRasterCoords(corr)
    coord <- coord - coord[round(nrow(coord)/2), round(ncol(coord)/2)]
    shift <- coord[which.max(corr)] + epcx
    shifts <- c(shifts, shift)
  }
  names(shifts) <- names(nnlist)

  final <-
    data.frame(
      'fov' = names(shifts),
      'shift_pixel' = as.complex(shifts)
    )

  if(returnTroubleShootPlots){
    if(verbose){ message('Generating QC plots for evaluating stitching...') }
    troubleshootPlots <<- new.env()
    troubleshootPlots[['STITCH_EVAL']] <-
      stitch_troubleshootPlots(
        imList = imList,
        stitchResults=final,
        alphaRange = c(0, 0.75),
        shiftColumn = 'shift_pixel',
        marginToPlot = NULL,
        stitchParams = get('stitchParams', envir = globalenv())
      )
  }

  return(final)
}
###

saveStitch <- function(
    stitchDF,
    stitchParams = get('stitchParams', envir=globalenv())
){
  if(!dir.exists(stitchParams$out_dir)){
    dir.create(stitchParams$out_dir)
  }
  data.table::fwrite(
    stitchDF,
    paste0(stitchParams$out_dir, 'STITCH_',
           stitchParams$current_fov, '.csv'),
    row.names = F)
}



