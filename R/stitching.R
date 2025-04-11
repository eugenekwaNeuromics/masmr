## This contains functions for stitching

require(data.table)

readImagesForStitch <- function(
    currentFOVName,
    subDirectory = 'IM',
    loadProcessedImages = F,
    registerTo = 1,
    imageFunction = function(im){ return(im) },
    params = get('params', envir=globalenv()),
    ...
){
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
  gcxi$block_dist <- rowSums(gcxi[,c('x', 'y')])
  gcxi$nn <- gcxi$block_dist==1

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
  chosen_fovs <- c( gcxi[gcxi$ref, 'fov'], gcxi[gcxi$nn, 'fov']  )
  fov_idx <- match(chosen_fovs, resolutions$fov_names)
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

  ## Read images
  if(loadProcessedImages){
    regims <- list.files( paste0(params$parent_out_dir, '/', subDirectory), pattern='REGIM_', full.names = T)
    if(length(regims) > 0){
      names(regims) <- strsplit2(regims, '_')[,ncol(strsplit2(regims, '_'))]
      names(regims) <- strsplit2(names(regims), '[.]')[,1]
      regims <- regims[match(chosen_fovs, names(regims))]
      names(regims) <- chosen_fovs
      if( all(is.na(regims)) ){
        warning('Unable to find any processed images: loaded raw images instead')
        loadProcessedImages = F
      }
      if(loadProcessedImages){
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
                       ...)
        }
      }
    }else{
      warning('Unable to load processed images! Loading raw images...')
      loadProcessedImages = F
    }
  }

  if(!loadProcessedImages){
    imList <- readImageList(
      fileNames = fs,
      safeLoad = F,
      params = tmp_param,
      subset = list( c = channel_index ),
      ...
    )
  }

  tmp_param$processed_images_loaded <- loadProcessedImages
  stitchParams <<- tmp_param

  ## Process each image
  message('\nProcessing images...')

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

stitchImages <- function(
    imList,
    registerTo = 1,
    stitchParams = get('stitchParams', envir=globalenv())
){
  if(length(registerTo) != 1 | is.na(as.numeric(registerTo))){
    stop('Invalid registerTo parameter!')
  }
  registerTo <- as.numeric(registerTo)
  if( (registerTo > length(imList)) | (registerTo <= 0) ){
    stop( paste0('registerTo needs to between 1 and ', length(imList), '!'))
  }
  refx <- imList[[registerTo]]
  nnlist <- imList[-registerTo]
  if(length(nnlist)==0){
    message('No neighbours found...Returning empty dataframe...')
    return(data.frame())
  }

  gcxi <- stitchParams$global_coords
  if(!all(names(nnlist) %in% gcxi$fov)){
    stop('Names of imList needs to be in stitchParams$global_coords$fov!')
  }
  per_pixel_microns <- as.numeric(stitchParams$resolutions$per_pixel_microns)

  message('\nAligning adjacent FOVs...')
  shifts <- c()
  # coord <- getRasterCoords(refx) - sum(dim(refx) * c(1, 1i)) + (1+1i)
  for(imidxi in 1:length(nnlist)){
    message(paste0(imidxi, ' of ', length(nnlist), '...'), appendLF = F)
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
      'shift_pixel' = shifts
    )

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



