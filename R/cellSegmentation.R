## This script is for cell segmentation functions

require(data.table)
require(EBImage)
require(imager)
require(reticulate)

establishCellSegModel <- function(
    model = 'cellpose',
    gpu = F,
    cellposePreTrainedModel = NULL,
    cellposeModelType = 'nuclei',
    stardistModelType = '2D_versatile_fluo',
    params = get('params', envir = globalenv())
){
  ## Check that entered parameters are ok
  model = tolower(model)
  acceptable_models <- c('cellpose', 'stardist')
  if(!all(model %in% acceptable_models)){
    stop('Invalid model choice!')
  }
  if(!is.null(cellposePreTrainedModel)){
    if(!file.exists(cellposePreTrainedModel)){
      stop('cellposePreTrainedModel specified does not exist!')
    }
  }
  
  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose 
  }else{
    verbose = T
  }

  ## Check that python is available
  t <- suppressWarnings( try( reticulate::use_python(params$python_location) ) )
  if(inherits(t, 'try-error')){
    stop('Invalid params$python_location specified!')
  }
  # use_python(params$python_location)
  reticulate::use_condaenv(params$python_location, required = TRUE)

  ## Create a new environment
  cell_seg_env <- new.env()

  message('\nPreparing cell segmentation model(s)...')

  ## If cellpose
  if('cellpose' %in% model){
    cellpose <- try(reticulate::import('cellpose'))
    if(inherits(cellpose, 'try-error')){
      stop('Cellpose not available in your environment!')
    }
    if(!is.null(cellposePreTrainedModel)){
      message('User-pretrained Cellpose model...')
      cell_seg_model <- try( cellpose$models$CellposeModel(
        gpu = gpu, pretrained_model = cellposePreTrainedModel) )
    }else{
      message('Cellpose packaged model...')
      cell_seg_model <- try( cellpose$models$CellposeModel(
        gpu = gpu, model_type = cellposeModelType) )
    }
    if(inherits(cell_seg_model, 'try-error')){
      stop('Invalid parameterisation of cell segmentation model!')
    }
    cell_seg_output <- c('masks', 'flows', 'styles', 'diams')
    cell_seg_env[['CELLPOSE']] <- list(
      'model' = cell_seg_model,
      'outputs' = cell_seg_output
    )
  }

  ## If stardist
  if('stardist' %in% model){
    stardist <- try(reticulate::import('stardist'))
    if(inherits(stardist, 'try-error')){
      stop('Stardist not available in your environment!')
    }
    message('Stardist packaged model...')
    cell_seg_model <- try( stardist$models$StarDist2D$from_pretrained(stardistModelType) )
    if(inherits(cell_seg_model, 'try-error')){
      stop('Invalid parameterisation of cell segmentation model!')
    }
    cell_seg_model$config$use_gpu <- gpu
    cell_seg_output <- c('masks', 'details')
    cell_seg_env[['STARDIST']] <- list(
      'model' = cell_seg_model,
      'outputs' = cell_seg_output
    )
  }

  cellSeg <<- cell_seg_env
  if(verbose){ message(paste0('Cell segmentation with ', paste(ls(cell_seg_env), collapse=', '), '...')) }
}

##

runCellSegModel <- function(
    im,
    masksOnly = T,
    recordLastImage = T,
    cellSeg = get('cellSeg', envir = globalenv()),
    params = get('params', envir = globalenv()),
    ...
){

  results <- list()

  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose 
  }else{
    verbose = T
  }
  
  ## If im is a list of images
  if(is.list(im)){
    for(i in 1:length(im)){
      if(verbose){ message(paste0(i, ' of ', length(im), '...')) }
      results[[i]] <-
        runCellSegModel(
          im[[i]],
          masksOnly = masksOnly,
          params = params,
          ...
        )
    }
    names(results) <- names(im)
    return(results)
  }

  ## If im is a single image
  if( recordLastImage ){
    params$cell_seg_image <<- im
  }
  for( modelType in ls(cellSeg) ){
    if(verbose){ message(paste0('Running ', modelType, '...')) }

    model <- cellSeg[[modelType]]$model
    outputs <- cellSeg[[modelType]]$outputs
    if( modelType=='CELLPOSE' ){
      res <- model$eval(im, ...)
    }
    if( modelType=='STARDIST' ){
      res <- model$predict_instances(im, ...)
    }
    names(res) <- outputs[1:length(res)]

    if(masksOnly){
      res <- res$masks
    }
    results[[modelType]] <- res
  }

  if(length(results)==1){
    results <- results[[1]]
  }

  return(results)
}

###

consolidateCellSegMasks <- function(
    maskList,
    minFlags = 1,
    cellSeg = get('cellSeg', envir = globalenv()),
    params = get('params', envir = globalenv())
){
  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose 
  }else{
    verbose = T
  }
  if(!is.list(maskList)){
    return(maskList)
  }
  while(is.list(maskList)){
    maskListOld <- maskList
    maskList <- unlist(maskList, recursive = F, use.names = T)
  }
  maskList <- Reduce('+', lapply(maskListOld, function(im){
    hasMask <- im>0
    grad <- imager::imgradient(imager::as.cimg(im))
    gradmag <- as.matrix( sqrt(grad$x^2+grad$y^2) )
    maskOutlines <- gradmag > 0 #Note that this shrinks cells by 1 pixel at least
    hasMask[maskOutlines] <- 0
    return(hasMask)
    }) )
  if(!is.matrix(maskList)){
    stop('Unable to parse maskList! Expecting a matrix after successive unlist()!')
  }
  final <- EBImage::bwlabel(maskList>=minFlags)
  return(final)
}

###

saveCellSegMasks <- function(
    masks,
    im = NULL,
    currentFOVName = NULL,
    cellSeg = get('cellSeg', envir = globalenv()),
    params = get('params', envir = globalenv())
){
  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose 
  }else{
    verbose = T
  }
  if(!is.matrix(masks)){
    stop('Expecting a matrix for masks!')
  }
  if(!is.null(im)){
    if(!is.matrix(im)){
      stop('Expecting a matrix for im!')
    }
  }else{
    if( is.matrix(params$cell_seg_image) ){
      im <- params$cell_seg_image
    }
  }
  if(is.null(currentFOVName)){
    currentFOVName <- params$current_fov
    if(is.null(currentFOVName)){
      stop('currentFOVName not specified and params$current_fov does not exist!')
    }
  }
  cellSegType <- paste(sort(ls(cellSeg)), collapse='')
  df <- as.data.frame(imager::as.cimg(masks))
  df$value[is.na(df$value) | is.infinite((df$value))] <- 0
  df <- df[df$value>0,]
  
  if(nrow(df)==0){
    #warning('No mask found...Will not save...')
    if(verbose){message('No mask found...Will not save...')}
    return()
  }

  ## Prepare dataframe to be saved
  gcx <- params$global_coords
  gcx <- gcx[gcx$fov == currentFOVName,]
  gcx <- as.numeric(unlist(gcx[!duplicated(gcx$fov),c('x_microns', 'y_microns')]))
  per_pixel_microns <- as.numeric(params$resolutions$per_pixel_microns)
  if( !is.null(gcx) & length(gcx)==2 & !is.null(per_pixel_microns) ){
    outdf <- data.frame(
      'Xm' = df$x * per_pixel_microns[1] + gcx[1],
      'Ym' = df$y * per_pixel_microns[2] + gcx[2],
      'fov' = currentFOVName,
      'WX' = df$x, 'WY' = df$y, 'cell' = df$value
    )
  }else{
    warning('Unable to parse global coordinate information: returning Xm and Ym = NA for now')
    outdf <- data.frame(
      'Xm' = NA,
      'Ym' = NA,
      'fov' = currentFOVName,
      'WX' = df$x, 'WY' = df$y, 'cell' = df$value
    )
  }

  ## Save mask compactly
  out_file <- paste0(
    params$out_dir,
    'CELLSEG_', cellSegType, '_',  currentFOVName, '.csv.gz'
  )
  data.table::fwrite(outdf, out_file, row.names = F )


  ## Save image
  if(!is.null(im)){
    out_file <- paste0(
      params$out_dir,
      'REGIM_', currentFOVName, '.png'
    )
    imager::save.image(imager::as.cimg(im), out_file, quality = 1)
  }

}




