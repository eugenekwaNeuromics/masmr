## Function for reading images.

require(data.table)
require(RBioFormats)

readImageMetaData <- function(
    dir = 'im_dir',
    params = get('params', envir = globalenv())
){

  writeTemplate <- function(){
    resolutions <- list(
      'per_pixel_microns' = c(0.13,0.13), #Assume square
      'xydimensions_pixels' = c(2056, 2056), #
      'xydimensions_microns' = c(2056 * 0.13, 2056 * 0.13),
      'channel_order' = c('cy7', 'cy5', 'alexa594', 'cy3'),
      'fov_names' = params$fov_names,
      'zslices' = 1
    )
    global_coord <- data.frame(
      'x_microns' = -700,
      'y_microns' = -700,
      'z_microns' = NA,
      'image_file' = gsub(params$meta_format, params$im_format, f),
      'fov' = 1,
      'bit_name_full' = sort(params$codebook_colnames),
      'bit_name_alias' = sort(params$codebook_colnames),
      'bit_name' = sort(params$codebook_colnames)
    )
    lapply(1:(length(resolutions)), function(idx){
      write(names(resolutions)[idx], paste0(params$out_dir, 'META_RESOLUTION', '.txt'), append = TRUE,  ncolumns=max(length(unlist(resolutions))))
      write(resolutions[[idx]], paste0(params$out_dir, 'META_RESOLUTION', '.txt'), append = TRUE,  ncolumns=max(length(unlist(resolutions))))
    })
    write.csv(global_coords, paste0(params$out_dir, 'GLOBALCOORD', '.csv'), row.names = F)
  }

  checkColBoolAlready = F

  ## Check that dir exists in params
  if(object.size(get(dir, envir=params))==0){
    stop('Invalid dir specified!')
  }
  dirname <- toupper(gsub('_|dir', '', tolower(dir)))
  params$out_dir <<- paste0(params$out_dir, '/', dirname, '/')
  if(!dir.exists(params$out_dir)){
    dir.create(params$out_dir)
  }

  ## Check if files exist, and if resumeMode = T, do not recreate
  createFiles = T
  if(
    length(list.files(params$out_dir, pattern='GLOBALCOORD'))>0 &
    length(list.files(params$out_dir, pattern='META_RESOLUTION'))>0
  ){
    if(params$resumeMode){
      createFiles = F
    }
  }

  # If need to create files
  # TO UPDATE WHEN EXPERIENCE MORE DATA FORMATS
  if(createFiles){

    message('Creating GLOBALCOORD and META_RESOLUTION files...')

    # Getting params from first meta data in first image directory
    f <- sort(list.files( params[[dir]][1], pattern=params$meta_format, full.names = T ))[1]

    # Try to read with RBioformats
    t <- try(RBioFormats::read.metadata(f), silent = T)
    if(inherits(t, 'try-error')){

      # Cannot read with RBioFormats
      if(!(params$meta_format %in% c('.xml'))){
        warning('Unable to read metadata with RBioFormats')
        writeTemplate()
        stop('GLOBALCOORD and META_RESOLUTION templates provided. Kindly fill in manually.')
      }

      # Niche case of .dax with .xml file - this is acquired for Dory
      if(params$meta_format == '.xml'){
        warning('Assuming this data has been acquired by Dory')
        fs <- sort(list.files( params[[dir]], pattern=params$meta_format, full.names = T ))
        names(fs) <- gsub(params$meta_format, '', basename(fs))
        resolutions <- list()
        global_coords <- list()

        for(fov in params$fov_names){
          bool <- grepl(paste0('_', fov, '$'), names(fs))
          fsub <- fs[bool]
          metas <- setNames(lapply(1:length(fsub), function(x){
            fx <- fsub[x]
            df <- readLines(fx)
            return(df)
          }), names(fsub))
          if(!all(unlist(lapply(metas, length)) == length(metas[[1]]))){
            stop("Inconsistent number of lines in meta data files!")
          }

          stage_pos_idx = 6 #From reading the metadata
          pixel_size_idx <- c(30,34) #From reading the metadata
          pp <- 1 / 7.678 #Per pixel resolution for Dory, not provided in meta data
          full_res <- metas[[1]][pixel_size_idx]
          full_res <- as.integer(strsplit2(full_res, '>|<')[,3])
          pos <- strsplit2(metas[[1]][stage_pos_idx], '>|<')[,3]
          x_coord <- -as.numeric(strsplit2(pos, ',')[,2]) #Very odd for Dory, but this is apparently the direction
          y_coord <- -as.numeric(strsplit2(pos, ',')[,1])
          bit_name_full <- gsub(paste0('_',fov, '$'), '', names(fsub))
          bit_name_alias <- tolower(bit_name_full)
          channels <- sort(unique(strsplit2(names(fsub), '_')[,1]))
          n_channels <- length(channels)
          resolutions <- list( # A bit redundant: generates a new resolutions every iteration
            'per_pixel_microns' = c(pp, pp),
            'xydimensions_pixels' = as.numeric(full_res),
            'xydimensions_microns' = as.numeric(full_res) * pp,
            'channel_order' = channels,
            'fov_names' = params$fov_names,
            'zslices' = 1
          )
          fsub <- gsub(params$meta_format, params$im_format, fsub)
          if(!all(file.exists(fsub))){
            stop( paste0('Corresponding image for meta data file of ', fov, ' does not exist!'))
          }
          global_coord <- data.frame(
            'x_microns' = x_coord,
            'y_microns' = y_coord,
            'z_microns' = NA,
            'image_file' = fsub,
            'fov' = fov,
            'bit_name_full' = bit_name_full,
            'bit_name_alias' = bit_name_alias,
            'bit_name' = bit_name_alias
            )
          global_coords[[fov]] <- global_coord
          fs <- fs[!bool]
        }
    }

    }else{

      if(!(params$meta_format %in% c('.ome.tif'))){
        warning('Metadata can be read with RBioFormats, but parsing currently not implemented.')
        writeTemplate()
        stop('GLOBALCOORD and META_RESOLUTION templates provided. Kindly fill in manually.')
      }

      if(params$meta_format == '.ome.tif'){
        warning('Assuming this data has been acquired by Triton')
        fs <- sort(list.files( params[[dir]], pattern=params$meta_format, full.names = T ))
        names(fs) <- gsub(params$meta_format, '', basename(fs))
        resolutions <- list()
        global_coords <- list()

        ## ome.tif is a bit unique in that a single file has information for all images
        f <- fs[[1]]
        meta <- RBioFormats::read.metadata(f)
        if(length(meta) != length(params$fov_names)){
          stop('The number of images not matching with the number of fov_names!')
        }
        full_res <<- meta[[1]]$coreMetadata[c('sizeX', 'sizeY')]
        n_channels <<- meta[[1]]$coreMetadata[c('sizeC')]
        ome <- RBioFormats::read.omexml(f, filter.metadata = TRUE)
        ome <- strsplit(ome, '<StageLabel Name=')[[1]][-1]
        if(length(ome) != length(params$fov_names)){
          stop('Wrong number of images found in metadata!')
        }
        names(ome) <- gsub('"', '', strsplit2(ome, ' ')[,1])
        tmp <- unlist(strsplit(tolower(ome), 'physicalsize'))
        tmp <- tmp[grepl('^x=|^y=', tolower(tmp))]
        tmp <- tmp[!duplicated(tmp)]
        tmp <- sort(tmp)
        pp <- as.numeric(gsub('\"| |x=|y=', '', tmp))
        channels <- unlist(strsplit(tolower(ome), 'channel id=|samplesperpixel='))
        channels <- channels[grepl('channel:', tolower(channels))]
        channels <- channels[order(channels)]
        channels <- as.character(sapply(channels, function(x) strsplit(x, ' name=')[[1]][2]))
        channels <- unique(as.character(gsub('\"| |"', '', channels)))
        fov_name_order <- sapply(names(ome), function(x) params$fov_names[grepl(x, params$fov_names)])

        resolutions <- list( # A bit redundant: generates a new resolutions every iteration
          'per_pixel_microns' = pp,
          'xydimensions_pixels' = as.numeric(full_res),
          'xydimensions_microns' = as.numeric(full_res) * pp,
          'channel_order' = channels,
          'fov_names' = fov_name_order, #Important for later
          'zslices' = meta[[1]]$coreMetadata$sizeZ
        )

        for(fov in params$fov_names){
          fsub <- fs[grepl(fov, names(fs))]
          names(fsub) <- gsub(paste0('_', fov, '$'), '', names(fsub))
          altnames <- as.character(as.integer(factor(names(fsub)))-1)
          omesub <- ome[which(sapply(names(ome), function(x) grepl(x, fov)))]
          tmp <- strsplit(tolower(omesub), 'position')[[1]]
          tmp <- suppressWarnings(na.omit(as.numeric(gsub('\"| |x=|y=|z=', '', tmp))))
          x_coord <- tmp[c(TRUE, FALSE, FALSE)]
          y_coord <- tmp[c(FALSE, TRUE, FALSE)]
          z_coord <- tmp[c(FALSE, FALSE, TRUE)]
          if( length(unique(x_coord))==1 ){
            x_coord <- unique(x_coord)
          }else{
            stop('More than one unique X tile value!')
          }
          if( length(unique(y_coord))==1 ){
            y_coord <- unique(y_coord)
          }else{
            stop('More than one unique Y tile value!')
          }
          if( length(unique(z_coord))==1 ){
            z_coord <- unique(z_coord)
          }
          fsub <- gsub(params$meta_format, params$im_format, fsub)
          if(!all(file.exists(fsub))){
            stop( paste0('Corresponding image for meta data file of ', fov, ' does not exist!'))
          }

          global_coord <- list()
          for( fsubidx  in 1:length(fsub) ){
            fsubi <- fsub[fsubidx]
            bit_name_full <- paste0(names(fsubi), '_', channels)
            bit_name_alias <- tolower( paste0(altnames[fsubidx], '_', channels) )
            for( zidx in 1:length(z_coord)){
              fsubidf <- suppressWarnings( data.frame(
                'x_microns' = x_coord,
                'y_microns' = y_coord,
                'z_microns' = z_coord[zidx],
                'image_file' = fsubi,
                'fov' = fov,
                'bit_name_full' = bit_name_full,
                'bit_name_alias' = bit_name_alias,
                'bit_name' = bit_name_alias
              ))
              global_coord[[length(global_coord) + 1]] <- fsubidf
            }
          }
          global_coord <- suppressWarnings( do.call(rbind, global_coord) )
          global_coords[[fov]] <- global_coord
        }
      }

    }

    global_coords <- suppressWarnings( data.frame(( do.call(rbind, global_coords)), check.names=F ) )
    check_colnames <- params$codebook_colnames[params$codebook_colnames!='gene']
    if(!all(check_colnames %in% global_coords$bit_name)){
      global_coords$bit_name <- global_coords$bit_name_full
      if(!all(check_colnames %in% global_coords$bit_name) & !params$inferred_codebookColnames){
        checkColBoolAlready = T
        warning('bit_name not matching codebookColumnNames: either edit GLOBALCOORD.csv or params$codebook_colnames (ignore this warning if codebook not needed, e.g. DAPI)')
      }
    }

    lapply(1:(length(resolutions)), function(idx){
      write(names(resolutions)[idx], paste0(params$out_dir, 'META_RESOLUTION', '.txt'), append = TRUE,  ncolumns=max(length(unlist(resolutions))))
      write(resolutions[[idx]], paste0(params$out_dir, 'META_RESOLUTION', '.txt'), append = TRUE,  ncolumns=max(length(unlist(resolutions))))
    })
    write.csv(global_coords, paste0(params$out_dir, 'GLOBALCOORD', '.csv'), row.names = F)
  }

  ## Load files
  resfs <- list.files(params$out_dir, pattern='META_RESOLUTION|GLOBALCOORD', full.names=TRUE)
  current_res_info <- readLines(resfs[grepl('META_RESOLUTION', resfs)])
  n <- current_res_info[c(TRUE, FALSE)]
  val <- current_res_info[c(FALSE, TRUE)]
  val <- lapply(val, function(v) strsplit(v, ' ')[[1]])
  resolutions <- setNames(val, n)
  global_coords <- data.table::fread(resfs[grepl('GLOBALCOORD', resfs)], data.table=FALSE)
  if(nrow(global_coords) < length(params$fov_names) * length(unique(global_coords$bit_names))){
    warning("GLOBALCOORD file seems too small - more rows expected.")
  }

  params$resolutions <<- resolutions
  params$global_coords <<- global_coords
  if( any(params$fov_names != params$resolutions$fov_names) ){
    warning('Updating params$fov_names')
    params$fov_names <- params$resolutions$fov_names
  }

  if( any(!(unique(params$global_coords$bit_name) %in% params$codebook_colnames[params$codebook_colnames!='gene'])) ){

    if( params$inferred_codebookColnames ){
      nChannels <- length(resolutions$channel_order)
      replace_channels <- setNames(
        resolutions$channel_order,
        paste0('CHANNEL', sprintf(paste0('%0', nchar(nChannels), 'd'), 1:nChannels)))
      for( cix  in 1:nChannels ){
        str_to_replace <- names(replace_channels)[cix]
        val_to_replace <- replace_channels[cix]
        params$codebook_colnames[grepl(paste0(str_to_replace, '$'), params$codebook_colnames)] <-
          gsub( paste0(str_to_replace, '$'), val_to_replace, params$codebook_colnames[grepl(paste0(str_to_replace, '$'), params$codebook_colnames)] )
      }
    }

    if( any(!(unique(params$global_coords$bit_name) %in% params$codebook_colnames[params$codebook_colnames!='gene'])) ){
      warning('Unsuccessful attempt at correcting placeholder codebook column names (params$codebook_colnames)')
      if( !checkColBoolAlready ){
        warning('bit_name not matching codebookColumnNames: either edit GLOBALCOORD.csv or params$codebook_colnames (ignore this warning if codebook not needed, e.g. DAPI)')
      }
    }else{
      warning('Placeholder codebook column names (params$codebook_colnames) found -- but successful at correcting ')
      colnames(params$raw_codebook) <- params$codebook_colnames[params$codebook_colnames!='gene']
    }
  }

  # return(params)
}

###

readImage <- function(
    fileName,
    nrows = NULL,
    ncols = NULL,
    endian = 'little',
    nbytes = 2,
    signed = T,
    normalise = T,
    ... ){

  ## check endian
  endian = tolower(endian)
  if(length(endian) > 1){
    warning('Only using the first endian value')
    endian <- endian[1]
  }
  if( !(endian %in% c('little', 'l', 'big', 'b')) ){
    stop('Invalid endian value: only accepts "little"/"l" or "big"/"b"')
  }
  if( endian %in% c('b', 'big') ){
    endian = 'big'
  }else{
    endian = 'little'
  }

  ## Check file exists
  file_size = file.size(fileName)
  if(is.na(file_size)){
    stop('File does not exist!')
  }
  if(!is.null(nrows) & !is.null(ncols)){
    nrows = as.numeric(nrows)
    ncols = as.numeric(ncols)
    file_size = nrows * ncols
    if(is.na(file_size)){
      stop('Invalid nrows and/or ncols!')
    }
  }

  ## Check if compatible with RBioFormats (read only 1 pixel for speed)
  t <- try(RBioFormats::read.image(fileName, subset = list(c=1, x=1, y=1), series=1), silent = T)
  if(inherits(t, 'try-error')){
    warning('Unable to read with RBioFormats')
    if( !grepl('.dax$', fileName) ){
      stop('Image type not supported!')
    }
    f <- file(fileName, 'rb')
    rawim <- readBin(fileName, endian=endian, what='integer',
                     size=nbytes, n=file_size, signed = signed) #16 bit little endian
    close(f)
    if(normalise){
      rawim <- rawim / (2^(nbytes * 8))
    }
    if(is.null(nrows)){
      ncols = nrows = sqrt(length(rawim)) # Assume square image
      if( nrows %% 1 != 0 ){
        stop('Not a square image: please provide nrows!')
      }
      warning('Assuming a square 2D image with a single channel')
    }
    if( nrows %% 1 != 0 ){
      stop('nrows needs to be an integer!')
    }
    im <- matrix(rawim, nrow=nrows, ncol=ncols)
  }else{
    im <- RBioFormats::read.image(fileName, normalize = normalise, ...)@.Data
    unloadNamespace('RBioFormats') #Stupid but fastest way to clear Java cache
  }
  return(im)
}


###

readImageList <- function(
    chosenFOV = NULL,
    fileNames = NULL,
    safeLoad = T,
    waitSeconds = 1,
    params = get('params', envir = globalenv()),
    ...
  ){

  ###
  if(is.factor(chosenFOV)){
    chosenFOV = as.character(chosenFOV)
  }
  if(!all(is.na(suppressWarnings(as.numeric(chosenFOV))))){
    chosenFOV = as.numeric(chosenFOV)
  }
  width <- height <- NULL
  if( !is.null(params$resolutions$xydimensions_pixels) ){
    width <- as.numeric(params$resolutions$xydimensions_pixels)[1]
    height <- as.numeric(params$resolutions$xydimensions_pixels)[2]
    if(is.na(width)){ width <- NULL }
    if(is.na(height)){ height <- NULL }
  }
  ###

  ###
  if(!exists('global_coords', envir = params)){
    stop('params$global_coords not found! Ensure readImageMetaData has been run!')
  }
  gcx <- params$global_coords
  if(is.null(fileNames)){
    if(is.null(chosenFOV)){
      stop('Provide either fileNames or chosenFOV!')
    }
    gcx <- gcx[gcx$fov %in% chosenFOV,]
  }else{
    gcx <- gcx[gcx$image_file %in% fileNames,]
  }
  fileNames <- unique(gcx$image_file)
  fovs <- unique(gcx$fov)
  ###

  ###
  imLists <- imList <- list()
  if(length(fovs)==1){
    message('\nReading images...')
    for(i in 1:length(fileNames)){

      message(paste0(i, ' of ', length(fileNames), '...'), appendLF = F)
      bitnames <- gcx[gcx$image_file==fileNames[i],'bit_name']
      bitindices <- params$global_coords[params$global_coords$image_file==fileNames[i], 'bit_name']
      bitindices <- which(bitindices %in% bitnames)

      if( params$im_format == '.ome.tif' ){
        fov_index <- which(params$fov_names %in% gcx[gcx$image_file==fileNames[i],'fov'])
        imListIntermediate <- list()
        if(safeLoad){
          for( j in bitindices ){
            imi <- readImage(
              fileNames[i],
              nrows = width,
              ncols = height,
              series = fov_index, #IMPORTANT FOR OME.TIF, otherwise will load ALL images
              subset = list(c = j), #Single channel
              ... )
            if( length(dim(imi)) != 2){
              stop('Problem loading single channel')
            }
            imListIntermediate[[ bitnames[j] ]] <- imi
          }
          imList <- c(imList, imListIntermediate)
        }else{
          imi <- readImage(
            fileNames[i],
            nrows = width,
            ncols = height,
            series = fov_index, #IMPORTANT FOR OME.TIF, otherwise will load ALL images
            ... )
          imList[[ basename(fileNames[i]) ]] <- imi
          Sys.sleep(waitSeconds)
        }
      }else{
        imi <- readImage(
          fileNames[i],
          nrows = width,
          ncols = height,
          ... )
        imList[bitnames[i]] <- imi
      }
    }
    imLists <- imList
  }else{
    for(j in 1:length(fovs)){
      # if(j > 1){ message('') } #Skip new line
      message(paste0('\nProcessing FOV ', fovs[j], ' (', j, ' of ', length(fovs), ')...'))
      fileNameSub <- unique(gcx$image_file[gcx$fov==fovs[j]])
      imLists[[fovs[j]]] <-
        readImageList(
          chosenFOV = NULL,
          fileNames = fileNameSub,
          safeLoad = safeLoad,
          params = params,
          waitSeconds = waitSeconds,
          ...
        )
    }
  }

  ###
  return(imLists)
}







