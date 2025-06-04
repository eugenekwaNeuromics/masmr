## Functions for prepping environment.

require(reticulate)
require(data.table)

strsplit2 <- function(x, split, ...){
  # From limma
  x <- as.character(x)
  n <- length(x)
  s <- strsplit(x, split = split, ...)
  nc <- unlist(lapply(s, length))
  out <- matrix("", n, max(nc))
  for (i in 1:n) {
    if (nc[i])
      out[i, 1:nc[i]] <- s[[i]]
  }
  return( out )
}

###

establishParams <- function(
    ## User to fill in
    imageDirs,
    imageFormat,
    outDir,
    codebookFileName,
    dapiDir,

    ## Will attempt automatic fill
    metaFormat,
    dapiFormat,
    pythonLocation,
    codebookColumnNames,

    ## Automatically filled in
    seed = 12345,
    verbose = T,
    resumeMode = T,

    ## Optional
    fpkmFileName = NULL,
    nProbesFileName = NULL
){

  ## Create a environment to house parameters
  ## Ignore DAPI for now
  ## <<- to ensure global environment gets this too
  params <<- new.env()
  params$seed <<- as.integer( seed )

  ## Specify imageDir
  im_dir = imageDirs
  if(length(imageDirs)==1){
    im_dir = imageDirs

    if( length(list.dirs(im_dir)) > 1){

      ## DEFAULT: Just stay in the specified directory
      ## ALT 1: There are no other images in the specified directory,
      ## Assume immediate daughter folders are the focus
      if(
        length( list.files(im_dir, pattern = imageFormat, recursive = F) )==0 &
        length( list.files(im_dir, pattern = imageFormat, recursive = T) )>0
      ){
        im_dir = list.dirs(im_dir)[-1]
      }

      if(
        length( list.files(im_dir, pattern = imageFormat, recursive = F) )==0 &
        length( list.files(im_dir, pattern = imageFormat, recursive = T) )==0
      ){
        stop('Invalid image directory OR invalid image format ( ',
             imageFormat,
             ' ) provided!')
      }
    }

  }
  if(length(im_dir)==0){
    stop('Invalid image directory provided!')
  }
  params$im_dir <<- im_dir

  ## Check if image format exists in all the folders specified
  if( missing(metaFormat) ){
    availformats <- list.files(im_dir)
    availformats <- table(sapply(availformats, function(x) gsub(strsplit2(x, '[.]')[,1], '', x) ))
    availformats <- names(availformats)[order(availformats, decreasing = T)]

    if(missing(metaFormat)){
      if(length(availformats)==1){
        metaFormat <- availformats
      }
      if('.xml' %in% availformats & imageFormat=='.dax'){
        metaFormat <- '.xml'
      }
      if(imageFormat=='.ome.tif'){
        metaFormat <- '.ome.tif'
      }
      message( paste0('Assuming ', metaFormat, ' is extension for meta files...') )
    }
  }
  if( !all(sapply(params$im_dir, function(x) any(grepl( paste0(imageFormat, '$'), list.files(x) )))) ){
    stop('Invalid image format provided!')
  }
  params$im_format <<- imageFormat

  ## Check if meta format exists in all the folders specified
  if( !all(sapply(params$im_dir, function(x) any(grepl( paste0(metaFormat, '$'), list.files(x) )))) ){
    ## Check if xml provided for .dax
    if(imageFormat=='.dax' & metaFormat!='.xml'){
      warning(
        paste0('Metadata format for .dax is likely to be .xml, but is currently ',
               metaFormat)
      )
    }
    stop('Invalid meta data format provided!')
  }
  params$meta_format <<- metaFormat


  if(!missing(dapiDir)){
    ## Check the DAPI directory
    if(missing(dapiFormat)){
      dapiFormat <- imageFormat # assume same as image format
    }
    ## Check that dapiformat exists in dapi directory
    if( !all(sapply(dapiDir, function(x) any(grepl( paste0(dapiFormat, '$'), list.files(x) )))) ){
      stop('Invalid DAPI image format provided!')
    }
    params$dapi_dir <- dapiDir
    params$dapi_format <- dapiFormat
  }



  ## Get directory names
  dir_names <- basename(im_dir)
  all_fov_names <- setNames(lapply(im_dir, function(f){
    fname <- list.files(f, full.names=FALSE, pattern=params$im_format)
    fname <- gsub(paste(c( paste0(dir_names, '_'), params$im_format), collapse='|'), '', fname)
  }), dir_names)
  fov_names <- Reduce(intersect, all_fov_names)
  if(!all(unlist(lapply(all_fov_names, function(x) all(x %in% fov_names))))){
    stop('Missing image!')
  }
  if(length(fov_names)==length(unlist(all_fov_names))){
    ## Assume this is a mistake: file names should be duplicated across channels
    warning('Unclear file organisation: assuming file name specified by the final underscore')
    fov_names <- unlist( sapply(fov_names, function(x){
      y <- strsplit(x, '_')[[1]]
      y <- y[length(y)]
      return(y)
    }) )
    fov_names <- sort(unique(fov_names))
  }
  params$dir_names <<- dir_names
  params$fov_names <<- fov_names


  ## Check codebookFileName
  if(!file.exists(codebookFileName)){
    stop('Codebook file name does not exist!')
  }
  params$codebook_file <<- codebookFileName


  ## Check that number of columns matches
  codebook <- data.table::fread(params$codebook_file, header = F, data.table = F)
  ngenes <- nrow(codebook)
  skiprows = 0
  while( !all(colnames(codebook)==paste0('V', 1:ncol(codebook))) ){
    # Check that first row is not a title
    skiprows = skiprows + 1
    if(skiprows==ngenes){
      stop('Uncertain which row codebook begins!')
    }
    codebook <- data.table::fread(params$codebook_file, skip = skiprows, data.table = F)
  }
  params$inferred_codebookColnames <<- F
  if( missing(codebookColumnNames) ){
    warning('codebookColumnNames not provided : will provide placeholders...')
    nChannels <- floor(ncol(codebook) / length(im_dir))
    codebookColumnNames <- paste0( rep( basename(im_dir), each=nChannels ), '_CHANNEL', rep(sprintf(paste0('%0', nchar(nChannels), 'd'), 1:nChannels), length(im_dir)))
    params$inferred_codebookColnames <<- T
  }
  if(ncol(codebook) == length(codebookColumnNames)+1){
    ## Figure out which column is the gene column
    newCodebookColNames <- rep('gene', ncol(codebook))
    isgene <- rep(F, ncol(codebook))
    for(i in 1:ncol(codebook)){
      isgene[i] <- all(suppressWarnings(is.na(as.numeric(codebook[,i]))))
    }
    if(sum(isgene) != 1){
      stop('Uncertain which column in codebook contains gene!')
    }
    newCodebookColNames[!isgene] <- codebookColumnNames
    codebookColumnNames <- newCodebookColNames
  }
  if(ncol(codebook) != length(codebookColumnNames)){
    stop('Number of codebook column names do not match number of columns in codebook!')
  }
  params$codebook_colnames <- codebookColumnNames
  colnames(codebook) <- codebookColumnNames
  rownames(codebook) <- codebook$gene
  codebook$gene <- NULL
  params$raw_codebook <<- codebook
  nbits <- unique(rowSums(codebook))
  if(length(nbits) > 1){
    warning('Unusual bit design')
  }
  params$nbits <<- nbits


  ## Make sure that resumeMode is a boolean
  if(is.logical(resumeMode[1])){
    params$resumeMode <<- resumeMode[1]
  }else{
    params$resumeMode <<- F
  }
  
  ## Make sure that resumeMode is a boolean
  if(is.logical(verbose[1])){
    params$verbose <<- verbose[1]
  }else{
    params$verbose <<- T
  }


  ## Create out directory
  if(!dir.exists(outDir)){
    dir.create(outDir)
  }
  params$parent_out_dir <<- outDir
  params$out_dir <<- outDir


  ## Python environment
  if(missing(pythonLocation)){
    pythonLocation = reticulate::miniconda_path()
  }
  params$python_location <<- pythonLocation


  ## Optional params
  if(!is.null(fpkmFileName)){
    if(file.exists(fpkmFileName)){
      params$fpkm_file = fpkmFileName
    }
  }
  if(!is.null(nProbesFileName)){
    if(file.exists(nProbesFileName)){
      params$n_probes_file = nProbesFileName
    }
  }

  if(params$verbose){

    message(
      paste0(
        'resumeMode = ', params$resumeMode, '...'
      )
    )
    message(
      paste0(
        length(params$fov_names), ' images across ', ncol(params$raw_codebook),
        ' cycles ( ', paste(params$nbits, collapse=' & '), ' ON-bit codes )...'
      )
    )
  }
}

##

establishCleanSlate <- function(){
  cleanSlate <- ls(all.names = T, envir = globalenv())
  cleanSlate <<- c(cleanSlate, 'cleanSlate')
  return(cleanSlate) #Print values in cleanSlate
}

##

cleanUp <- function( cleanSlate=get('cleanSlate', envir = globalenv()) ){
  cleanSlate <- unlist(cleanSlate)
  forRemoval <- ls( envir = globalenv() )
  forRemoval <- forRemoval[!(forRemoval %in% cleanSlate)]
  rm(list=forRemoval, envir = globalenv())
}

##

checkReadyLoop <- function(){

  if(!exists('params', envir=globalenv())){
    warning(
    'No params environment detected:
    Either need to run establishParams(); or
    ensure params environment specified when needed for all downstream functions')
    return(F)
  }

  params = get('params', envir = globalenv())
  if(!exists('resolutions', envir=params)){
    warning(
      'params$resolutions not found:
    Ensure readImageMetaData() has been run!')
    return(F)
  }
  if(!exists('ordered_codebook', envir=params)){
    warning(
    'params$ordered_codebook not found:
    Ensure prepareCodebook() has been run!')
    return(F)
  }
  if(!exists('current_fov', envir = params)){
    warning(
      'params$current_fov not found:
      Ensure that fovNames are specified as needed downstream')
    return(F)
  }
  if(!exists('brightness_min', envir = params)){
    warning(
      'params$brightness_min not found:
      If desired, ensure getAnchorParams() has been run')
    return(F)
  }
  if(!exists('brightness_max', envir = params)){
    warning(
      'params$brightness_max not found:
      If desired, ensure getAnchorParams() has been run')
    return(F)
  }

  return(T)
}
