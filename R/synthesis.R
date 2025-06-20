## Function for synthesising data

require(data.table)

synthesiseData <- function(
    referenceFOV,
    stitchChosenColumn = 1,
    nLeadingZeroes = 5,
    cellOverlapFraction = 0.25,
    gZip = TRUE,
    removeRedundancy = TRUE,
    subsetFOV = NULL,
    customOutDirName = NULL,
    params = get('params', envir = globalenv())
){

  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose
  }else{
    verbose = T
  }

  ## First, check that append=True is working for data.table
  test_file <- paste0(params$parent_out_dir, 'DATATABLE_TEST_FILE.csv', ifelse(gZip, '.gz', ''))
  data.table::fwrite(data.frame('x'=1), test_file, row.names = F, showProgress = F)
  data.table::fwrite(data.frame('x'=1), test_file, append = T, row.names = F, showProgress = F)
  test_df <- data.table::fread(test_file, data.table = F, showProgress = F)
  if( nrow(test_df) != 2 ){
    warning('Unable to append to .csv.gz with data.table! Attempting to update data.table package...')
    oldversion = utils::packageVersion('data.table')
    detach(package:data.table, unload = T)
    install.packages("data.table")
    newversion = utils::packageVersion('data.table')
    if(newversion==oldversion){
      stop(paste0(
        '\nUnsuccessful update: try updating manually, using an older version of data.table (1.14.8), or setting gZip = FALSE',
        ifelse(newversion=='1.17.0',
        '\nVersion 1.17.0 detected: there is a known bug with this version (https://github.com/Rdatatable/data.table/issues/6863)',
        '')
        ))
    }
    ## Check again
    data.table::fwrite(data.frame('x'=1), test_file, row.names = F, showProgress =F)
    data.table::fwrite(data.frame('x'=1), test_file, append = T, row.names = F, showProgress =F)
    test_df <- data.table::fread(test_file, data.table = F, showProgress =F)
    if(nrow(test_df)==2){
      if(verbose){ message( paste0('Successful update of data.table: ', oldversion, ' --> ', newversion, '...') ) }
    }else{
      stop('Successful update but unable to correct problem: try using an older version of data.table (1.14.8), or setting gZip = FALSE')
    }
  }
  file.remove(test_file)
  rm(test_df, test_file)

  ## Read information
  if(is.null(customOutDirName)){
    params$out_dir <<- gsub('[/][/]', '/', paste0(params$parent_out_dir, '/OUT/'))
  }else{
    params$out_dir <<- gsub('[/][/]', '/', paste0(params$parent_out_dir, '/', customOutDirName, '/'))
  }
  if(!dir.exists(params$out_dir)){
    dir.create(params$out_dir, recursive = T)
  }
  if(is.null(subsetFOV)){
    fileChecks <- c(
      'OUT_GLOBALCOORD.csv',
      paste0('OUT_SPOTCALL_PIXELS.csv', ifelse(gZip, '.gz', '')),
      paste0('OUT_CELLSEG_PIXELS.csv', ifelse(gZip, '.gz', '')),
      paste0('OUT_CELLEXPRESSION.csv', ifelse(gZip, '.gz', '')),
      paste0('OUT_CELLS.csv', ifelse(gZip, '.gz', ''))
      )
    fileChecks <- fileChecks[fileChecks %in% list.files(params$out_dir)]
    if( length(fileChecks) > 0 ){
      warning('Overwriting existing files')
      t <- suppressWarnings(try(file.remove( paste0(params$out_dir, fileChecks) )))
      if(inherits(t, 'try-error')){
        stop('Unable to delete existing files!')
      }
    }
  }

  ## Load files
  fs <- list.files( params$parent_out_dir, full.names = T, recursive = T, all.files = T)
  fs <- fs[!grepl(params$out_dir, fs)]

  ## Read the first GLOBALCOORD and META_RESOLUTION file
  global_coords <- try( data.table::fread(fs[grepl('GLOBALCOORD.csv', fs)][1], data.table = F) )
  if( inherits(global_coords, 'try-error') ){
    stop('Unable to load a GLOBALCOORD.csv file!')
  }
  global_coords <- global_coords[!duplicated(global_coords$fov),]
  current_res_info <- try( readLines(fs[grepl('META_RESOLUTION.txt', fs)][1]) )
  if( inherits(current_res_info, 'try-error') ){
    stop('Unable to load a META_RESOLUTION file!')
  }
  n <- current_res_info[c(TRUE, FALSE)]
  val <- current_res_info[c(FALSE, TRUE)]
  val <- lapply(val, function(v) strsplit(v, ' ')[[1]])
  resolutions <- setNames(val, n)
  fov_names <- resolutions$fov_names
  params$all_fovs <<- fov_names

  ## Subset FOVs...
  if( !is.null(subsetFOV) ){
    if(verbose){ message('Subsetting FOVs...') }
    if( is.logical(subsetFOV) ){
      if( length(subsetFOV) != length(fov_names) ){
        stop("If providing a boolean vector for subsetFOV, it must be the same length as params$all_fovs!")
      }
      subsetFOV <- fov_names[subsetFOV]
    }
    subsetFOV <- as.character(subsetFOV)
    if( any(!(subsetFOV %in% fov_names)) ){
      stop( 'Your subsetFOV contains FOVs not found in your dataset!' )
    }
    fov_names <- fov_names[fov_names %in% subsetFOV]
  }

  if(length(fov_names)<=1){
    stop('Insufficient FOVs to attempt data synthesis! Ensure all FOVs have been processed, or that subsets have >1 FOV!')
  }

  params$fov_names <<- fov_names
  params$global_coords <<- global_coords[global_coords$fov %in% fov_names,]
  params$resolutions <<- resolutions

  ## Read codebook for gene names
  codebook <- try(data.table::fread( fs[grepl('ORDEREDCODEBOOK.csv', fs)][1], data.table = F))
  if( inherits(codebook, 'try-error') ){
    stop('Unable to load ORDEREDCODEBOOK file!')
  }
  if( all(rownames(params$raw_codebook) %in% codebook[,1]) ){
    rownames(codebook) <- codebook[,1]
    codebook[,1] <- NULL
  }
  g <- rownames(codebook)

  ## Prepare file subsets
  fovfs <- fs[Reduce('+', lapply(fov_names, function(x) grepl(x, fs))) > 0]
  spotcallfs <- fovfs[grepl('^SPOTCALL_', basename(fovfs))]
  cellsegfs <- fovfs[grepl('^CELLSEG_', basename(fovfs))]
  stitchfs <- fovfs[grepl('^STITCH_', basename(fovfs))]
  if( length(spotcallfs)== 0 ){
    stop('Unable to find SPOTCALL_{FOV}.csv.gz files!')
  }
  if( length(cellsegfs)== 0 ){
    stop('Unable to find CELLSEG_{FOV}.csv.gz files!')
  }
  if( length(spotcallfs)== 0 ){
    stop('Unable to find STITCH_{FOV}.csv files!')
  }

  ## Pick a reference FOV
  if(missing(referenceFOV)){
    if(verbose){ message(paste0('\nReference FOV not specified...Choosing FOV with most spots (approximated from file size)...')) }
    # nspots <- lapply(spotcallfs, function(fx){
    #   nrow(data.table::fread(fx))
    # })
    nspots <- lapply(spotcallfs, function(x) file.info(x)$size) #Faster
    referenceFOV <- basename( spotcallfs[which.max(unlist(nspots))] )
    referenceFOV <- gsub('SPOTCALL_|[.]csv|[.]gz', '', referenceFOV)
  }
  if( !(referenceFOV %in% fov_names) ){
    stop( paste0('FOV chosen (', referenceFOV, ' ) not in fov_names! Please check params$fov_names!'))
  }
  if(verbose){ message( paste0('FOV chosen: ', referenceFOV, '...') ) }

  if( !is.null(subsetFOV) ){
    ## Update out_dir name
    if(is.null(customOutDirName)){
      params$out_dir <<- gsub('[/][/]', '/', paste0(params$parent_out_dir, '/OUT/SUBSET_', referenceFOV, '/'))
    }else{
      params$out_dir <<- gsub('[/][/]', '/', paste0(params$parent_out_dir, '/', customOutDirName, '/'))
    }
    if(!dir.exists(params$out_dir)){
      dir.create(params$out_dir, recursive = T)
    }
    fileChecks <- c(
      'OUT_GLOBALCOORD.csv',
      paste0('OUT_SPOTCALL_PIXELS.csv', ifelse(gZip, '.gz', '')),
      paste0('OUT_CELLSEG_PIXELS.csv', ifelse(gZip, '.gz', '')),
      paste0('OUT_CELLEXPRESSION.csv', ifelse(gZip, '.gz', '')),
      paste0('OUT_CELLS.csv', ifelse(gZip, '.gz', ''))
    )
    fileChecks <- fileChecks[fileChecks %in% list.files(params$out_dir)]
    if( length(fileChecks) > 0 ){
      warning('Overwriting existing files')
      t <- suppressWarnings(try(file.remove( paste0(params$out_dir, fileChecks) )))
      if(inherits(t, 'try-error')){
        stop('Unable to delete existing files!')
      }
    }
  }

  ## Update global_coords
  if(verbose){ message('\nUpdating global coordinates...') }
  if( length(stitchChosenColumn) > 1 ){
    warning('More than one stitch dataframe column specified: taking mean of vectors...')
  }
  gcx <- global_coords[,c('x_microns', 'y_microns', 'z_microns', 'fov')]
  stitchResults <- setNames(lapply(stitchfs, function(fx){
    dfx <- data.table::fread(fx, data.table = F)
    newdf <- NULL
    if( all(is.character(stitchChosenColumn)) ){
      if( all(stitchChosenColumn %in% colnames(dfx)) ){
        if(length(stitchChosenColumn)==1){
          shiftx <- as.complex(dfx[,stitchChosenColumn])
          newdf <- data.frame('fov' = dfx$fov, 'shift_pixel' = shiftx)
        }else{
          shiftdf <- dfx[,stitchChosenColumn]
          shiftx <- rep(0+0i, nrow(shiftdf))
          for(j in 1:ncol(shiftdf)){
            shiftx <- shiftx + as.complex(shiftdf[,j])
          }
          shiftx <- shiftx / ncol(shiftdf)
          newdf <- data.frame('fov' = dfx$fov, 'shift_pixel' = shiftx)
        }
      }else{
        stop('stitchChosenColumn is a character, but does not match any columns in the stitch dataframe!')
      }
    }
    if( all(is.numeric(stitchChosenColumn)) ){
      if( length(stitchChosenColumn)== 1){
        shiftx <- as.complex(dfx[,colnames(dfx)[colnames(dfx)!='fov'][stitchChosenColumn]])
        newdf <- data.frame('fov' = dfx$fov, 'shift_pixel' = shiftx)
      }else{
        shiftdf <- dfx[,colnames(dfx)[colnames(dfx)!='fov'][stitchChosenColumn]]
        shiftx <- rep(0+0i, nrow(shiftdf))
        for(j in 1:ncol(shiftdf)){
          shiftx <- shiftx + as.complex(shiftdf[,j])
        }
        shiftx <- shiftx / ncol(shiftdf)
        newdf <- data.frame('fov' = dfx$fov, 'shift_pixel' = shiftx )
      }
    }
    if(is.null(newdf)){
      stop('Unable to parse desired column in stitch dataframe to use!')
    }
    if( any(is.na(newdf)) ){
      stop('NA stitch vectors returned!')
    }
    return(newdf)
  }), stitchfs)

  ## Propagate new coordinates wrt reference FOV
  if(verbose){ message('\nPropagating new coordinates...') }
  fovs_processed <- unique(do.call(rbind, stitchResults)$fov)
  new_gcx <- list()
  new_gcx[[referenceFOV]] <- sum(gcx[gcx$fov==referenceFOV,c('x_microns', 'y_microns')] * c(1, 1i))
  STOP = F
  problemCounter = 0
  while(length(new_gcx)<length(fovs_processed) & !STOP){
    current_length = length(new_gcx)
    subStitchResults <- stitchResults[ grepl( paste(paste0(names(new_gcx),'.csv'), collapse='|'), names(stitchResults)) ]
    for(j in 1:length(subStitchResults)){
      ref_cx <- sapply(names(new_gcx), function(x) grepl(paste0(x, '.csv'), names(subStitchResults)[j]) )
      if( !is.logical(ref_cx) ){ next }
      if( all(!ref_cx) ){ next }
      ref_cx <- new_gcx[[names(ref_cx)[ref_cx]]]
      df <- subStitchResults[[j]]
      df <- df[!(df$fov %in% names(new_gcx)),]
      if(nrow(df)==0){next}
      df$shift_micron <-
        Re(as.complex(df$shift_pixel)) * as.numeric(resolutions$per_pixel_microns)[1] +
        1i * Im(as.complex(df$shift_pixel)) * as.numeric(resolutions$per_pixel_microns)[2]
      df$shift_micron <- df$shift_micron + ref_cx
      for( rowi in 1:nrow(df) ){
        if(is.null(new_gcx[[ df$fov[rowi] ]])){
          new_gcx[[ df$fov[rowi] ]] <- df$shift_micron[rowi]
        }else{
          new_gcx[[ df$fov[rowi] ]] <- mean( c(df$shift_micron[rowi], new_gcx[[ df$fov[rowi] ]]) )
        }
      }
    }
    new_length = length(new_gcx)
    if( length(new_gcx)<length(fovs_processed) & (new_length == current_length) ){
      problemCounter = problemCounter + 1 #Technically, this needs only fail once, but to be safe, let it loop 5 times
      if( problemCounter > 5 ){
        STOP = T
      }
    }
  }

  if(length(new_gcx) < length(fovs_processed)){
    stop(
    '
    Stitched coordinates cannot be propagated to some FOVs!
    Ensure that FOVs are adjacent in your STITCH folder.
    If you have disjoint FOVs (e.g. multiple samples),
    use the subsetFOVs field to restrict to just FOVs belonging to each intact sample.
    ')
  }

  gcx <- gcx[gcx$fov %in% names(new_gcx),]
  if(nrow(gcx)==0){
    stop('FOV names not consistently preserved!')
  }
  colnames(gcx)[colnames(gcx) %in% c('x_microns', 'y_microns')] <- paste0(colnames(gcx)[colnames(gcx) %in% c('x_microns', 'y_microns')], '_raw')
  gcx$x_microns <- Re( unlist(new_gcx)[match(gcx$fov, names(new_gcx))] )
  gcx$y_microns <- Im( unlist(new_gcx)[match(gcx$fov, names(new_gcx))] )
  params$global_coords_final <<- gcx
  data.table::fwrite(gcx, paste0(params$out_dir, 'OUT_GLOBALCOORD.csv'), row.names = F, showProgress = F)

  ## For redundancy, create a winLose temp file
  winLoseFile <- paste0( params$out_dir, 'TEMP_WINLOSE.csv' )
  if( file.exists(winLoseFile) ){
    warning( paste0(
      'Found an existing temporary file: ', winLoseFile,
      '...\nIf intending a fresh run this should be deleted, otherwise will read from this file...' )
      )
  }
  cellDropFile <- paste0( params$out_dir, 'TEMP_CELLDROP.csv' )
  if( file.exists(cellDropFile) ){
    warning( paste0(
      'Found an existing temporary file: ', cellDropFile,
      '...\nIf intending a fresh run this should be deleted, otherwise will read from this file...' )
    )
  }

  ## Begin looping through files to synthesise data
  if(verbose){ message('\nSynthesising data...') }
  finalcellseg <- finaldf <- data.frame()
  NCHAR_NAME = nLeadingZeroes #Assume a max of 99999 cells per FOV
  for( fovName in fov_names ){

    if(verbose){
      message('')
      message( paste0(which(fov_names == fovName), ' of ', length(fov_names), '...'), appendLF = F )
    }
    spotcallfx <- spotcallfs[grepl( paste0(fovName, '.csv'), spotcallfs )]
    if(length(spotcallfx)==0){
      if(verbose){ message('No spotcall file found...Skipping...', appendLF = F) }
      next
    }
    if(length(spotcallfx)>1){
      stop('More than one spotcall file specified!')
    }
    isThereStitch <- any(grepl(paste0('STITCH_', fovName, '.csv'), names(stitchResults)))
    if( !isThereStitch ){
      if(verbose){ message('\nNo stitch file found...Assuming spots called in this FOV are spurious...Skipping...', appendLF = F) }
      next
    }

    ## Update coordinates
    spotcalldf <- try(data.table::fread(spotcallfx, data.table = F, showProgress = F))
    if(inherits(spotcalldf, 'try-error')){
      stop('Unable to read spotcall file!')
    }
    coords <-
      spotcalldf$WX * as.numeric(resolutions$per_pixel_microns[1])  +
      1i * spotcalldf$WY * as.numeric(resolutions$per_pixel_microns[2]) +
      gcx[match(spotcalldf$fov, gcx$fov), 'x_microns'] +
      1i * gcx[match(spotcalldf$fov, gcx$fov), 'y_microns']
    spotcalldf$Xm <- Re(coords)
    spotcalldf$Ym <- Im(coords)
    spotcalldf$IDX <- spotcalldf$WX + (spotcalldf$WY-1) * as.numeric(resolutions$xydimensions_pixels[1])

    ## Load neighbours
    if(removeRedundancy){

      if(verbose){ message('\nRemoving redundant spot calls...', appendLF = F) }
      fovNeighbours <- stitchResults[grepl(paste0('STITCH_', fovName, '.csv'), names(stitchResults))][[1]]
      colnames(fovNeighbours) <- tolower(colnames(fovNeighbours))
      if( sum(colnames(fovNeighbours)=='fov') != 1){
        stop('STITCH files need to have a single column named "fov"!')
      }
      fovNeighbours <- as.vector( fovNeighbours[,'fov'] )

      ## But also get the first order neighbours for each neighbour in fovNeighbours
      firstDegreeNeighbours <- unlist( lapply( fovNeighbours, function(x){
        tmpdf <- stitchResults[grepl(paste0('STITCH_', x, '.csv'), names(stitchResults))][[1]]
        colnames(tmpdf) <- tolower(colnames(tmpdf))
        if( sum(colnames(tmpdf)=='fov') != 1){
          stop('STITCH files need to have a single column named "fov"!')
        }
        tmpdf <- as.vector( tmpdf[,'fov'] )
        tmpdf <- tmpdf[tmpdf!=fovName]
        return(tmpdf)
      }) )
      fovNeighbours <- c(fovNeighbours, firstDegreeNeighbours)
      fovNeighbours <- unique(fovNeighbours)

      redundancyWinnersLosers <- NULL
      if(file.exists(winLoseFile)){
        redundancyWinnersLosers <- data.table::fread(winLoseFile, data.table = F, showProgress = F)
      }

      ## Not strict enough
      # ## NEWER ##
      # ## Load all data
      # allNeighbours <- list()
      # for( i in 1:length(fovNeighbours) ){
      #   neighbourfx <- spotcallfs[grepl( paste0(fovNeighbours[i], '.csv'), spotcallfs )]
      #   if( length(neighbourfx)==0 ){ next }
      #   neighbourdf <- try(data.table::fread(neighbourfx, data.table = F))
      #   if(inherits(neighbourdf, 'try-error')){
      #     stop('Unable to read neighbouring spotcall file!')
      #   }
      #
      #   neighbour_coords <-
      #     neighbourdf$WX * as.numeric(resolutions$per_pixel_microns[1])  +
      #     1i * neighbourdf$WY * as.numeric(resolutions$per_pixel_microns[2]) +
      #     gcx[match(neighbourdf$fov, gcx$fov), 'x_microns'] +
      #     1i * gcx[match(neighbourdf$fov, gcx$fov), 'y_microns']
      #   neighbourdf$Xm <- Re(neighbour_coords)
      #   neighbourdf$Ym <- Im(neighbour_coords)
      #
      #   neighbourdf <- neighbourdf[
      #     neighbourdf$Xm >= min(spotcalldf$Xm)
      #     & neighbourdf$Xm < max(spotcalldf$Xm)
      #     & neighbourdf$Ym >= min(spotcalldf$Ym)
      #     & neighbourdf$Ym < max(spotcalldf$Ym),
      #   ]
      #   if( nrow(neighbourdf) == 0 ){ next }
      #   allNeighbours[[ length(allNeighbours) + 1 ]] <- neighbourdf
      # }
      # ## Now split overlaps even further
      # for(i in 1:length(allNeighbours) ){
      #   ref <- allNeighbours[[i]]
      #   everyoneElse <- c(1:length(allNeighbours))[-i]
      #   if(length(everyoneElse)==0){ next }
      #   for(j in everyoneElse){
      #     quer <- allNeighbours[[j]]
      #     querbool <- (
      #       quer$Xm >= min(ref$Xm)
      #       & quer$Xm < max(ref$Xm)
      #       & quer$Ym >= min(ref$Ym)
      #       & quer$Ym < max(ref$Ym)
      #     )
      #     quer <- quer[querbool,]
      #     if(nrow(quer)==0){ next }
      #     boundingWindow <- c(range(quer$Xm), range(quer$Ym) )
      #     underContention <-
      #       (ref$Xm >= boundingWindow[1]
      #        & ref$Xm <  boundingWindow[2]
      #        & ref$Ym >= boundingWindow[3]
      #        & ref$Ym <  boundingWindow[4])
      #
      #     ## If there are more spots in window from neighbour, filter out the reference FOV
      #     if(sum(underContention) < ( nrow(quer) )){
      #       allNeighbours[[length(allNeighbours) + 1]] <- quer
      #       allNeighbours[[i]] <- ref[!underContention,]
      #       allNeighbours[[j]] <- allNeighbours[[j]][!querbool,]
      #     }else{
      #       allNeighbours[[length(allNeighbours) + 1]] <- ref[underContention,]
      #       allNeighbours[[i]] <- ref[!underContention,]
      #       allNeighbours[[j]] <- allNeighbours[[j]][!querbool,]
      #     }
      #   }
      # }
      # ## Now compare against our spotcalldf
      # for(i in 1:length(allNeighbours)){
      #   neighbourdf <- allNeighbours[[i]]
      #   if(nrow(neighbourdf)==0){ next }
      #   boundingWindow <- c(range(neighbourdf$Xm), range(neighbourdf$Ym) )
      #   underContention <-
      #     (spotcalldf$Xm >= boundingWindow[1]
      #      & spotcalldf$Xm <  boundingWindow[2]
      #      & spotcalldf$Ym >= boundingWindow[3]
      #      & spotcalldf$Ym <  boundingWindow[4])
      #   ## If there are more spots in window from neighbour, filter out the current FOV
      #   if(sum(underContention) < ( nrow(neighbourdf) )){
      #     spotcalldf <- spotcalldf[!underContention,]
      #   }
      # }
      # ## NEWER ##


      ## OLDER ##
      for( i in 1:length(fovNeighbours) ){

        ## First, check if fovNeighbours[i] is already competed against
        checkCompetition <- NULL
        if(!is.null(redundancyWinnersLosers)){
          checkCompetition <- redundancyWinnersLosers[
            redundancyWinnersLosers$reference==fovNeighbours[i]
            & redundancyWinnersLosers$query==fovName,]
        }

        neighbourfx <- spotcallfs[grepl( paste0(fovNeighbours[i], '.csv'), spotcallfs )]
        if( length(neighbourfx)==0 ){ next }
        neighbourdf <- try(data.table::fread(neighbourfx, data.table = F, showProgress =F))
        if(inherits(neighbourdf, 'try-error')){
          stop('Unable to read neighbouring spotcall file!')
        }

        ## Instead of pixel wise-info, we will use micron info here [both should work]
        neighbour_coords <-
          neighbourdf$WX * as.numeric(resolutions$per_pixel_microns[1])  +
          1i * neighbourdf$WY * as.numeric(resolutions$per_pixel_microns[2]) +
          gcx[match(neighbourdf$fov, gcx$fov), 'x_microns'] +
          1i * gcx[match(neighbourdf$fov, gcx$fov), 'y_microns']
        neighbourdf$Xm <- Re(neighbour_coords)
        neighbourdf$Ym <- Im(neighbour_coords)

        neighbourdf <- neighbourdf[
          neighbourdf$Xm >= min(spotcalldf$Xm)
          & neighbourdf$Xm <= max(spotcalldf$Xm)
          & neighbourdf$Ym >= min(spotcalldf$Ym)
          & neighbourdf$Ym <= max(spotcalldf$Ym),
        ]

        competedBefore = F
        if(!is.null(checkCompetition)){
          if( nrow(checkCompetition) == 1 ){
            competedBefore = T
            if( !as.logical(checkCompetition$refWon) ){
              previousBound <- c(
                checkCompetition$queryBound1,
                checkCompetition$queryBound2,
                checkCompetition$queryBound3,
                checkCompetition$queryBound4
              )
              neighbourdf <- neighbourdf[
                !(neighbourdf$Xm >= previousBound[1]
                & neighbourdf$Xm <= previousBound[2]
                & neighbourdf$Ym >= previousBound[3]
                & neighbourdf$Ym <= previousBound[4]),
              ]
            }
          }
        }
        if( nrow(neighbourdf) == 0 ){ next }
        boundingWindow <- c(range(neighbourdf$Xm), range(neighbourdf$Ym) )
        underContention <-
          (spotcalldf$Xm >= boundingWindow[1]
           & spotcalldf$Xm <=  boundingWindow[2]
           & spotcalldf$Ym >= boundingWindow[3]
           & spotcalldf$Ym <=  boundingWindow[4])

        ## NEW ## Actually, this might shrink the bounding window, which we don't want
        ## Double check that bounding window is correctly sized
        # newBounding <- c(
        #   range(spotcalldf$Xm[underContention]), range(spotcalldf$Ym[underContention])
        # )
        # neighbourdf <- neighbourdf[
        #   neighbourdf$Xm >= newBounding[1]
        #   & neighbourdf$Xm <= newBounding[2]
        #   & neighbourdf$Ym >= newBounding[3]
        #   & neighbourdf$Ym <= newBounding[4],
        # ]
        # if( nrow(neighbourdf) == 0 ){ next }
        # boundingWindow <- c(range(neighbourdf$Xm), range(neighbourdf$Ym) )
        # underContention <-
        #   (spotcalldf$Xm >= boundingWindow[1]
        #    & spotcalldf$Xm <=  boundingWindow[2]
        #    & spotcalldf$Ym >= boundingWindow[3]
        #    & spotcalldf$Ym <=  boundingWindow[4]) & underContention #From before
        ## NEW ##

        ## First check if competition has been done before
        ## If there are more spots in window from neighbour, filter out the current FOV
        if( competedBefore ){
          if( as.logical(checkCompetition$refWon) ){
            spotcalldf <- spotcalldf[!underContention,]
          }
          next
        }

        if( (sum(underContention) < nrow(neighbourdf)) & !competedBefore){
          spotcalldf <- spotcalldf[!underContention,]
          newWinLose <- data.frame(
            'reference' = fovName,
            'query' = fovNeighbours[i],
            'refWon' = FALSE,
            'queryBound1' = boundingWindow[1],
            'queryBound2' = boundingWindow[2],
            'queryBound3' = boundingWindow[3],
            'queryBound4' = boundingWindow[4]
          )
        }else{
          newWinLose <- data.frame(
            'reference' = fovName,
            'query' = fovNeighbours[i],
            'refWon' = TRUE,
            'queryBound1' = boundingWindow[1],
            'queryBound2' = boundingWindow[2],
            'queryBound3' = boundingWindow[3],
            'queryBound4' = boundingWindow[4]
          )
        }

        if(is.null(redundancyWinnersLosers)){
          redundancyWinnersLosers <- newWinLose
        }else{
          redundancyWinnersLosers <- rbind( redundancyWinnersLosers, newWinLose )
        }
        data.table::fwrite(redundancyWinnersLosers, file = winLoseFile, row.names = F, append = F, quote = F, showProgress =F)
      }
      ## OLDER ##
    }
    if(nrow(spotcalldf)==0){
      if(verbose){ message(paste0('No spots left in this FOV...'), appendLF = F) }
      next
    }

    ## Load cell segment information
    cellsegfx <- cellsegfs[grepl( paste0(fovName, '.csv'), cellsegfs )]
    if(length(cellsegfx)==0){
      if(verbose) { message('No cell segmentation file found...Skipping cell annotation...', appendLF = F) }
      spotcalldf$CELLNAME <- NA
    }
    if(length(cellsegfx)>1){
      warning( paste0('More than one cell segmentation file specified for this FOV! Picking first: ', basename(cellsegfx[1]), '...'))
      cellsegfx <- cellsegfx[1]
    }
    if(length(cellsegfx)==1){
      cellsegdf <- try(data.table::fread(cellsegfx, data.table = F, showProgress =F))
      if(inherits(cellsegdf, 'try-error')){
        stop('Unable to read cell segmentation file!')
      }
      coords <-
        cellsegdf$WX * as.numeric(resolutions$per_pixel_microns[1])  +
        1i * cellsegdf$WY * as.numeric(resolutions$per_pixel_microns[2]) +
        gcx[match(cellsegdf$fov, gcx$fov), 'x_microns'] +
        1i * gcx[match(cellsegdf$fov, gcx$fov), 'y_microns']
      cellsegdf$Xm <- Re(coords)
      cellsegdf$Ym <- Im(coords)
      cellsegdf$IDX <- cellsegdf$WX + (cellsegdf$WY-1) * as.numeric(resolutions$xydimensions_pixels[1])

      NCHAR_NAME <- pmax(nLeadingZeroes, nchar(max(cellsegdf$cell[!is.na(cellsegdf$cell)])) + 2)
      cellsegdf$CELLNAME <- paste0(cellsegdf$fov, '_', sprintf( paste0('%0', NCHAR_NAME, 'd'), cellsegdf$cell ))

      ## Load neighbours
      if(removeRedundancy){
        if(verbose){ message('\nRemoving redundant cell calls...', appendLF = F) }
        fovNeighbours <- stitchResults[grepl(paste0('STITCH_', fovName, '.csv'), names(stitchResults))][[1]]

        ## Also propagate to diagonal neighbours
        ## Because of earlier checks, can make things simpler here

        firstDegreeNeighbours <- do.call(rbind, lapply( fovNeighbours$fov, function(x){
          tmpdf <- stitchResults[grepl(paste0('STITCH_', x, '.csv'), names(stitchResults))][[1]]
          colnames(tmpdf) <- tolower(colnames(tmpdf))
          if( sum(colnames(tmpdf)=='fov') != 1){
            stop('STITCH files need to have a single column named "fov"!')
          }
          tmpdf <- tmpdf[tmpdf[,'fov']!=fovName,]
          tmpdf[,'shift_pixel'] = tmpdf[,'shift_pixel'] + fovNeighbours[fovNeighbours$fov==x,'shift_pixel']
          return(tmpdf)
        }) )
        fovNeighbours <- rbind(fovNeighbours, firstDegreeNeighbours)
        fovNeighbours <- fovNeighbours[!duplicated(fovNeighbours$fov),]

        cellDropCheck <- NULL
        if(file.exists(cellDropFile)){
          cellDropCheck <- data.table::fread(cellDropFile, data.table = F, showProgress =F)$CELLNAME
        }
        cellsegdf <- cellsegdf[!(cellsegdf$CELLNAME %in% cellDropCheck),]

        for( i in 1:nrow(fovNeighbours) ){

          ## Load neighbour cell segment calls
          neighbourfx <- cellsegfs[grepl( paste0(fovNeighbours$fov[i], '.csv'), cellsegfs )]
          if( length(neighbourfx)==0 ){ next }
          neighbourdf <- try(data.table::fread(neighbourfx, data.table = F, showProgress =F))
          if(inherits(neighbourdf, 'try-error')){
            stop('Unable to read neighbouring cell segmentation file!')
          }

          neighbourdf$WX <- neighbourdf$WX + Re( round(fovNeighbours[i, 'shift_pixel']) )
          neighbourdf$WY <- neighbourdf$WY + Im( round(fovNeighbours[i, 'shift_pixel']) )
          neighbourdf$IDX <- neighbourdf$WX + (neighbourdf$WY-1) * as.numeric(resolutions$xydimensions_pixels[1])
          NEIGHBOUR_CHAR_NAME <- pmax(nLeadingZeroes, nchar(max(neighbourdf$cell[!is.na(neighbourdf$cell)])) + 2)
          neighbourdf$CELLNAME <- paste0(neighbourdf$fov, '_', sprintf( paste0('%0', NEIGHBOUR_CHAR_NAME, 'd'), neighbourdf$cell ))
          if(
            !all(colnames(neighbourdf)==colnames(cellsegdf))
            ){
            stop('Column names between neighbour cell segmentation and current FOV dataframes not matching!')
          }
          # neighbourdf <- neighbourdf[,colnames(cellsegdf)]

          ## Drop cells that have lost before:
          neighbourdf <- neighbourdf[!(neighbourdf$CELLNAME %in% cellDropCheck),]

          ## These contain just the pixels that overlap with our FOV
          neighbourdfWindowed <- neighbourdf[
            neighbourdf$WX >= min(cellsegdf$WX)
            & neighbourdf$WX <= max(cellsegdf$WX)
            & neighbourdf$WY >= min(cellsegdf$WY)
            & neighbourdf$WY <= max(cellsegdf$WY),
          ]

          neighbourCells <- unique( neighbourdfWindowed$cell  )
          if( length(neighbourCells) == 0 ){ next }

          ## Plot to check
          # library(ggplot2)
          # ggplot() +
          #   geom_point(data=cellsegdf, aes(x=WX, y=WY), colour='red', pch='.', alpha=0.5) +
          #   geom_point(data=neighbourdf, aes(x=WX, y=WY), colour='blue', pch='.', alpha=0.5) +
          #   scale_y_reverse() +
          #   theme_minimal(base_size=20)

          ## This contains all pixels belonging to relevant cells
          neighbourdf <- neighbourdf[neighbourdf$cell %in% neighbourCells,]

          ## List cells we have to care about
          boundingWindow <- c(
            range(neighbourdfWindowed$WX[neighbourdfWindowed$cell %in% neighbourCells]),
            range(neighbourdfWindowed$WY[neighbourdfWindowed$cell %in% neighbourCells]) )
          underContentionCells <-
            cellsegdf[
              (cellsegdf$WX >= boundingWindow[1]
               & cellsegdf$WX <=  boundingWindow[2]
               & cellsegdf$WY >= boundingWindow[3]
               & cellsegdf$WY <=  boundingWindow[4]),
              'cell'
            ]
          underContentionCells <- sort(unique(underContentionCells))

          ## Plot to check
          # library(ggplot2)
          # ggplot() +
          #   geom_point(data=cellsegdf[cellsegdf$cell %in% underContentionCells,], aes(x=WX, y=WY), colour='red', pch=16, alpha=0.5) +
          #   # geom_point(data=neighbourdf, aes(x=WX, y=WY), colour='blue', pch=16, alpha=0.5) +
          #   scale_y_reverse() +
          #   theme_minimal(base_size=20)

          
          for( celli in 1:length(underContentionCells) ){
            ## Pick the bigger of two cells
            contentionCellIDX <- cellsegdf[cellsegdf$cell==underContentionCells[i],'IDX']
            competeingCell <- sort( unique(neighbourdfWindowed[neighbourdfWindowed[,'IDX'] %in% contentionCellIDX,'cell']) )
            competeingCellSizes <- table(neighbourdf[neighbourdf$cell %in% competeingCell, 'cell'])
            isThereConsiderableOverlap <- sapply(competeingCell, function(x){
              pixLocs <- neighbourdf[neighbourdf$cell==x,'IDX']
              AinB = length(intersect(pixLocs,contentionCellIDX)) / length(pixLocs) #Jaccard for competing cell in contention cell
              BinA = length(intersect(pixLocs,contentionCellIDX)) / length(contentionCellIDX) #Jaccard for vice versa
              considerableOverlap <- (AinB > cellOverlapFraction) | (BinA > cellOverlapFraction)
              return(considerableOverlap)
            })

            ## Plot to check
            # library(ggplot2)
            # ggplot() +
            #   geom_point(data=cellsegdf[cellsegdf$IDX %in% contentionCellIDX,], aes(x=WX, y=WY), colour='red', pch=16, alpha=0.5) +
            #   geom_point(data=neighbourdf[neighbourdf$cell %in% competeingCell,], aes(x=WX, y=WY), colour='blue', pch=16, alpha=0.5) +
            #   scale_y_reverse() +
            #   theme_minimal(base_size=20)

            ## When overlap is small, we ignore pixels
            if( any(!isThereConsiderableOverlap) ){
              ignorePixels <- neighbourdfWindowed[neighbourdfWindowed$cell %in% competeingCell[!isThereConsiderableOverlap], 'IDX']
              cellsegdf <- cellsegdf[!(cellsegdf[,'IDX'] %in% ignorePixels),]
            }

            ## When overlap is big, have to choose which cell to drop: drop the smaller nucleus, or if the nucleus has multiple overlaps
            if( any(isThereConsiderableOverlap) ){
              competeingCellSizes <- competeingCellSizes[names(competeingCellSizes) %in% competeingCell[isThereConsiderableOverlap]]
              correspondingCELLNAMES <- neighbourdfWindowed[match(names(competeingCellSizes), neighbourdfWindowed$cell),'CELLNAME']
              if( length(competeingCellSizes) > 1 | any(competeingCellSizes > length(contentionCellIDX)) ){
                filtout <- (cellsegdf$cell == underContentionCells[i])
                cellDropCheck <- unique( c(cellDropCheck, unique(cellsegdf[cellsegdf$cell == underContentionCells[i],'CELLNAME'])) )
                cellsegdf <- cellsegdf[!filtout,]
                cellsegdf <- rbind(
                  neighbourdfWindowed[neighbourdfWindowed$cell %in% competeingCell[isThereConsiderableOverlap],],
                  cellsegdf)
                cellsegdf <- cellsegdf[!duplicated(cellsegdf[,'IDX']),]
              }else{
                cellDropCheck <- unique(c(cellDropCheck, unique(correspondingCELLNAMES)))
              }
            }
          }
        }
        
        ## Remember which cells to drop in future
        data.table::fwrite(data.frame('CELLNAME' = unique(cellDropCheck)), cellDropFile,
                           row.names = F, quote = F, append = F, showProgress = F)
      }

      ## Transfer cellnames to spotcalldf
      ## IDX for cross referencing with spotcalldf, CIDX for cross referencing with other cellsegdf
      spotcalldf$CELLNAME <- cellsegdf[match(spotcalldf$IDX, cellsegdf$IDX),'CELLNAME']
    }

    ## Append information with fwrite
    data.table::fwrite(
      spotcalldf,
      paste0(params$out_dir, 'OUT_SPOTCALL_PIXELS.csv', ifelse(gZip, '.gz', '')),
      row.names = F, append = T, showProgress = F)
    data.table::fwrite(
      cellsegdf,
      paste0(params$out_dir, 'OUT_CELLSEG_PIXELS.csv', ifelse(gZip, '.gz', '')),
      row.names = F, append = T, showProgress = F)

    ## Summarise info per cell
    if(!all(is.na(spotcalldf$CELLNAME)) & nrow(cellsegdf) > 0){
      cellNames <- sort(unique(cellsegdf$CELLNAME))
      cellMeta <- data.frame(
        'CELLNAME' = cellNames,
        'Xm' = as.numeric(by(cellsegdf$Xm, factor(cellsegdf$CELLNAME, levels=cellNames), median )),
        'Ym' = as.numeric(by(cellsegdf$Ym, factor(cellsegdf$CELLNAME, levels=cellNames), median )),
        'nPixels' = as.numeric(table(factor(cellsegdf$CELLNAME, levels=cellNames)))
      )
      spotcalldf <- spotcalldf[!is.na(spotcalldf$CELLNAME) & (spotcalldf$CELLNAME %in% cellMeta$CELLNAME),]
      cellExp <- do.call(rbind, by( factor(spotcalldf$g, levels=g), factor(spotcalldf$CELLNAME, cellNames), table ) )
      cellExp <- data.frame( cellExp[match(cellMeta$CELLNAME, rownames(cellExp)),], check.names = F)
      rownames(cellExp) <- cellMeta$CELLNAME
      colnames(cellExp) <- g
      cellExp[is.na(cellExp)] <- 0
      nCounts <- rowSums(cellExp)
      cellMeta$nCounts <- as.numeric(nCounts)[match(cellMeta$CELLNAME, names(nCounts))]
      cellMeta$nCounts[is.na(cellMeta$nCounts)] <- 0
      data.table::fwrite(
        cellExp,
        paste0(params$out_dir, 'OUT_CELLEXPRESSION.csv', ifelse(gZip, '.gz', '')),
        row.names = T, append = T, showProgress = F)
      data.table::fwrite(
        cellMeta,
        paste0(params$out_dir, 'OUT_CELLS.csv', ifelse(gZip, '.gz', '')),
        row.names = F, append = T, showProgress = F)
    }

    if(verbose){ message('Done!', appendLF = F) }
  }

  ## Note that above tends to create duplicate entries in CELLEXPRESSION, have to load and edit
  if(verbose){ message('\nCleaning up...') }

  if(file.exists(winLoseFile)){ file.remove( winLoseFile ) }
  # if(file.exists(cellDropFile)){ file.remove( cellDropFile ) }

  cellMetaClean <- cellMeta <- try( data.table::fread(
    paste0(params$out_dir, 'OUT_CELLS.csv', ifelse(gZip, '.gz', '')),
    data.table = F, showProgress = F) )
  if( inherits(cellMeta, 'try-error') ){
    warning('OUT_CELLS file was not created / could not be read! Unable to perform clean up!')
  }else{
    if( any(duplicated(cellMeta$CELLNAME)) ){
      duplicatedCells <- sort(unique( cellMeta[duplicated(cellMeta$CELLNAME), 'CELLNAME'] ))
      cellMetaClean <- cellMeta[!(cellMeta$CELLNAME %in% duplicatedCells), ]
      dupMeta <- cellMeta[(cellMeta$CELLNAME %in% duplicatedCells), ]
      for( dupCelli in duplicatedCells ){
        npix <- sum( as.numeric(dupMeta[dupMeta$CELLNAME==dupCelli, 'nPixels']) )
        weightedX <- sum( dupMeta[dupMeta$CELLNAME==dupCelli, 'nPixels'] * dupMeta[dupMeta$CELLNAME==dupCelli, 'Xm'] / npix )
        weightedY <- sum( dupMeta[dupMeta$CELLNAME==dupCelli, 'nPixels'] * dupMeta[dupMeta$CELLNAME==dupCelli, 'Ym'] / npix )
        res <- data.frame(
          'CELLNAME' = dupCelli,
          'Xm' =  weightedX,
          'Ym' = weightedY,
          'nPixels' = npix,
          'nCounts' = sum(as.numeric(dupMeta[dupMeta$CELLNAME==dupCelli, 'nCounts']))
        )
        cellMetaClean <- rbind(cellMetaClean, res)
      }
    }
    cellMeta <- cellMetaClean
    data.table::fwrite(
      cellMeta,
      paste0(params$out_dir, 'OUT_CELLS.csv', ifelse(gZip, '.gz', '')),
      row.names = T, append = F, showProgress = F)
    suppressWarnings(rm(list=c('cellMetaClean', 'npix', 'weightedX', 'weightedY', 'res', 'duplicatedCells')))
  }

  cellExp <- try( data.table::fread(
    paste0(params$out_dir, 'OUT_CELLEXPRESSION.csv', ifelse(gZip, '.gz', '')),
    data.table = F, showProgress = F) )
  if( inherits(cellExp, 'try-error') ){
    warning('OUT_CELLEXPRESSION file was not created / could not be read! Unable to perform clean up!')
  }else{
    duplicatedCells <- sort(unique(cellExp[duplicated(cellExp[,1]),1]))
    dupCellExp <- cellExp[(cellExp[,1] %in% duplicatedCells),]
    cellExpClean <- cellExp[!(cellExp[,1] %in% duplicatedCells),]
    rownames(cellExpClean) <- cellExpClean[,1]
    cellExpClean[,1] <- NULL
    if( nrow(dupCellExp) > 0){
      newVals <- data.frame(matrix(0, nrow=length(duplicatedCells), ncol=ncol(cellExp)-1), row.names = duplicatedCells)
      colnames(newVals) <- colnames(cellExp)[-1]
      for( dupCelli in duplicatedCells ){
        newVals[dupCelli,] <- colSums( dupCellExp[dupCellExp[,1]==dupCelli,-1] )
      }
      cellExpClean <- rbind(cellExpClean, newVals)
    }
    cellExpClean <- cellExpClean[match(cellMeta$CELLNAME, rownames(cellExpClean)),]
    cellExpClean[is.na(cellExpClean)] <- 0
    rownames(cellExpClean) <- cellMeta$CELLNAME
    data.table::fwrite(
      cellExpClean,
      paste0(params$out_dir, 'OUT_CELLEXPRESSION.csv', ifelse(gZip, '.gz', '')),
      row.names = T, append = F, showProgress = F)
  }
}

