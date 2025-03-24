## This script contains functions for preparing the codebook.

require(viridis)
require(data.table)

prepareCodebook <- function(
    params = get('params', envir = globalenv()),
    exhaustiveBlanks = T,
    hammingDistanceThreshold = 2
){

  if(!file.exists(paste0(params$out_dir, 'ORDEREDCODEBOOK.csv.gz')) | !params$resumeMode){

    ## Load existing codebook
    if(!exists('raw_codebook', envir = params)){
      stop('Codebook not found! Ensure that establishParams and readImageMetaData has been run!')
    }
    codebook <- params$raw_codebook

    ## Capture order
    if(!exists('global_coords', envir = params)){
      stop('params$global_coords not found! Ensure readImageMetaData has been run!')
    }
    global_coords <- params$global_coords
    check_colnames <- params$codebook_colnames[params$codebook_colnames!='gene']
    if(!all(check_colnames %in% global_coords$bit_name)){
      global_coords$bit_name <- global_coords$bit_name_full
      if(!all(check_colnames %in% global_coords$bit_name)){
        stop('bit_name not matching codebookColumnNames: either edit GLOBALCOORD.csv or params$codebook_colnames. Then re-run readImageMetaData.')
      }
    }
    gcx <- params$global_coords
    nimages <- table(gcx$fov)
    chosenFOV <- names(nimages)[which.max(nimages)][1]
    gcx <- gcx[gcx$fov==chosenFOV,]
    gcx <- gcx[!duplicated(gcx$bit_name),]
    gcx <- gcx[gcx$z_microns == gcx$z_microns[1],]
    capture_order <- gcx$bit_name_alias
    capture_order <- order(strsplit2(capture_order, '_')[,1], factor(strsplit2(capture_order, '_')[,2], levels=params$resolutions$channel_order) )
    capture_order <- gcx[capture_order, 'bit_name']
    params$capture_order <<- capture_order

    ## Ensure that capture_order matches colnames of codebook
    codebook <- codebook[,match(capture_order, colnames(codebook))]
    if(!all(colnames(codebook) == capture_order )){
      stop('Unable to order codebook!')
    }

    if(exhaustiveBlanks){
      message('\nFinding more blanks...')
      universe <- expand.grid(lapply(1:ncol(codebook), function(x) return(c(0,1)) ))
      intcodes <- Reduce('+', lapply(1:ncol(codebook), function(x){
        val <- universe[,x] * 2^(x-1)
        return(val)
      }))
      donecodes <-  Reduce('+', lapply(1:ncol(codebook), function(x){
        val <- codebook[,x] * 2^(x-1)
        return(val)
      }))
      newblanks <- universe[!(intcodes %in% donecodes),]
      newblanks <- newblanks[rowSums(newblanks)==params$nbits,]

      message('\nCalculating hamming distances to existing codes...')
      hdists <- rep(Inf, nrow(newblanks))
      closest <- rep(NA, nrow(newblanks))
      for(i in 1:nrow(codebook)){
        message(paste0(i, ' of ', nrow(codebook), '...'), appendLF = F)
        ref <- do.call(rbind, lapply(1:nrow(newblanks), function(x) return( codebook[i,] )))
        hdistx <- rowSums(abs(newblanks - ref))
        bool <- (hdists - hdistx) > 0
        hdists[ bool ] <- hdistx[ bool ]
        closest[bool] <- rownames(codebook)[i]
      }
      valid <- (hdists > hammingDistanceThreshold)
      newblanks <- newblanks[valid,]

      message(paste0('\nAdding ', sum(valid), ' new blanks...'))
      if(sum(valid) > 0){
        rownames(newblanks) <- paste0('BLANK-NEW-', sprintf(paste0('%0', nchar(nrow(newblanks)), 'd'), 1:nrow(newblanks)))
        colnames(newblanks) <- colnames(codebook)
        codebook <- rbind(codebook, newblanks)
      }
    }

    write.csv(codebook, paste0(params$out_dir, 'ORDEREDCODEBOOK.csv.gz'), row.names = T)
  }

  codebook <- data.table::fread(paste0(params$out_dir, 'ORDEREDCODEBOOK.csv.gz'), data.table=FALSE)
  if(colnames(codebook)[1] == 'V1'){
    rownames(codebook) <- codebook$V1
    codebook$V1 <- NULL
  }
  params$ordered_codebook <<- codebook
  params$capture_order <<- colnames(codebook)
  params$isblank <<- grepl('^blank-', tolower(rownames(codebook)))
  genePalette <- setNames( viridis::turbo(nrow(codebook), alpha = 1), rownames(codebook) )
  params$genePalette <<- genePalette
  message(paste0('\nCodebook has ', sum(params$isblank), ' blanks and ', sum(!params$isblank), ' genes...' ))
}
