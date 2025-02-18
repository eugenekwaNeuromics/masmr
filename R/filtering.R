# Function for dataframe filtering

require(data.table)

filterDF <- function(
    spotcalldf,
    filterOut,
    logBase = exp(1),
    params = get('params', envir = globalenv())
){

  if(length(filterOut) != nrow(spotcalldf)){
    stop('filterOut length does not match nrows of the dataframe!')
  }
  filterOut <- as.logical(filterOut)
  if(any(is.na(filterOut))){
    stop('filterOut is not a boolean vector!')
  }
  keep <- !filterOut

  message('Filtering...')

  ## Report number of pixels filtered
  before <- length(keep)
  after <- sum(keep)
  before_denom <- before
  after_denom <- before
  before_perc <- round(100 * before / before_denom, digits = 2)
  after_perc <- round(100 * after / after_denom, digits = 2)
  message(paste0(
    'N pixels: ',
    before, ' ( ', before_perc, '% ) --> ',
    after, ' ( ', after_perc, '% )'))

  ## Report number of blanks filtered
  if('BLANK' %in% colnames(spotcalldf)){
    isblank <- as.logical(spotcalldf$BLANK)
    if(!any(is.na(isblank))){
      before <- sum(isblank)
      after <- sum(isblank[keep])
      before_denom <- length(keep)
      after_denom <- sum(keep)
      before_perc <- round(100 * before / before_denom, digits = 2)
      after_perc <- round(100 * after / after_denom, digits = 2)
      message(paste0(
        'N blanks: ',
        before, ' ( ', before_perc, '% ) --> ',
        after, ' ( ', after_perc, '% )'))
    }
  }

  ## Report FPKM correlation
  if( exists('fpkm_file', envir = params)){
    if( !exists('fpkm', envir = params) | object.size(params$fpkm)==0 ){
      fpkm <- try( data.table::fread(params$fpkm_file, data.table = F) )
      if(!inherits(fpkm, 'try-error')){
        if(ncol(fpkm)==0 | ncol(fpkm)>2){
          warning('Unable to parse FPKM file')
          fpkm <- NULL
        }
        if( any(ncol(fpkm)==1) ){
          if(all( tolower(rownames(fpkm)) %in% tolower(rownames(params$ordered_codebook)))){
            fpkm <- setNames(fpkm[,1], rownames(fpkm))
          }
        }
        if( any(ncol(fpkm)==2) ){
          gCol <- sapply(1:ncol(fpkm), function(i){
            all(tolower(fpkm[,i]) %in% tolower(rownames(params$ordered_codebook)))
          })
          fpkm <- setNames(as.numeric(fpkm[,!gCol]), as.character(fpkm[,gCol]))
        }
        if(!is.null(fpkm)){
          fpkm <- fpkm[ match(tolower(rownames(params$ordered_codebook)), tolower(names(fpkm))) ]
          fpkm[is.na(fpkm)] <- 0
          names(fpkm) <- rownames(params$ordered_codebook)
          params$fpkm <<- fpkm
        }
      }
    }
    fpkm <- params$fpkm
    if(!is.null(fpkm)){
      lfpkm <- log(1 + fpkm, base = logBase)
      ng <- table(factor(spotcalldf$g, levels=rownames(params$ordered_codebook)))
      ng[is.na(ng)] <- 0
      ng <- log(1+ng, base = logBase)
      before <- cor(lfpkm, ng)

      ng <- table(factor(spotcalldf$g[keep], levels=rownames(params$ordered_codebook)))
      ng[is.na(ng)] <- 0
      ng <- log(1+ng, base = logBase)
      after <- cor(lfpkm, ng)

      message(paste0(
        'LogFPKM Corr: ',
        round(before, digits=5), ' --> ',
        round(after, digits=5) ))
    }
  }

  ## N probes correlation
  if( exists('n_probes_file', envir = params)){
    if( !exists('n_probes', envir = params) | object.size(params$n_probes)==0 ){
      n_probes <- try( data.table::fread(params$n_probes_file, data.table = F) )
      if(!inherits(n_probes, 'try-error')){
        if(ncol(n_probes)==0 | ncol(n_probes)>2){
          warning('Unable to parse FPKM file')
          n_probes <- NULL
        }
        if( any(ncol(n_probes)==1) ){
          if(all( tolower(rownames(n_probes)) %in% tolower(rownames(params$ordered_codebook)))){
            n_probes <- setNames(n_probes[,1], rownames(n_probes))
          }
        }
        if( any(ncol(n_probes)==2) ){
          gCol <- sapply(1:ncol(n_probes), function(i){
            all(tolower(n_probes[,i]) %in% tolower(rownames(params$ordered_codebook)))
          })
          n_probes <- setNames(as.numeric(n_probes[,!gCol]), as.character(n_probes[,gCol]))
        }
        if(!is.null(n_probes)){
          n_probes <- n_probes[ match(tolower(rownames(params$ordered_codebook)), tolower(names(n_probes))) ]
          n_probes[is.na(n_probes)] <- 0
          names(n_probes) <- rownames(params$ordered_codebook)
          params$n_probes <<- n_probes
        }
      }
    }
    n_probes <- params$n_probes
    if(!is.null(n_probes)){
      ln <- log(1 + n_probes, base = logBase)
      ng <- table(factor(spotcalldf$g, levels=rownames(params$ordered_codebook)))
      ng[is.na(ng)] <- 0
      ng <- log(1+ng, base = logBase)
      before <- cor(ln, ng)

      ng <- table(factor(spotcalldf$g[keep], levels=rownames(params$ordered_codebook)))
      ng[is.na(ng)] <- 0
      ng <- log(1+ng, base = logBase)
      after <- cor(ln, ng)

      message(paste0(
        'LogN Probes Corr: ',
        round(before, digits=5), ' --> ',
        round(after, digits=5) ))
    }
  }

  return(spotcalldf[keep,])
}
