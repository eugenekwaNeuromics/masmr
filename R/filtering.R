# Function for dataframe filtering

require(rlang)
require(data.table)
require(ggplot2)
require(viridis)

spotcall_troubleshootPlots <- function(
    spotcalldf,
    imList,
    chosenCoordinate,
    plotWindowRadius = 10,
    decodeMetric = 'DECODE',
    params = get('params', envir = globalenv())
){
  if(!('g' %in% colnames(spotcalldf)) ){
    stop('spotcalldf must have a g column: ensure decoding has been done!')
  }
  if(!missing(chosenCoordinate)){
    if(!is.complex(chosenCoordinate)){
      stop('Please specify chosenCoordinate as a complex number: with real being X and imaginary being Y!')
    }
  }
  if(missing(imList)){
    ## Get decode
    if( !exists('imMetrics', envir=globalenv()) ){
      stop('Unable to plot: imMetrics does not exist and imageList not provided!')
    }
    imMetrics <- get('imMetrics', envir = globalenv())
    if( !exists(decodeMetric, envir=imMetrics) ){
      if(length(names(imMetrics))<1 | is.null(names(imMetrics))){
        stop('Unable to plot: cannot find relevant object in imMetrics to reference!')
      }
      warning( paste0('Unable to find ', decodeMetric, ': replacing with ', names(imMetrics)[1], '...'))
      decodeMetric <- names(imMetrics)[1]
    }
    imList <- get(decodeMetric, envir=imMetrics)
  }
  shifts <- params$shifts
  if(is.null(shifts)){
    warning('Unable to find params$shifts: will not adjust coordinates')
    shifts <- rep(0+0i, length(imList))
  }
  if( !all( c('WX', 'WY') %in% colnames(spotcalldf)) ){
    stop('Unable to plot: cannot find columns WX and WY in your spotcalldf!')
  }

  if(is.null( params$genePalette )){
    genePalette <- setNames( viridis::turbo(nrow(params$ordered_codebook), alpha = 1), rownames(params$ordered_codebook) )
    params$genePalette <<- genePalette
  }
  genePalette <- params$genePalette

  if(!missing(chosenCoordinate)){
    chosen_cx <- chosenCoordinate
  }else{
    pixelIndex <- 1
    chosen_cx <- spotcalldf[pixelIndex,'WX'] + 1i * spotcalldf[pixelIndex,'WY']
  }
  coordM <- getRasterCoords(imList[[ which(shifts==(0+0i))[1] ]]) + 1 + 1i
  MAXX <- max(Re(coordM))
  MAXY <- max(Im(coordM))

  plotList <- list()
  for(j in 1:length(chosen_cx)){
    cxx <- chosen_cx[j]
    bbox <- c(
      Re(cxx) - plotWindowRadius,
      Re(cxx) + plotWindowRadius,
      Im(cxx) - plotWindowRadius,
      Im(cxx) + plotWindowRadius
    )
    subimList <- do.call(rbind, lapply(1:length(imList), function(i){
      shifti <- shifts[i]
      subIm <- imList[[i]]
      dfsub <- as.data.frame(imager::as.cimg(subIm))
      dfsub$x <- dfsub$x + Re(shifti)
      dfsub$y <- dfsub$y + Im(shifti)
      dfsub <- dfsub[
        dfsub$x >= bbox[1]
        & dfsub$x <= bbox[2]
        & dfsub$y >= bbox[3]
        & dfsub$y <= bbox[4],
      ]
      bname <- names(imList)[i]
      if(is.null(bname)){
        bname = i
        lvls <- 1:length(imList)
      }else{
        lvls <- names(imList)
      }
      dfsub$bit_name <- factor(bname, levels = lvls)
      return(dfsub)
    }) )
    spdf_bool <-
      spotcalldf$WX >= bbox[1] &
      spotcalldf$WX <= bbox[2] &
      spotcalldf$WY >= bbox[3] &
      spotcalldf$WY <= bbox[4]
    glevels = sort(unique(as.character(spotcalldf[spdf_bool, 'g'])))
    gene_pal <- genePalette[names(genePalette) %in% glevels]

    p <-
      ggplot2::ggplot() +
      ggplot2::geom_tile(
        data=subimList,
        ggplot2::aes(x=x, y=y, fill=value)) +
      ggplot2::geom_point(
        data=spotcalldf[spdf_bool,],
        ggplot2::aes(x=WX, y=WY, colour=factor(g, levels = glevels) ),
        shape=1, alpha = 1, stroke=1) +
      ggplot2::scale_colour_manual(name = '', values=gene_pal, na.value='none') +
      ggplot2::scale_fill_gradient(low='black', high='white', na.value = 'black') +
      ggplot2::facet_wrap( ~factor(bit_name) ) +
      ggplot2::theme_void(base_size=14) +
      ggplot2::coord_fixed() +
      ggplot2::xlab('') + ggplot2::ylab('') +
      ggplot2::theme(
        plot.background = ggplot2::element_rect( fill = 'black' ),
        text = ggplot2::element_text(colour = 'white'),
        strip.text = ggplot2::element_text(colour = 'white', size=16),
        axis.text = ggplot2::element_text(colour = 'grey', size=12)) +
      ggplot2::guides(
        fill='none',
        colour = ggplot2::guide_legend(
          override.aes = list(
            size = 10,
            alpha = 1, shape = 16
          )))
    spot_name <- paste0('X', Re(cxx), '_Y', Im(cxx))
    plotList[[spot_name]] <- p
  }

  return(plotList)
}

###

filterDF <- function(
    spotcalldf,
    filterOut,
    returnTroubleShootPlots = FALSE,
    troubleShootCoordinates = NULL,
    logBase = exp(1),
    params = get('params', envir = globalenv()),
    ...
){
  filterExpression <- rlang::as_label(rlang::enquo( filterOut ))
  message( paste0('\nApplying filter: ', filterExpression, '...') )

  if(length(filterOut) != nrow(spotcalldf)){
    stop('filterOut length does not match nrows of the dataframe!')
  }
  filterOut <- as.logical(filterOut)
  if(any(is.na(filterOut))){
    stop('filterOut is not a boolean vector!')
  }
  keep <- !filterOut

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

  if( returnTroubleShootPlots ){
    troubleshootPlots <<- new.env()
    if(is.null(troubleShootCoordinates)){
      refcx <- spotcalldf[keep,]
      refcx <- refcx[which.min(refcx$COS),c('WX','WY')]
      refcx <- refcx$WX + refcx$WY * 1i
      plist <- spotcall_troubleshootPlots(
        spotcalldf=spotcalldf,
        params=params,
        chosenCoordinate=refcx,
        ...)
      troubleshootPlots[['BEFORE']] <<- plist
      plist <-  spotcall_troubleshootPlots(
        spotcalldf=spotcalldf[keep,],
        params=params,
        chosenCoordinate=refcx,
        ...)
      troubleshootPlots[['AFTER']] <<- plist
    }else{
      plist <- spotcall_troubleshootPlots(
        spotcalldf=spotcalldf,
        params=params,
        chosenCoordinate = troubleShootCoordinates,
        ...)
      troubleshootPlots[['BEFORE']] <<- plist
      plist <-  spotcall_troubleshootPlots(
        spotcalldf=spotcalldf[keep,],
        params=params,
        chosenCoordinate = troubleShootCoordinates,
        ...)
      troubleshootPlots[['AFTER']] <<- plist
    }
  }

  return(spotcalldf[keep,])
}
