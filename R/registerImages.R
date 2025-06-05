# This script contains the function for image registration.

require(ggplot2)
require(imager)

register_troubleshootPlots <- function(
    imList,
    chosenCoordinate,
    plotWindowRadius = 50,
    params = get('params', envir = globalenv())
){
  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose
  }else{
    verbose = T
  }
  if(!missing(chosenCoordinate)){
    if(!is.complex(chosenCoordinate)){
      stop('Please specify chosenCoordinate as a complex number: with real being X and imaginary being Y!')
    }
  }
  shifts <- params$shifts
  if(is.null(shifts)){
    warning('Unable to find params$shifts: will not adjust coordinates')
    shifts <- rep(0+0i, length(imList))
  }

  if(!missing(chosenCoordinate)){
    chosen_cx <- chosenCoordinate
  }else{
    imsub <- imList[[ which(shifts==0+0i)[1] ]]
    imsub <- maxIntensityProject(imsub)
    chosen_cx <- getRasterCoords(imsub)[which.max(imsub)]
  }
  
  for(chosen_cxi in chosen_cx){
    dfp <- do.call(rbind, lapply(1:length(imList), function(ix){
      imi <- imList[[ix]]
      imName <- names(imList)[ix]
      if(length(dim(imi))>2){
        imi <- maxIntensityProject(imi)
      }
      if(is.null(imName)){
        imName <- sprintf( paste0('%0', nchar(length(imList)), 'd'), ix)
      }
      dfx <- as.data.frame(imager::as.cimg(imi))
      dfx[,c('orig_x', 'orig_y')] <- dfx[,c('x', 'y')]
      dfx$x <- dfx$x + Re(shifts[ix])
      dfx$y <- dfx$y + Im(shifts[ix])
      
      dfx <- dfx[dfx$x > (Re(chosen_cxi) - plotWindowRadius),]
      dfx <- dfx[dfx$x <= (Re(chosen_cxi) + plotWindowRadius),]
      dfx <- dfx[dfx$y > (Im(chosen_cxi) - plotWindowRadius),]
      dfx <- dfx[dfx$y <= (Im(chosen_cxi) + plotWindowRadius),]
      dfx$imName <- imName
      
      dfx$norm <- dfx$value - min(dfx$value)
      dfx$norm <- dfx$norm / max(dfx$norm)
      
      return(dfx)
    }) )
    refdfp <- names( which(shifts==0+0i)[1] )
    if(is.null(refdfp)){
      refdfp <- sprintf( paste0('%0', nchar(length(imList)), 'd'), which(shifts==0+0i)[1] )
    }
    refdfp <- dfp[dfp$imName==refdfp,]
    refdfp$imName <- NULL
    
    p_before <-
      ggplot2::ggplot() +
      ggplot2::geom_raster(data=refdfp, ggplot2::aes(x=x, y=y, alpha=norm), fill='blue') +
      ggplot2::geom_raster(data=dfp, ggplot2::aes(x=orig_x, y=orig_y, alpha=norm), fill='red') +
      ggplot2::facet_wrap(~imName) +
      ggplot2::scale_alpha(range=c(0, 0.75)) +
      ggplot2::theme_minimal(base_size=12) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::scale_y_reverse() +
      ggplot2::xlab('') + ggplot2::ylab('') +
      ggplot2::coord_fixed()
    
    p_after <-
      ggplot2::ggplot() +
      ggplot2::geom_raster(data=refdfp, ggplot2::aes(x=x, y=y, alpha=norm), fill='blue') +
      ggplot2::geom_raster(data=dfp, ggplot2::aes(x=x, y=y, alpha=norm), fill='red') +
      ggplot2::facet_wrap(~imName) +
      ggplot2::scale_alpha(range=c(0, 0.75)) +
      ggplot2::theme_minimal(base_size=12) +
      ggplot2::theme(legend.position = 'none') +
      ggplot2::scale_y_reverse() +
      ggplot2::xlab('') + ggplot2::ylab('') +
      ggplot2::coord_fixed()
    
    plotList[[as.character(chosen_cxi)]] <- list('PRE_REGISTER' = p_before, 'POST_REGISTER' = p_after)
  }
  
  return(plotList)
}

registerImages <- function(
    imList,
    registerTo = 1,
    chosenZslice = NULL,
    params = get('params', envir = globalenv()),
    currentFovName = NULL,
    maxAcceptableShiftDistance = Inf,
    gaussianWeightedAlignment = T,
    gaussianWeightAmplitude = NULL,
    gaussianWeightSigmaX = NULL,
    gaussianWeightSigmaY = NULL,
    tolerableBestAlignmentPercentage = 50,
    standardiseShifts = F,
    saveShifts = T,
    returnTroubleShootPlots = F,
    ...
){

  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose
  }else{
    verbose = T
  }

  ## Check that chosenZslice is correctly specified
  if(!is.null(chosenZslice)){
    if(length(chosenZslice) > 1){
      warning('Can only specify 1 chosenZslice...Defaulting to first element...')
    }
    chosenZslice <- as.integer(chosenZslice)[1]
  }

  # Check that the current FOV is known
  if(is.null(currentFovName)){
    currentFovName = params$current_fov #Try to see if there exists current FOV in params
    if(is.null(currentFovName) & saveShifts){
      warning('currentFovName NOT specified: will not save shifts as their own files!')
      saveShifts = F
    }
  }

  # Specify image to align against
  if(registerTo <= 0 | registerTo > length(imList)){
    stop('Invalid registerTo value: needs to refer to an image in imList!')
  }
  ref_im <- imList[[registerTo]]
  if( any(dim(ref_im)==1) ){
    ref_im <- array(as.vector(ref_im), dim=dim(ref_im)[dim(ref_im)!=1] )
  }
  if( length(dim(ref_im)) > 2 ){

    ## Image stack detected
    if(is.null(chosenZslice)){
      warning('Image stack detected but chosenZslice not specified...Performing maxIntensityProject...')
      ref_im <- maxIntensityProject(ref_im)
      ref_im <- array(as.vector(ref_im), dim=dim(quer_im)[dim(ref_im)!=1] )
    }else{
      warning('Image stack detected...Selecting im[,,chosenZslice]...')
      ref_im <- ref_im[,,chosenZslice]
      ref_im <- array(as.vector(ref_im), dim=dim(ref_im)[dim(ref_im)!=1] )
    }
  }
  coord <- getRasterCoords(ref_im) - sum(dim(ref_im) * c(1, 1i)) + (1+1i)

  # Prepare the Gaussian weight
  if(is.null(gaussianWeightAmplitude)){
    gaussianWeightAmplitude = 1
  }
  if(is.null(gaussianWeightSigmaX)){
    gaussianWeightSigmaX = 0.2 * as.numeric(params$resolutions$xydimensions_pixels)[1]
  }
  if(is.null(gaussianWeightSigmaY)){
    gaussianWeightSigmaY = 0.2 * as.numeric(params$resolutions$xydimensions_pixels)[2]
  }
  w <- gaussianWeightAmplitude * exp(-( (Re(coord)^2/(2*gaussianWeightSigmaX^2) + (Im(coord)^2/(2*gaussianWeightSigmaY^2) ))))
  acceptable_coords <- coord[Mod(coord) <= maxAcceptableShiftDistance]

  if(verbose){ message('\nRegistering...') }
  shifts <- rep(0+0i, length(imList))
  for(i in 1:length(imList)){
    if(verbose){ message(paste0(i, ' of ', length(imList), '...'), appendLF = F) }
    if(i==registerTo){ next }
    quer_im <- imList[[i]]

    if( any(dim(quer_im)==1) ){
      quer_im <- array(as.vector(quer_im), dim=dim(quer_im)[dim(quer_im)!=1] )
    }

    if( length(dim(quer_im)) > 2 ){

      ## Image stack detected
      if(is.null(chosenZslice)){
        quer_im <- maxIntensityProject(quer_im)
        quer_im <- array(as.vector(quer_im), dim=dim(quer_im)[dim(quer_im)!=1] )
      }else{
        quer_im <- quer_im[,,chosenZslice]
        quer_im <- array(as.vector(quer_im), dim=dim(quer_im)[dim(quer_im)!=1] )
      }

    }

    ## FAST version: correlate images without padding
    corra <- crossCorrelate2D(ref_im, quer_im, normalized=FALSE, pad = F)
    corrb <- crossCorrelate2D(quer_im, ref_im, normalized=FALSE, pad = F)
    if(gaussianWeightedAlignment){
      corra <- corra * w
      corrb <- corrb * w
    }
    shifta <- coord[which.max(corra)]
    shiftb <- coord[which.max(corrb)]
    shift <- #This check is not really necessary, but better to be safe
      ifelse( abs(Re(shifta)) < (median(abs(Re(coord)))), Re(shifta) + 0i, 0+0i ) +
      ifelse( abs(Im(shifta)) < (median(abs(Im(coord)))), 0 + 1i * Im(shifta), 0+0i ) +
      ifelse( abs(Re(shiftb)) < (median(abs(Re(coord)))), -Re(shiftb) + 0i, 0+0i ) +
      ifelse( abs(Im(shiftb)) < (median(abs(Im(coord)))), 0 - 1i * Im(shiftb), 0+0i )

    ## Case for when multiple values have high correlations (e.g. largely empty image)
    if(length(shift)>1){
      warning('Multiple coordinates have max cross-correlation value! Picking coordinate with smallest shift...')
      shift <- shift[order(Mod(shift), decreasing = FALSE)][1]
    }

    ## Case for when the distance shifted is more than what the user will accept
    if(Mod(shift) > maxAcceptableShiftDistance){
      warning('Best match outside of maximum pixel shift search area...Will restrict search...')
      if(max(corra)==0 | max(corrb)==0){
        warning('Maximum cross-correlation is 0...Will not shift...')
        shift <- 0 + 0i
      }else{
        # Search for highest cross correlation value WITHIN maxAcceptableShiftDistance
        # Only accept if this second highest value is at least x% of highest correlation value

        best_matches <- setNames(corra[coord %in% acceptable_coords] / max(corra), coord[coord %in% acceptable_coords])
        best_matches <- best_matches[order(best_matches, decreasing = TRUE)]
        best_match_perc = round(100*(best_matches)[1], digits = 2)
        if(best_match_perc < tolerableBestAlignmentPercentage){
          if(verbose){ message('Match within search area <', tolerableBestAlignmentPercentage, '% of maximum...Will not shift...\n') }
          shifta <- 0 + 0i
        }else{
          if(verbose){ message( paste0('Shift A: Picking ', names(best_matches[1]), ' with a value ', best_match_perc, '% of maximum...\n')) }
          shifta <- as.complex(names(best_matches[1]))
        }

        best_matches <- setNames(corrb[coord %in% acceptable_coords] / max(corrb), coord[coord %in% acceptable_coords])
        best_matches <- best_matches[order(best_matches, decreasing = TRUE)]
        best_match_perc = round(100*(best_matches)[1], digits = 2)
        if(best_match_perc < tolerableBestAlignmentPercentage){
          if(verbose){ message('Match within search area <', tolerableBestAlignmentPercentage, '% of maximum...Will not shift...\n') }
          shiftb <- 0 + 0i
        }else{
          if(verbose){ message( paste0('Shift B: Picking ', names(best_matches[1]), ' with a value ', best_match_perc, '% of maximum...\n')) }
          shiftb <- as.complex(names(best_matches[1]))
        }

        shift <- #This check is not really necessary, but better to be safe
          ifelse( abs(Re(shifta)) < (median(abs(Re(coord)))), Re(shifta) + 0i, 0+0i ) +
          ifelse( abs(Im(shifta)) < (median(abs(Im(coord)))), 0 + 1i * Im(shifta), 0+0i ) +
          ifelse( abs(Re(shiftb)) < (median(abs(Re(coord)))), -Re(shiftb) + 0i, 0+0i ) +
          ifelse( abs(Im(shiftb)) < (median(abs(Im(coord)))), 0 - 1i * Im(shiftb), 0+0i )
      }
    }

    shifts[i] <- shift
  }
  if(verbose){ message('') }
  names(shifts) <- names(imList)

  # The end result is a vector of shifts, expressed as complex coordinates
  # Now, user decides if want additional normalisations
  new_shifts <- shifts
  if( standardiseShifts & any(is.null(names(shifts))) ){
    warning('Invalid image names found: will not perform standardiseShifts')
    standardiseShifts = F
  }else{
    channels <- unique(strsplit2(names(shifts), '_')[,1])
    ims <- unique(strsplit2(names(shifts), '_')[,2])
    if(!all(channels %in% params$resolutions$channel_order)){
      if(all(ims %in% params$resolutions$channel_order)){
        tmp <- ims
        ims <- channels
        channels <- ims
      }else{
        if(standardiseShifts){
          warning('Invalid image names found: will not perform standardiseShifts')
          standardiseShifts = F
        }
      }
    }
  }
  if(standardiseShifts){
    if(verbose){ message('\nStandardising shift across channels and images...') }
    shift_bychannel <- by(shifts, channels, function(x){
      val <- mean(x, na.rm=T)
      val <- round(Re(val)) + (1i * round(Im(val)))
      return(val)
      })
    shift_byim <- by(shifts, ims, function(x){
      val <- mean(x, na.rm=T)
      val <- round(Re(val)) + (1i * round(Im(val)))
      return(val)
    })
    chanshift <- as.complex(shift_bychannel)[match(channels, names(shift_bychannel))]
    imshift <- as.complex(shift_byim)[match(ims, names(shift_byim))]
    if(length(chanshift) != length(imshift) | length(chanshift) != length(shifts)){
      new_shifts <- setNames(chanshift + imshift, names(shifts))
    }
    if(saveShifts){
      write.csv(as.data.frame(shifts), paste0(params$out_dir, 'RAWSHIFTS_', currentFovName, '.csv.gz'))
    }
  }

  if(saveShifts){
    write.csv(as.data.frame(new_shifts), paste0(params$out_dir, 'SHIFTS_', currentFovName, '.csv.gz'))
  }
  params$shifts <<- new_shifts

  ## Calculate intersecting window
  window <- do.call(rbind, lapply(1:length(shifts), function(idx){
    # cat(paste0(idx, ' of ', length(shifts), '...'))
    df <- suppressWarnings( as.data.frame(imager::as.cimg(imList[[idx]])) )
    df$x <- df$x + Re(shifts[idx])
    df$y <- df$y + Im(shifts[idx])
    val <- c(min(df$x), max(df$x), min(df$y), max(df$y))
    return(val)
  }))
  window <- c(max(window[,1]), min(window[,2]), max(window[,3]), min(window[,4]))
  params$intersecting_window <<- window
  if(verbose){ message(paste0('\nIntersecting window is ', paste(window, collapse = ' x '), '...')) }

  if(saveShifts){
    out_file <- paste0(params$out_dir, 'REGIM_', currentFovName, '.png')
    imager::save.image(imager::as.cimg(ref_im), out_file, quality = 1)
  }

  if( returnTroubleShootPlots ){
    new_imList <- setNames( lapply(imList, function(x){
      if(!is.null(chosenZslice)){
        y <- x[,,chosenZslice]
      }else{
        y <- maxIntensityProject(y)
      }
      return(y)
    }), names = names(imList) )
    troubleshootPlots <<- new.env()
    pList <- register_troubleshootPlots(imList = new_imList, ...)
    troubleshootPlots[['REGISTRATION_EVAL']] <<- pList
  }
}
