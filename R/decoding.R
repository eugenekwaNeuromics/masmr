## Code for decoding

require(parallel)
require(Rcpp)
require(RcppEigen)
require(RcppArmadillo)

## For fast cosine distance calculations
# Rcpp::cppFunction('
#   SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
#                       Eigen::Map<Eigen::MatrixXd> B,
#                       int n_cores){
#   Eigen::setNbThreads(n_cores);
#   Eigen::MatrixXd C = A * B;
#   return Rcpp::wrap(C);
#   }', depends = c('RcppArmadillo', 'RcppEigen')  )

##

subsetDF <- function(
    spotcalldf,
    prefix
){
  NORM <- as.matrix(spotcalldf[,grepl(paste0('^', prefix, '_\\d+$'), colnames(spotcalldf))])
  return(NORM)
}

##

decode <- function(
    spotcalldf,
    nCores,
    currentFovName,
    decodeMetric = 'DECODE',
    params = get('params', envir = globalenv())
){
  suppressMessages(
    Rcpp::cppFunction('
      SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                          Eigen::Map<Eigen::MatrixXd> B,
                          int n_cores){
      Eigen::setNbThreads(n_cores);
      Eigen::MatrixXd C = A * B;
      return Rcpp::wrap(C);
      }', depends = c('RcppEigen'), verbose = F  )
  )

  if(missing(nCores)){
    nCores = pmax(1, parallel::detectCores() - 1)
  }

  if(missing(currentFovName)){
    currentFovName <- params$current_fov
  }
  if(is.null(currentFovName)){
    stop('currentFovName NOT specified and params$current_fov does not exist!')
  }


  if( !(decodeMetric %in% params$imageMetrics) ){
    stop('decodeMetric absent from params$imageMetrics!')
  }
  NORM <- subsetDF(spotcalldf, decodeMetric)

  ## Remove values where decoding is not possible
  filtout <- rowSums(NORM) == 0
  spotcalldf <- spotcalldf[!filtout,]
  NORM <- NORM[!filtout,]
  codebook <- params$ordered_codebook
  nbits <- params$nbits
  if(is.null(codebook)){
    stop("Unable to find ordered codebook!")
  }

  ## Generate metrics for decoding
  message('Decoding...')
  l2norm <- (NORM / rowSums(NORM))
  tc <- t(apply(codebook, 2, as.numeric))
  nA <- sqrt(rowSums(NORM^2))
  B <- sqrt(colSums(tc^2))
  nextbestL <- nextbestE <- nextbestN <- bestL <- bestE <- bestN <- rep(Inf, nrow(spotcalldf))
  nlabL <- nlabE <- nlabN <- labL <- labE <- labN <- rep(0, nrow(spotcalldf))
  for(idx in 1:ncol(tc)){
    message(paste0(idx, ' of ', ncol(tc), '...'), appendLF = F)
    nscores <- eigenMapMatMult2(NORM, matrix(tc[,idx]), n_cores = nCores) #Dot product
    nscores <- 1 - (as.numeric(nscores) / (nA * B[idx])) #Cosine distance (arg min). This is half of the L2 Norm Euclidean described in Moffit.
    ndiff <- bestN - nscores
    nextbestN[ndiff>0] <- bestN[ndiff>0]
    nlabN[ndiff>0] <- labN[ndiff>0]
    bestN[ndiff>0] <- nscores[ndiff>0]
    labN[ndiff>0] <- idx

    escores <- colSums((t(NORM) - tc[,idx])^2) #Raw euclidean distance (arg min).
    ediff <- bestE - escores
    nextbestE[ediff>0] <- bestE[ediff>0]
    nlabE[ediff>0] <- labE[ediff>0]
    bestE[ediff>0] <- escores[ediff>0]
    labE[ediff>0] <- idx

    lscores <- colSums((t(l2norm) - (tc[,idx]/nbits) )^2) #L2Norm euclidean distance (arg min).
    lscores <- sqrt(lscores)
    ldiff <- bestL - lscores
    nextbestL[ldiff>0] <- bestL[ldiff>0]
    nlabL[ldiff>0] <- labL[ldiff>0]
    bestL[ldiff>0] <- lscores[ldiff>0]
    labL[ldiff>0] <- idx
  }
  spotcalldf <- data.frame(
    spotcalldf,

    'COS' = bestN,
    'COSLAB' = labN,
    'EUC' = bestE,
    'EUCLAB' = labE,
    'L2E' = bestL,
    'L2ELAB' = labL,

    'BCOS' = nextbestN,
    'BCOSLAB' = nlabN,
    'BEUC' = nextbestE,
    'BEUCLAB' = nlabE,
    'BL2E' = nextbestL,
    'BL2ELAB' = nlabL,

    check.names = F
  )

  gcx <- params$global_coords
  gcx <- gcx[gcx$fov == currentFovName,]
  gcx <- as.numeric(unlist(gcx[!duplicated(gcx$fov),c('x_microns', 'y_microns')]))
  per_pixel_microns <- as.numeric(params$resolutions$per_pixel_microns)
  if( !is.null(gcx) & length(gcx)==2 & !is.null(per_pixel_microns) ){
    globalcoorddf <- data.frame(
      'Xm' = spotcalldf$WX * per_pixel_microns[1] + gcx[1],
      'Ym' = spotcalldf$WY * per_pixel_microns[2] + gcx[2],
      'fov' = currentFovName,
      'g' = rownames(codebook)[spotcalldf$COSLAB]
    )
  }else{
    warning('Unable to parse global coordinate information: returning Xm and Ym as NA for now')
    globalcoorddf <- data.frame(
      'Xm' = NA,
      'Ym' = NA,
      'fov' = currentFovName,
      'g' = rownames(codebook)[spotcalldf$COSLAB]
    )
  }
  spotcalldf <- cbind(globalcoorddf, spotcalldf)

  return(spotcalldf)
}


###

bitsFromDecode <- function(
    spotcalldf,
    trainIndex,
    metricToMaximise = 'f1',
    fBeta = 2,
    nThresholdTests = 100,
    saveThresholds = T,
    decodeMetric = 'DECODE',
    params = get('params', envir = globalenv())
){
  ## Set up parameters
  if(!is.numeric(fBeta)){
    stop('Invalid fBeta value!')
  }
  if(length(metricToMaximise) > 1){
    warning('Considering only the first value in metricToMaximise')
  }
  metricToMaximise <- tolower(metricToMaximise[1])
  valid_metrics <- c('precision', 'recall', 'f1', 'fb')
  if( !(metricToMaximise %in% valid_metrics) ){
    stop('Invalid metricToMaximise!')
  }
  if( !all(c('g', 'fov') %in% colnames(spotcalldf)) ){
    stop('Ensure decode() has been run first!')
  }
  if( !(decodeMetric %in% params$imageMetrics) ){
    stop('decodeMetric absent from params$imageMetrics!')
  }
  if( any(grepl('^BIT_\\d+', colnames(spotcalldf))) ){
    warning('Will be overwriting BIT_XX columns in dataframe')
  }
  g <- spotcalldf$g
  codebook <- params$ordered_codebook
  NORM <- subsetDF(spotcalldf, decodeMetric)
  CALL <- data.frame( matrix(0, nrow=nrow(NORM), ncol=ncol(NORM)) )
  colnames(CALL) <- paste0(
    'BIT_',  sprintf(paste0('%0', nchar(ncol(CALL)), 'd'), 1:ncol(CALL)))
  currentFovName <- unique( spotcalldf$fov )[1]
  if(ncol(NORM) != ncol(codebook)){
    stop('Invalid subsetting of dataframe! Check decodeMetric_XXX is specific!')
  }
  if(missing(trainIndex)){
    trainIndex <- T
  }

  ## Initialise
  message('Finding optimal intensity thresholds per bit...')


  ## Recommend thresholds
  thresholds <- rep(0, ncol(codebook))
  for(i in 1:ncol(codebook)){
    message( paste0(i, ' of ', ncol(codebook), '...'), appendLF = F )
    expectedBit <- codebook[g,i]
    intensityVals <- NORM[,i]

    ## Make sure valid numbers
    expectedBit <- expectedBit[!is.na(intensityVals) & !is.infinite(intensityVals)]
    intensityVals <- intensityVals[!is.na(intensityVals) & !is.infinite(intensityVals)]

    ## Subset values
    Y <- expectedBit[trainIndex]
    X <- intensityVals[trainIndex]

    ## Goal is to find threshold for X that maximises accuracy to Y
    ## Owing to possible class imbalance, will opt for precision-recall, and F1 metrics
    thresholdsToTry <- seq(min(X), max(X), length.out = nThresholdTests)
    scores <- sapply(thresholdsToTry, function(tx){
      Yhat <- as.integer(X > tx) #1 if brighter, 0 if not
      TP <- sum( (Y+Yhat)==2 ) #True positives
      TN <- sum( (Y+Yhat)==0 ) #True negatives
      FP <- sum( ((Yhat==1)+(Y==0))==2 ) #False positives
      FN <- sum( ((Yhat==0)+(Y==1))==2 ) #False negatives
      precision <- TP / (TP + FP)
      recall <- TP / (TP + FN)

      if(metricToMaximise == 'precision'){
        return(precision)
      }
      if(metricToMaximise == 'recall'){
        return(recall)
      }
      if(metricToMaximise == 'f1'){
        f1 <- (2 * precision * recall) / ( precision + recall )
        return(f1)
      }
      if(metricToMaximise == 'fb'){
        fb <- ( (1 + fBeta^2) * precision * recall ) /
          ( (fBeta^2 * precision) + recall )
        return(fb)
      }
      return(NA)
    })

    thresholdsToTry <- thresholdsToTry[!is.na(scores)]
    scores <- scores[!is.na(scores)]
    if(length(scores)==0){
      stop('Unable to calculate suitable threshold!')
    }
    threshx <- thresholdsToTry[which.max(scores)]
    thresholds[i] <- threshx

    CALL[,i] <- as.integer(NORM[,i] > threshx)
  }
  spotcalldf[,colnames(CALL)] <- CALL
  spotcalldf$HD <- rowSums( abs(CALL - codebook[g,]) > 0)

  if(saveThresholds){
    write.csv( data.frame(thresholds),
               paste0(params$out_dir,
                      'THRESHBITCALL_',
                      toupper(metricToMaximise), '_',
                      toupper(decodeMetric), '_',
                      currentFovName, '.csv'),
               row.names = F)
  }

  return(spotcalldf)
}


