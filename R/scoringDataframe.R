## This script is for calculating additional metrics

scoreClassifier <- function(
    spotcalldf,
    labels,
    variablesExcluded,
    variablesIncluded,
    trainIndex,
    modelFormula,
    classifyFunction,
    LHS = 'OUT',
    seed = 12345,
    predictType = 'response',
    returnModel = F,
    ...
){
  if(exists('params', envir = globalenv())){
    params <- get('params', envir = globalenv())
    if(exists('seed', envir = params)){ seed = params$seed }
  }
  set.seed(seed)

  ## Labels can be a column name, an index, or a vector
  if(missing(variablesExcluded)){
    variablesExcluded <- c()
  }
  if(length(labels)==1){
    if(is.character(labels)){
      variablesExcluded <- c(variablesExcluded, labels)
    }
    if(is.numeric(labels)){
      variablesExcluded <- c(variablesExcluded, colnames(spotcalldf)[labels])
    }
    labels <- spotcalldf[,labels]
  }
  if(length(labels)!=nrow(spotcalldf)){
    stop('Invalid labels parameter specified!')
  }
  if(!is.logical(labels)){
    lvls <- levels(factor(labels))
    if(length(lvls) == 2 & !any(is.na(as.logical(labels)))){
      labels <- as.logical(labels)
    }else{
      labels <- factor(labels)
    }
  }
  if( sum(labels)==0 | sum(labels)==length(labels) ){
    stop('Unable to separate labels into two classes!')
  }
  if(LHS %in% colnames(spotcalldf)){
    warning('Overwriting the ', LHS, ' column in spotcalldf!')
  }
  spotcalldf[[LHS]] <- labels

  ## Default is logistic regression
  if(missing(classifyFunction)){
    classifyFunction <- glm
    formals(classifyFunction)$family <- 'binomial'
  }

  ## If no training index specified, train against all
  if(missing(trainIndex)){
    trainIndex <- T
  }

  ## Check which covariates have missing, non.numeric, or Inf values
  if(missing(modelFormula)){
    message('\nChecking which covariates to exclude...')
    drop_covar <- c(LHS, variablesExcluded)
    for(i in 1:ncol(spotcalldf)){
      if(colnames(spotcalldf)[i] %in% drop_covar){ next }
      if( all(spotcalldf[,i] == spotcalldf$OUT) ){
        drop_covar <- c(drop_covar, colnames(spotcalldf)[i])
        next
      }
      if(is.character(spotcalldf[,i])){
        drop_covar <- c(drop_covar, colnames(spotcalldf)[i])
        next
      }
      nacheck <- sum(is.na(spotcalldf[,i]))
      infcheck <- sum(is.infinite(spotcalldf[,i]))
      if(nacheck > 0 | infcheck > 0){
        drop_covar <- c(drop_covar, colnames(spotcalldf)[i])
        next
      }
      if(var(as.numeric(spotcalldf[,i]))==0){
        drop_covar <- c(drop_covar, colnames(spotcalldf)[i])
        next
      }
    }
    if(missing(variablesIncluded)){
      variablesIncluded <- colnames(spotcalldf)
    }
    chosen_covars <- variablesIncluded[!(variablesIncluded %in% drop_covar)]

    ## Build formula
    modelFormula <- paste0(
      LHS,
      ' ~ ',
      paste(chosen_covars, collapse = ' + ')
    )
  }else{
    formsplit <- strsplit2(modelFormula, '~')
    if(ncol(formsplit) != 2){
      stop('Invalid formula format!')
    }
    modelFormula <- paste0(LHS, ' ~ ', formsplit[,2])
  }
  message( paste0('\nFormula: ', modelFormula ))
  t <- try( classifyFunction(as.formula(modelFormula), data = spotcalldf[trainIndex,], ...), silent = T )
  if( inherits(t, 'try-error') ){
    modelFormula <<- modelFormula
    stop('classifyFunction error! Returning formula instead (modelFormula)...')
  }

  ## Build model
  lmodel <- t
  if( returnModel ){ classifierModel <<- lmodel }

  ## Predict
  t <- try( predict(lmodel, newdata = spotcalldf, type = predictType ), silent = T )
  if(inherits(t, 'try-error')){
    classifierModel <<- lmodel
    stop('predict function incompatible with model! Returning model instead (classifierModel)...')
  }
  score <- t
  if(length(score) != nrow(spotcalldf)){
    warning('Vector length does not match nrows of the dataframe, but returning vector regardless...')
  }

  return(score)
}

##

scoreDeltaF <- function(
    spotcalldf,
    prefixes,
    params = get('params', envir = globalenv())
){
  ## Get required data
  if(missing(prefixes)){
    prefixes <- params$imageMetrics
  }
  if(is.null(prefixes)){
    stop('Prefixes for calculating dF are unspecified!')
  }
  codebook <- params$ordered_codebook
  if(is.null(codebook)){
    stop('Ordered codebook not found!')
  }
  if(!('g' %in% colnames(spotcalldf))){
    stop('spotcalldf$g column not found! Ensure decode has been run!')
  }
  geneMatrix <- codebook[spotcalldf$g,]
  if( any(dim(geneMatrix) != c(nrow(spotcalldf), ncol(codebook))) ){
    stop('Unable to generate geneMatrix from spotcalldf$g!')
  }

  message('\nCalculating delta F metrics...')
  for(i in 1:length(prefixes)){
    prefixi <- prefixes[i]
    message(paste0(prefixi, ' ( ', i, ' of ', length(prefixes), ' )...'), appendLF = F)
    NORM <- subsetDF(spotcalldf, prefixi)
    if(ncol(NORM) != ncol(codebook)){
      warning( paste0('Incorrect number of columns subset for ', prefixi, '...Skipping...') )
      next
    }
    gM <- geneMatrix #ON=1, OFF=0
    fon <- rowSums(NORM * gM) / rowSums(gM) #Mean
    fonVar <- rowSums((NORM - fon)^2 * gM) / (rowSums(gM)-1) #Variance
    gM <- 1-geneMatrix #ON=0, OFF=1
    foff <- rowSums(NORM * gM) / rowSums(gM)
    foffVar <- rowSums((NORM - foff)^2 * gM) / (rowSums(gM)-1)
    dffo <- (fon-foff)/foff
    dffo[is.infinite(dffo) | foff==0] <- max(dffo[!(is.infinite(dffo) | foff==0)], na.rm = T)
    fres <- data.frame(
      'F1' = fon,
      'F0' = foff,
      'DF' = (fon-foff),
      'DFF0' = dffo, #When F0 is 0, maximum non-infinite value assigned
      'F1VAR' = fonVar,
      'F0VAR' = foffVar,
      'F1CV' = sqrt(fonVar) / fon,
      'F0CV' = sqrt(foffVar) / foff
    )
    colnames(fres) <- paste0(colnames(fres), '_', prefixi)
    if(nrow(fres) != nrow(spotcalldf)){
      stop('Number of rows not consistent across dataframes!')
    }
    if( any(colnames(fres) %in% colnames(spotcalldf)) ){
      spotcalldf[,colnames(fres)] <- fres
    }else{
      spotcalldf <- cbind(spotcalldf, fres)
    }
  }
  return(spotcalldf)
}

##

scoreBlanks <- function(
    spotcalldf,
    criteria = c(
      'gBlank',
      'gbPossible',
      'consistentLabels'
      ),
    distanceMetrics = c('COS', 'EUC', 'L2E'),
    params = get('params', envir = globalenv())
){

  ## Checks that required data to reference is there
  criteria <- toupper(criteria)
  valid_criteria <- c(
      'gBlank',
      'bPossible',
      'gbPossible',
      'consistentLabels'
    )
  valid_criteria <- toupper(valid_criteria)
  if(any(!(criteria %in% valid_criteria))){
    stop('Invalid criteria specified!')
  }
  valid_distances <- c('COS', 'EUC', 'L2E')
  distanceMetrics <- toupper(distanceMetrics)
  if(any(!(distanceMetrics %in% valid_distances))){
    stop('Invalid distance metric specified!')
  }
  if( !all(distanceMetrics %in% colnames(spotcalldf)) ){
    stop('One or more distanceMetrics not in colnames of dataframe!')
  }
  codebook <- params$ordered_codebook
  if(is.null(codebook)){
    stop('params$ordered_codebook not found: make sure prepareCodebook() has been run!')
  }
  if( !('g' %in% colnames(spotcalldf)) ){
    stop('g column not found in dataframe: ensure decode() has been run!')
  }

  pixelIsBlank <- rep(0, nrow(spotcalldf))

  # Check: that 'blank' is in the name of the gene
  if( toupper('gBlank') %in% criteria){
    message('gBlank test...')
    # params$isblank from perpareCodebook
    x <- spotcalldf$g %in% rownames(params$ordered_codebook)[params$isblank]
    pixelIsBlank = pixelIsBlank + as.integer(x)
  }

  # Check: that next closest entry is possible for all distance metrics in distanceMetrics
  if( toupper('bPossible') %in% criteria){
    message('bPossible test...')
    labels <- spotcalldf[,paste0(distanceMetrics, 'LAB')]
    blabels <- spotcalldf[,paste0('B', distanceMetrics, 'LAB')]

    # Find which values are possible for a given blank
    if(!exists('cognates')){
      cognates <- list()
      for(i in 1:nrow(codebook)){
        ref <- do.call(rbind, lapply(1:nrow(codebook), function(x) return(codebook[i,])))
        ref <- codebook - ref
        dists <- rowSums(abs(ref))
        dists[dists==0] <- NA
        cognates[[i]] <- which(dists==min(dists, na.rm = T))
      }
    }

    # Check if each gene has a cognate neighbour
    y <- rep(0, nrow(spotcalldf))
    for(i in 1:length(distanceMetrics)){
      message(paste0('For ', distanceMetrics[i], '...' ), appendLF = F)
      labi <- labels[,i]
      blabi <- blabels[,i]
      for(j in 1:length(cognates)){
        x <- (labi==j) & !(blabi %in% cognates[[j]])
        y = y + as.integer(x)
      }
    }
    message('')
    pixelIsBlank = pixelIsBlank + as.integer(y > 0)
  }

  # Check: similar to above, but for one distance metric only
  if( toupper('gbPossible') %in% criteria){
    message('gbPossible test...')
    labs <- spotcalldf[,grepl('LAB$', colnames(spotcalldf))]
    labs <- labs[,!grepl('^B', colnames(labs))]
    gvector <- match(spotcalldf$g, rownames(codebook))
    domDM <- colnames(labs)[sapply(1:ncol(labs), function(i){
      all(labs[,i] == gvector)
    })]
    if(length(domDM)==0){
      warning('Unable to find original g column: gbPossible test cannot be performed...Skipping...')
    }
    if(length(domDM)>1){
      if('COSLAB' %in% domDM){
        domDM <- 'COSLAB'
      }
      domDM <- domDM[1]

      message(paste0('Using ', domDM, '...'), appendLF = F)
      message('')
      labi <- spotcalldf[,domDM]
      blabi <- spotcalldf[,paste0('B',domDM)]

      # Find which values are possible for a given blank
      if(!exists('cognates')){
        cognates <- list()
        for(i in 1:nrow(codebook)){
          ref <- do.call(rbind, lapply(1:nrow(codebook), function(x) return(codebook[i,])))
          ref <- codebook - ref
          dists <- rowSums(abs(ref))
          dists[dists==0] <- NA
          cognates[[i]] <- which(dists==min(dists, na.rm = T))
        }
      }

      # Check if each gene has a cognate neighbour
      y <- rep(0, nrow(spotcalldf))
      for(j in 1:length(cognates)){
        x <- (labi==j) & !(blabi %in% cognates[[j]])
        y = y + as.integer(x)
      }
      pixelIsBlank = pixelIsBlank + as.integer(y > 0)
    }
  }

  # Check: consistent labels
  if( toupper('consistentLabels') %in% criteria ){
    if(length(distanceMetrics) == 1){
      warning('Cannot perform consistentLabels test if only one distance metric specified!')
    }else{
      message('consistentLabels test...')
      labels <- spotcalldf[,paste0(distanceMetrics, 'LAB')]
      x <- rowSums( apply(labels, 2, function(li){
        li!=labels[,1]
      }) )>0
      pixelIsBlank = pixelIsBlank + as.integer(x)
    }
  }

  return(pixelIsBlank)
}
