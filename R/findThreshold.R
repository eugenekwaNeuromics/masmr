## This script contains functions for thresholding given some value.

findThreshold <- function(
    values,
    thresholdMethod,
    ...
){

  if(missing(thresholdMethod)){
    thresholdMethod = thresh_Otsu
  }
  values <- values[!is.na(values) & !is.infinite(values)]
  thresholdChosen <- NA
  if(!is.function(thresholdMethod)){
    stop('Need to supply a function for thresholdMethod!')
  }
  thresholdChosen <- thresholdMethod(values, ...)
  if(is.na(thresholdChosen)){
    warning('Threshold is NA')
  }
  return(thresholdChosen)
}

##

thresh_Otsu <- function(
    values,
    thresholdsToTry = NULL,
    nbreaks = NULL,
    ...
){

  if(is.null(nbreaks)){
    nbreaks <- pmin(1000, round(length(values)/100))
  }
  if(is.null(thresholdsToTry)){
    thresholdsToTry <- seq(
      min(values), max(values),
      length.out=pmax(100, round(nbreaks * 0.1))
      )
  }
  h <- suppressWarnings( hist(values, breaks=nbreaks, plot=F, ...) )
  x <- h$breaks[-length(h$breaks)]
  y <- h$counts
  otsu_res <- sapply(thresholdsToTry, function(t){
    class1 <- y[x <= t]
    class2 <- y[x > t]

    ## If threshold excludes one class, return Inf
    if(length(class1)==0 | length(class2)==0){
      return(Inf)
    }

    w <- length(x)
    w1 <- length(class1) / w
    w2 <- length(class2) / w
    var1 <- var(class1)
    var2 <- var(class2)

    return( w1 * var1 + w2 * var2)
  })
  finalThreshold <- thresholdsToTry[which.min(otsu_res)]
  return(finalThreshold)

}


##

thresh_Elbow <- function(
    values,
    rightElbow = T,
    nbreaks = NULL,
    span = 0.5,
    ...
){
  if(is.null(nbreaks)){
    nbreaks <- pmin(1000, round(length(values)/100))
  }
  h <- suppressWarnings( hist(values, breaks=nbreaks, plot=F, ...) )
  x <- h$breaks[-length(h$breaks)]
  y <- h$counts

  ## Some smoothening with loess
  lmodel <- loess(y~x, span = span, ... )
  yhat <- predict(lmodel)

  ## Find peak
  peaki <- which.max(yhat)
  if( (peaki == length(y)) & rightElbow){
    warning('Right elbow not possible: peak value is on the right! Switching rightElbow to FALSE...')
    rightElbow = F
  }
  if( (peaki == 1) & !rightElbow){
    warning('Left elbow not possible: peak value is on the left! Switching rightElbow to TRUE...')
    rightElbow = T
  }

  if(rightElbow){
    x <- x[peaki:length(x)]
    yhat <- yhat[peaki:length(yhat)]
  }else{
    x <- x[peaki:1]
    yhat <- yhat[1:peaki]
  }

  df <- data.frame(
    'x' =  x,
    'y' = yhat
  )
  lmodel <- lm(y~x, data=df[1,nrow(df)])
  df$pred <- predict(lmodel, new.data = df)
  df$resid <- df$pred - df$y
  finalThreshold <- df$x[which.max(df$y)]
  return(finalThreshold)
}

##

thresh_Quantile <- function(
    values,
    labels,
    quantileFalse = NULL,
    quantileTrue = NULL,
    ...
){
  labels <- as.logical(labels)
  if(sum(labels)==0 | sum(labels)==length(labels)){
    stop('Labels needs to be a boolean with some TRUE values and some FALSE!')
  }

  if(is.null(quantileFalse)){
    quantileFalse = 0.5
  }
  if(is.null(quantileTrue)){
    quantileTrue = 0.5
  }

  threshF <- quantile(values[!labels], quantileFalse, ...)
  threshT <- quantile(values[labels], quantileTrue, ...)
  finalThreshold <- mean(c(threshF, threshT))
  return(finalThreshold)
}



