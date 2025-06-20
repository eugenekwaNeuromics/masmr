## This script contains functions for plotting QC

require(RcppML)
require(tripack)
require(ggplot2)
require(scales)
require(reshape2)

###

getPNNMatrix <- function(

  ## INPUT
  x, ## X coordinates
  y, ## Y coordinates
  label, ## Cell type labels

  ## Delaunay Triangulation Option: fast, and number of times run depends on delaunayNNDegrees
  delaunayTriangulation = T,
  delaunayDistanceThreshold = NULL,
  delaunayNNDegrees = c(1:3),

  ## Distances option: much slower, but only done once
  euclideanDistances = NULL,

  seed = 12345,
  verbose = TRUE,
  ...
){


  ###
  if(length(x) != length(y) | length(x) != length(label) ){
    stop('x, y, and labels must be all of same length!')
  }

  if(!delaunayTriangulation & is.null(euclideanDistances)){
    stop('delaunayTriangulation is FALSE but no euclideanDistances provided!')
  }
  ###

  set.seed( as.integer(seed) )

  ###
  CLB <- as.integer(factor(label))
  OUT <- list()
  all_coords <- x + (y * 1i)
  ###

  ###
  if( is.null(delaunayDistanceThreshold) ){
    if(verbose){ message('delaunayDistanceThreshold not provided...Will approximate...') }

    n_samples = 100
    sampledCx <- sample(all_coords, size=pmin(n_samples, length(all_coords)), replace = F)
    delaunayDistanceThreshold <- mean( sapply(sampledCx, function(cxi){
      distances <- sort( Mod(all_coords - cxi), decreasing = F )
      distances <- distances[1:pmin(100, length(distances))]
      return( mean( distances[-1] ) + 4 * sd( distances[-1] ) )
    }) )
    delaunayDistanceThreshold <- pmax(delaunayDistanceThreshold, 20)

    if(verbose){ message( paste0('Setting delaunayDistanceThreshold to ', delaunayDistanceThreshold, '...')) }
  }
  ###


  ###
  if(delaunayTriangulation){

    if(verbose){ message('Running delaunay triangulation approach for getting nearest neighbours...') }
    dt <- tripack::tri.mesh(
      x=x, y=y, ...
    )
    nns_raw <- tripack::neighbours(dt)
    convex_hull <- tripack::convex.hull(dt)

    if(verbose){ message('\nTabulating...') }

    DTS <- list()
    for(j in 1:length(delaunayNNDegrees)){
      DEGREE <- delaunayNNDegrees[j]
      if(verbose){ message(paste0('Examining neighbours of ', DEGREE, ' degree...')) }

      ## Iteratively add more neighbours
      di = 1
      while(di <= DEGREE){
        if(verbose){ message(paste0('Updating neighbours list with ', di, ' order neighbours...'), appendLF = F)}
        if(di == 1){
          nns <- nns_raw
          di = di + 1
          next
        }
        nns <- lapply(nns, function(x){
          return(unique(unlist(nns_raw[x])))
        })
        di = di + 1
      }
      if(DEGREE != 1){
        nns <- lapply(1:length(nns), function(j){
          return(nns[[j]][(nns[[j]]!=j)])
        })
      }


      if(verbose){ message('\nBuilding frequency matrix...') }
      frequency_matrix <- do.call(rbind, lapply(1:length(nns), function(i){

        # if(verbose){
        #   if(i %% 1000 == 0 | i == 1){
        #     message(paste0(i, ' of ', length(nns), '...'), appendLF = F)
        #   }
        # }

        ## Removal of convex hull
        nnvector <- nns[[i]]
        if(i %in% convex_hull$i){
          nnvector <- nnvector[!(nnvector %in% convex_hull$i)]
        }

        ## Cutting by distance
        if(!is.null(delaunayDistanceThreshold) &
           !is.na(delaunayDistanceThreshold) &
           !is.infinite(delaunayDistanceThreshold)){
          ref_coord <- all_coords[i]
          coords <- all_coords[nnvector]
          dists <- Mod(coords - ref_coord)
          nnvector <- nnvector[dists < delaunayDistanceThreshold]
        }

        values <- tabulate( CLB[nnvector] )
        if(length(values) < max(CLB)){
          values <- c(values, rep(0, max(CLB) - length(values)))
        }

        return(values)
      }))
      # if(verbose){ message('') }

      colnames(frequency_matrix) <-
        paste0('DT', sprintf(paste0('%0', max(nchar(delaunayNNDegrees)), 'd'), DEGREE), '_', levels(factor(label)))
      DTS[[j]] <- frequency_matrix
    }
    DTS <- do.call(cbind, DTS)
    OUT[['DT']] <- DTS
  }
  ###


  ###
  if(!is.null(euclideanDistances)){

    euDists <- sort(unique(euclideanDistances))
    if(verbose){ message('Getting nearest neighbours within specific euclidean distances...') }

    FM <- matrix(0, nrow = length(all_coords), ncol = max(CLB) * length(euDists))

    for(i in 1:length(all_coords)){
      if(verbose){
        if(i %% 1000 == 0 | i == 1){
          message(paste0(i, ' of ', length(all_coords), '...'), appendLF = F)
        }
      }
      dists <- Mod(all_coords - all_coords[i])
      reflab <- CLB[i] + c( (1:length(euDists)) - 1) * max(CLB)
      for(j in 1:length(reflab)){
        FM[,reflab[j]] <- FM[,reflab[j]] + (dists < euDists[j])
      }
    }

    colns <- sapply(1:length(euDists), function(j){
      paste0('EUD',
             sprintf( paste0('%0', max(nchar(euclideanDistances)), 'd'), euDists[j]),
             '_', levels(factor(label)))
    })
    colnames(FM) <- as.vector(colns)
    OUT[['EUD']] <- FM
  }
  ###

  OUT <- do.call(cbind, OUT)
  return(OUT)
}

getPNNLatentFactors <- function(

  ## INPUT
  PNNMatrix,

  ## Expected number of clusters
  nFactors = 4,

  mergeSimilarFactors = TRUE,
  correlationThreshold = 0.1,

  seed = 12345,
  verbose = T,
  ...
){

  ###
  set.seed(as.integer(seed))
  k = as.integer(nFactors)
  OUT <- list()
  ###

  ###
  if(verbose){ message('\nRunning NMF...') }
  nmfd <- RcppML::nmf(
    PNNMatrix,
    k=nFactors,
    verbose = verbose,
    seed = seed,
    ...)
  w <- nmfd$w
  h <- nmfd$h

  ## Min-max scale
  w_scaled <- apply(w, 2, function(x) (x-min(x)) / (max(x) - min(x)) )
  h_scaled <- apply(t(h), 2, function(x) (x-min(x)) / (max(x) - min(x)) )
  if( mergeSimilarFactors ){
    if(verbose){ message('Merging similar factors...') }

    ## Will just cluster using the feature matrix
    corw <- cor(w_scaled)
    hc <- hclust( dist( 1-corw ) )
    clusts <- cutree(hc, h=1-correlationThreshold)

    w_merged <- ( do.call(cbind, by(1:ncol(w_scaled), clusts, function(x){
      if(length(x) > 1){
        y <- rowSums(w_scaled[,x])
      }else{
        y <- w_scaled[,x]
      }
      y = (y - min(y)) / (max(y) - min(y))
      return(y)
    })) )

    h_merged <- ( do.call(cbind, by(1:ncol(h_scaled), clusts, function(x){
      if(length(x) > 1){
        y <- rowSums(h_scaled[,x])
      }else{
        y <- h_scaled[,x]
      }
      y = (y - min(y)) / (max(y) - min(y))
      return(y)
    })) )

    k <- ncol(w_merged)
    if(k != nFactors){
      if(verbose){ message( paste0(
        'Similar factors (>', correlationThreshold,  ' correlation) found...Returning ',
        k, ' instead of ', nFactors, ' factors...'))}
    }else{
      if(verbose){ message('No factor merging needed...Returning ', k, ' factors...')}
    }
    w_scaled <- w_merged
    h_scaled <- h_merged
  }
  h_scaled <- t(h_scaled)

  w_scaled <- data.frame(w_scaled)
  h_scaled <- data.frame(h_scaled)
  rownames(h_scaled) <- colnames(w_scaled) <- paste0('NMF', sprintf(paste0('%0', nchar(k), 'd'), 1:k))
  colnames(h_scaled) <- colnames(PNNMatrix)
  rownames(w_scaled) <- rownames(PNNMatrix)

  OUT[['n_factors']] <- k
  OUT[['point_scores']] <- w_scaled
  OUT[['coefficients']] <- h_scaled

  return(OUT)
}


###

plotQC <- function(
    includeNewBlanks = T,
    synthesisDir = NULL,
    params = get('params', envir = globalenv())
){

  ## Get verbosity
  if(is.logical(params$verbose)){
    verbose = params$verbose
  }else{
    verbose = T
  }
  seed <- params$seed
  if(is.null(seed)){
    seed = 12345L
  }

  if( is.null(synthesisDir) ){
    if(!dir.exists(paste0( params$out_dir ))){
      stop('synthesisDir unspecified and params$out_dir does not exist!')
    }
    synthesisDir <- params$out_dir
  }

  fs <- list.files( synthesisDir, full.names = T )
  if(length(fs)==0){
    stop('No files found! Ensure that synthesiseData() has been run!')
  }
  names(fs) <- gsub('[.]csv|[.]gz|OUT_', '', basename(fs))

  ## Modify params
  params$out_dir <<- gsub('[/][/]', '/', paste0( synthesisDir, '/QC/'))
  if(!dir.exists(params$out_dir)){
    dir.create(params$out_dir)
  }
  if( !exists('ordered_codebook', envir = params) ){
    fx <- list.files(params$parent_out_dir, pattern='ORDEREDCODEBOOK', recursive = T, full.names = T)[1]
    if( length(fx)==0 ){
      stop('ORDEREDCODEBOOK file not found!')
    }
    codebook <- data.table::fread(fx, data.table = F, showProgress = F)
    if(colnames(codebook)[1]=='V1'){
      rownames(codebook) <- codebook[,1]
      codebook[,1] <- NULL
    }
    params$ordered_codebook <<- codebook
  }
  codebook <- params$ordered_codebook
  if(!includeNewBlanks){
    codebook <- codebook[!grepl('^blank-new-\\d+$', tolower(rownames(codebook))),]
  }
  ## Read the first GLOBALCOORD and META_RESOLUTION file
  global_coords <- data.table::fread(fs['GLOBALCOORD'], data.table = F, showProgress = F)
  fx <- list.files(params$parent_out_dir, recursive = T, pattern = 'META_RESOLUTION', full.names = T)[1]
  current_res_info <- try( readLines(fx) )
  if( inherits(current_res_info, 'try-error') ){
    stop('Unable to load a META_RESOLUTION file!')
  }
  n <- current_res_info[c(TRUE, FALSE)]
  val <- current_res_info[c(FALSE, TRUE)]
  val <- lapply(val, function(v) strsplit(v, ' ')[[1]])
  resolutions <- setNames(val, n)
  fov_names <- resolutions$fov_names
  params$fov_names <<- fov_names
  params$global_coords <<- global_coords
  params$resolutions <<- resolutions
  ## For FPKM plot
  if( exists('fpkm_file', envir = params)){
    if( !exists('fpkm', envir = params) | object.size(params$fpkm)==0 ){
      fpkm <- try( data.table::fread(params$fpkm_file, data.table = F, showProgress = F) )
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
  }
  fpkm <- params$fpkm
  ## For nProbes
  if( exists('n_probes_file', envir = params)){
    if( !exists('n_probes', envir = params) | object.size(params$n_probes)==0 ){
      n_probes <- try( data.table::fread(params$n_probes_file, data.table = F, showProgress = F) )
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
  }
  n_probes <- params$n_probes

  ## Load data
  fx <- fs['SPOTCALL_PIXELS']
  if(is.na(fx)){
    stop('No SPOTCALL_PIXELS file found!')
  }
  spotcalldf <- data.table::fread(fx, data.table = F, showProgress = F)
  spotcalldf$CELLNAME[spotcalldf$CELLNAME==''] <- NA
  if(!includeNewBlanks){
    spotcalldf <- spotcalldf[!grepl('^blank-new-\\d+$', tolower(spotcalldf$g)),]
  }
  spotcalldf$spottype <-
    ifelse( !grepl('^blank', tolower(spotcalldf$g)), 'Gene',
            ifelse( grepl('^blank-new-\\d+$', tolower(spotcalldf$g)), 'New blank', 'Blank' ))
  spotcalldf$blank <- spotcalldf$spottype %in% c('Blank', 'New blank')


  ## Plot: correlation of FPKM with n spots
  if(!is.null(fpkm)){

    plotName <- 'FPKMCorrelation'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- data.frame( table( factor(spotcalldf$g, levels=rownames(codebook) ) ))
    df$fpkm <- as.numeric(fpkm[match(df$Var1, names(fpkm))])
    df$spottype <- ifelse( !grepl('^blank', tolower(df$Var1)), 'Gene',
                           ifelse( grepl('^blank-new-\\d+$', tolower(df$Var1)), 'New blank', 'Blank' ))
    df$blank <- df$spottype %in% c('New blank', 'Blank')
    # df$COS <- as.numeric( by(spotcalldf$COS, factor(spotcalldf$g, levels=rownames(codebook)), median) )
    plotTitle <- paste0(
      'Total nCounts: ', nrow(spotcalldf),

      ' | Log Correlation: ', round( cor(log10(1+df$Freq), log10(1+df$fpkm), method = 'pearson'), digits=3),

      '\nCellular %: ', round( 100 * sum(!is.na(spotcalldf$CELLNAME)) / nrow(spotcalldf), digits = 3 ),
      '% (', sum(!is.na(spotcalldf$CELLNAME)), ' counts)',

      ' | Blank %: ', round(100 * sum(df$Freq[df$blank])/sum(df$Freq), digits=3), '% (',
      sum(df$Freq[df$blank]), ' counts)'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot( df ) +
        ggplot2::geom_text( ggplot2::aes(
          x= log10(1+fpkm),
          y= log10(1+Freq),
          label = Var1,
          colour=spottype,
          alpha = blank
        ) ) +
        ggplot2::scale_y_continuous( labels = scales::math_format() ) +
        ggplot2::scale_x_continuous( labels = scales::math_format() ) +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::scale_colour_manual( name='Spot type', values=c('Gene'='black', 'Blank'='red', 'New blank'='blue')) +
        ggplot2::scale_alpha_manual( values=c('TRUE' = 0.25, 'FALSE'=1) ) +
        ggplot2::ylab('Log10 nSpots') + ggplot2::xlab( 'Log10 FPKM') +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5)) +
        ggplot2::guides(alpha='none')
      })

    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = 16, width = 30, units = 'cm')

  }

  ## Plot: Metric distributions
  fx <- fs['SPOTCALL_PIXELS']
  if(!is.na(fx)){

    plotName <- 'CosineDistanceDistribution'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- spotcalldf
    gOrder <- by(spotcalldf$COS, spotcalldf$g, median)
    gOrder <- names(gOrder)[order(gOrder, decreasing = F)]
    plotTitle <- paste0(
      'Distribution of Cosine Distances'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            y= COS,
            x= factor(g, levels = gOrder),
            colour=spottype
          )
        ) +
        ggplot2::geom_violin( fill='transparent' ) +
        ggplot2::geom_boxplot( fill='transparent', width=0.25 ) +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::scale_colour_manual( name='Spot type', values=c('Gene'='black', 'Blank'='red', 'New blank'='blue')) +
        ggplot2::ylab('Cosine distance') + ggplot2::xlab('') +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5),
                       axis.text = ggplot2::element_text(angle=90, hjust=1, vjust=0.5))
    })

    width = length(gOrder) * 0.75
    height = 0.1 * width
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## Spot distribution
  fx <- fs['SPOTCALL_PIXELS']
  if(!is.na(fx)){

    plotName <- 'SpotFOVDistribution'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- spotcalldf
    plotTitle <- paste0(
      'Spot spatial distributions'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot(
          spotcalldf,
          ggplot2::aes(
            y= WY,
            x= WX
          )
        ) +
        ggplot2::facet_wrap( ~spottype ) +
        ggplot2::geom_hex( bins=100 ) +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::scale_x_continuous( limits = c(0, as.numeric(resolutions$xydimensions_pixels[1])) ) +
        ggplot2::scale_y_continuous( limits = c(0, as.numeric(resolutions$xydimensions_pixels[2])) ) +
        ggplot2::ylab('Y pixel coordinate') + ggplot2::xlab('X pixel coordinate') +
        ggplot2::scale_fill_viridis_c(name='nSpots', option='turbo', trans='log10') +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5)) +
        ggplot2::coord_fixed()
    })

    width = 14 * 2
    height = 14
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## Cosine distance distribution
  fx <- fs['SPOTCALL_PIXELS']
  if(!is.na(fx)){

    plotName <- 'CosineSpatialDistribution'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- spotcalldf
    plotTitle <- paste0(
      'Mean cosine distance across spatial bins'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            y= WY,
            x= WX,
            z= COS
          )
        ) +
        ggplot2::facet_wrap( ~spottype ) +
        ggplot2::stat_summary_hex( bins=100 ) +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::scale_x_continuous( limits = c(0, as.numeric(resolutions$xydimensions_pixels[1])) ) +
        ggplot2::scale_y_continuous( limits = c(0, as.numeric(resolutions$xydimensions_pixels[2])) ) +
        ggplot2::ylab('Y pixel coordinate') + ggplot2::xlab('X pixel coordinate') +
        ggplot2::scale_fill_viridis_c(name='Mean Cosine', option='turbo') +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5)) +
        ggplot2::coord_fixed()
    })

    width = 14 * 2
    height = 14
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## Bit detection rate
  fx <- fs['SPOTCALL_PIXELS']
  if(!is.na(fx)){

    plotName <- 'BitDetectionRateAll'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- data.frame( table( factor(spotcalldf$g, levels=rownames(codebook) ) ))
    df$spottype <- ifelse( !grepl('^blank', tolower(df$Var1)), 'Gene',
                           ifelse( grepl('^blank-new-\\d+$', tolower(df$Var1)), 'New blank', 'Blank' ))
    df$blank <- df$spottype %in% c('New blank', 'Blank')
    df <- cbind(df, codebook[match(df$Var1, rownames(codebook)),])
    df <- reshape2::melt(df, measure.vars=colnames(codebook))
    df$variable <- factor(df$variable, levels=colnames(codebook))

    plotTitle <- paste0(
      'Influence of bit status on detection rate (all spots)'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            y= Freq,
            x= factor(value),
            colour=factor(value)
          )
        ) +
        ggplot2::facet_wrap( ~variable, scales='free' ) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::ylab('N spots per gene') + ggplot2::xlab('Bit status') +
        ggplot2::scale_colour_manual(values=c('0'='black', '1'='red')) +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5), legend.position = 'none') +
        ggplot2::scale_y_log10()
    })

    width = 20
    height = 24
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## Bit detection rate for non-blanks
  fx <- fs['SPOTCALL_PIXELS']
  if(!is.na(fx)){

    plotName <- 'BitDetectionRateNonBlank'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- data.frame( table( factor(spotcalldf$g, levels=rownames(codebook) ) ))
    df$spottype <- ifelse( !grepl('^blank', tolower(df$Var1)), 'Gene',
                           ifelse( grepl('^blank-new-\\d+$', tolower(df$Var1)), 'New blank', 'Blank' ))
    df$blank <- df$spottype %in% c('New blank', 'Blank')
    df <- cbind(df, codebook[match(df$Var1, rownames(codebook)),])
    df <- df[!df$blank,]
    df <- reshape2::melt(df, measure.vars=colnames(codebook))
    df$variable <- factor(df$variable, levels=colnames(codebook))

    plotTitle <- paste0(
      'Influence of bit status on detection rate (non-blanks only)'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            y= Freq,
            x= factor(value),
            colour=factor(value)
          )
        ) +
        ggplot2::facet_wrap( ~variable, scales='free' ) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::ylab('N spots per gene') + ggplot2::xlab('Bit status') +
        ggplot2::scale_colour_manual(values=c('0'='black', '1'='red')) +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5), legend.position = 'none') +
        ggplot2::scale_y_log10()
    })

    width = 20
    height = 22
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## Bit error rates
  fx <- fs['SPOTCALL_PIXELS']
  if( !is.na(fx) & any(grepl('^BIT_\\d+$', toupper(colnames(spotcalldf)))) ){

    plotName <- 'BitwiseErrorRates'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- subsetDF(spotcalldf, 'BIT')
    if(ncol(df) != ncol(codebook)){
      stop('Wrong subsetting of spotcall dataframe!')
    }
    df <- as.matrix(df) - as.matrix(codebook)[match(spotcalldf$g, rownames(codebook)),]
    df <- data.frame(df)
    df[,c('g', 'spottype', 'blank')] <- spotcalldf[,c('g', 'spottype', 'blank')]
    df <- reshape2::melt( df, id.vars=c('g', 'spottype', 'blank') )
    df$error_type <- ifelse(df$value==0, 'No error',
                            ifelse(df$value==-1, 'ON to OFF',
                                   ifelse(df$value==1, 'OFF to ON', NA)))

    plotTitle <- paste0(
      'Bitwise error rates'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot( df ) +
        ggplot2::facet_wrap( ~spottype, ncol=1, strip.position = 'right' ) +
        ggplot2::geom_bar( ggplot2::aes(x=as.integer(variable), fill=error_type),
                           position = 'fill' ) +
        ggplot2::theme_minimal(base_size=14) +
        ggplot2::ylab('Error Proportion') + ggplot2::xlab('Bit') +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5) ) +
        ggplot2::scale_fill_manual( name='', values=c('No error' = 'grey', 'OFF to ON' = 'blue', 'ON to OFF' = 'red'))
    })
    height = length(unique(df$spottype)) * 10
    width = 20
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## FOV metrics
  fx <- fs['SPOTCALL_PIXELS']
  if( !is.na(fx) ){

    plotName <- 'FOVMetrics'
    if(verbose){ message( paste0(plotName, '...') ) }

    df <- data.frame(table(spotcalldf$fov))
    df$nblank <- as.numeric(table(factor(spotcalldf$fov[spotcalldf$blank], levels=df$Var1)))
    df$perc_blank <- 100 * df$nblank / df$Freq
    df$medcos <- as.numeric(by(spotcalldf$COS, factor(spotcalldf$fov, levels=df$Var1), median))
    df$x <- as.integer(factor(global_coords[match(df$Var1, global_coords$fov),'x_microns_raw']))
    df$y <- as.integer(factor(global_coords[match(df$Var1, global_coords$fov),'y_microns_raw']))
    measureVars <- c(
      'nSpots' = 'Freq', '% Blank' = 'perc_blank',  'Median Cosine' = 'medcos')

    if( !is.null(fpkm) ){
      df$fpkmcor <- as.numeric(by(
        factor(spotcalldf$g, levels=names(fpkm)),
        factor(spotcalldf$fov, levels=df$Var1), function(x){
          cor( log10(1+table(x)), log10(1+fpkm), method = 'pearson' )
        }))
      measureVars <- c(measureVars, 'Log10 FPKM Correlation' = 'fpkmcor')
    }

    df <- reshape2::melt(df, measure.var = measureVars)
    df$varName <- factor(names(measureVars)[match(df$variable, measureVars)], levels=names(measureVars))
    df$minMax <- as.numeric(unlist(
      by(df$value, factor(df$variable, levels=df$variable[!duplicated(df$variable)]),
         function(x){ (x-min(x)) / (max(x) - min(x)) })))

    plotTitle <- paste0(
      'FOV-wise metrics'
    )

    suppressWarnings({
    p <-
      ggplot2::ggplot( df ) +
      ggplot2::facet_wrap( ~varName, nrow = 2 ) +
      ggplot2::geom_tile( ggplot2::aes(x=x, y=y, fill=minMax), colour='white' ) +
      ggplot2::geom_text( ggplot2::aes(x=x, y=y, label=round(value, digits=4), colour=minMax < 0.25 | minMax > 0.75) ) +
      ggplot2::theme_void(base_size=14) +
      ggplot2::ylab('') + ggplot2::xlab('') +
      ggplot2::ggtitle( plotTitle ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5), legend.position = 'none' ) +
      ggplot2::scale_fill_viridis_c( name='', option='turbo') +
      ggplot2::scale_colour_manual( name='', values=c('TRUE' = 'white', 'FALSE' = 'black')) +
      ggplot2::scale_x_continuous(limits=c(0.5, max(df$x) + 0.5)) +
      ggplot2::scale_y_reverse(limits=c(max(df$y) + 0.5, 0.5)) +
      ggplot2::coord_fixed()
    })
    width = length(measureVars)/2 * sqrt( length(unique(df$Var1)) ) * 2.5
    height =  2 * sqrt( length(unique(df$Var1)) ) * 2.5
    ggplot2::ggsave(paste0(params$out_dir, plotName, '.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## Cellular metrics
  fx <- fs['CELLS']
  if(!is.na(fx)){

    plotName <- 'CellMetrics'
    if(verbose){ message( paste0(plotName, '...') ) }

    cellMeta <- data.table::fread(fx, data.table = F)
    cellMeta <- cellMeta[cellMeta$nCounts>0,]
    cellMeta$fov <- gsub('_\\d+$', '', cellMeta$CELLNAME)
    if( !all(cellMeta$fov %in% params$fov_names) ){
      stop('Unable to parse cell names!')
    }

    plotTitle <- paste0(
      'Cell size spatial distribution'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot( cellMeta, ggplot2::aes(x=Xm, y=Ym, z=nPixels) ) +
        ggplot2::stat_summary_hex(fun = mean, bins=100) +
        ggplot2::theme_void(base_size=14) +
        ggplot2::scale_fill_viridis_c(name = 'Cell Size', option='turbo', trans='log10') +
        ggplot2::theme( legend.position = 'right' ) +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme( plot.title = ggplot2::element_text(hjust=0.5) ) +
        ggplot2::coord_fixed()
    })
    width = 14
    height = 14
    ggplot2::ggsave(paste0(params$out_dir, plotName, '_CellSize.png'), plot = p, height = height, width = width, units = 'cm')

    plotTitle <- paste0(
      'Cell counts spatial distribution'
    )
    suppressWarnings({
      p <-
        ggplot2::ggplot( cellMeta, ggplot2::aes(x=Xm, y=Ym, z=nCounts) ) +
        ggplot2::stat_summary_hex(fun = mean, bins=100) +
        ggplot2::theme_void(base_size=14) +
        ggplot2::scale_fill_viridis_c(name = 'nCounts', option='turbo', trans='log10') +
        ggplot2::theme( legend.position = 'right' ) +
        ggplot2::ggtitle( plotTitle ) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme( plot.title = ggplot2::element_text(hjust=0.5) ) +
        ggplot2::coord_fixed()
    })
    width = 14
    height = 14
    ggplot2::ggsave(paste0(params$out_dir, plotName, '_nCounts.png'), plot = p, height = height, width = width, units = 'cm')
  }

  ## SpatialPatternsPlot
  fx <- fs['SPOTCALL_PIXELS']
  if(!is.na(fx)){

    plotName <- 'SpatialPatterns'
    if(verbose){ message( paste0(plotName, '...') ) }

    nnm <- getPNNMatrix(
      x=spotcalldf$Xm, y=spotcalldf$Ym, label=spotcalldf$g,
      delaunayTriangulation = T,
      delaunayDistanceThreshold = 20,
      delaunayNNDegrees = c(1:3),
      euclideanDistances = NULL,
      seed = seed,
      verbose = verbose
      )

    ## Save the NNM
    if(verbose){ message('Saving PNNMatrix.csv...') }
    data.table::fwrite(
      data.frame(nnm, check.names = F),
      paste0(params$out_dir, 'PNNMatrix.csv'),
      row.names = F, showProgress = F
    )

    nmfResults <- try(getPNNLatentFactors(
      nnm,
      nFactors = 4,
      mergeSimilarFactors = TRUE,
      correlationThreshold = 0.1,
      seed = seed,
      verbose = verbose
    ) )

    if(inherits(nmfResults, 'try-error')){
      warning('Something went wrong with getPNNLatentFactors...Unable to return SpatialPatterns plots...')

    }else{

      df <- cbind( spotcalldf[,c('Xm', 'Ym', 'g', 'blank')], nmfResults[['point_scores']] )

      plotTitle <- paste0('Spatial patterns of gene expression (k=', nmfResults$n_factors, ')')
      suppressWarnings({
        p <-
          ggplot2::ggplot(
            reshape2::melt(df, measure.vars=colnames(nmfResults[['point_scores']])),
            ggplot2::aes(x=Xm, y=Ym, z=value)) +
            ggplot2::stat_summary_hex( fun = mean, bins=100 ) +
            ggplot2::facet_wrap(~variable) +
            ggplot2::scale_y_reverse() +
            ggplot2::scale_fill_viridis_c(name='Score', option='turbo') +
            ggplot2::theme_void(base_size=14) +
            ggplot2::ggtitle(plotTitle) +
            ggplot2::theme( plot.title = ggplot2::element_text(hjust=0.5) ) +
            ggplot2::coord_fixed()
      })
      width = 14
      height = 14
      ggplot2::ggsave(paste0(params$out_dir, plotName, '_SpatialPatterns.png'), plot = p, height = height, width = width, units = 'cm')

      plotTitle <- paste0('Mean membership to each pattern (k=', nmfResults$n_factors, ')')
      df <- df[!df$blank,]
      normdf <- meandf <- cvdf <- list()
      for( gi in sort(unique(as.character(df$g))) ){
        dfsub <- df[df$g==gi, colnames(nmfResults[['point_scores']])]
        dfcv <- rep(0, ncol(dfsub))
        if(nrow(dfsub) > 1){
          dfcv <- apply(dfsub, 2, function(x){
            sum(x>mean(x)) / length(x)
          })
          dfsub <- colMeans(dfsub)
        }
        dfsub <- setNames( as.numeric(dfsub), colnames(nmfResults[['point_scores']]) )
        # dfcv <- setNames(dfcv, colnames(nmfResults[['point_scores']])  )
        meandf[[gi]] <- ( dfsub - min(dfsub) ) / (max(dfsub) - min(dfsub))
        # cvdf[[gi]] <- dfcv
        dfsub_ordered <- sort(dfsub, decreasing = T)
        dfnorm <- dfsub - dfsub_ordered[2]
        dfnorm[dfnorm <0] <- 0
        normdf[[gi]] <- dfnorm
      }
      normdf <- do.call(rbind, normdf)
      gene_order <- c()
      for( ri in 1:ncol(normdf) ){
        subdf <- normdf[ order(normdf[,ri], decreasing=T), ]
        gene_order <- c(gene_order, rownames(subdf)[subdf[,ri] > 0] )
      }
      gene_order <- c( gene_order, rownames(normdf)[!(rownames(normdf) %in% gene_order)] )
      newdf <- reshape2::melt( data.frame('gene' = names(meandf), do.call(rbind, meandf), check.names = F), id.vars = 'gene', value.name = 'mean' )
      # newdf$variability <- reshape2::melt( data.frame('gene' = names(cvdf), do.call(rbind, cvdf), check.names = F), id.vars = 'gene')$value
      newdf$gene <- factor(newdf$gene, levels=gene_order)
      newdf$variable <- factor(newdf$variable, levels=rev(sort(unique(newdf$variable))))

      suppressWarnings({
        p <-
        ggplot2::ggplot(
          newdf, ggplot2::aes(y=variable, x=gene, fill=mean)
        ) +
          ggplot2::geom_raster() +
          ggplot2::xlab('') + ggplot2::ylab('') +
          ggplot2::theme_minimal(base_size=16) +
          ggplot2::theme( axis.text.x = ggplot2::element_text(angle=90, hjust=1, vjust=0.5) ) +
          ggplot2::scale_fill_viridis_c(name='Membership', option='rocket') +
          ggplot2::ggtitle( plotTitle ) +
          ggplot2::theme( plot.title = ggplot2::element_text(hjust=0.5))
      })
      width = length(gene_order) * 0.5 + 2
      height = ncol(nmfResults[['point_scores']]) * 2 + 2
      ggplot2::ggsave(paste0(params$out_dir, plotName, '_GeneMembership.png'), plot = p, height = height, width = width, units = 'cm')
    }
  }

  if(verbose){ message('Done!') }
}
