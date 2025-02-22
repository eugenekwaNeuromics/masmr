## This contains functions for spatial metrics

require(tripack)
require(igraph)

spatialHNNDist <- function(
    spotcalldf,
    seed = 12345,
    params = get('params', envir = globalenv())
){
  codebook <- params$ordered_codebook
  g <- spotcalldf$g
  coords <- spotcalldf$WX + (1i * spotcalldf$WY)
  if(exists('seed', envir=params)){
    seed <- params$seed
  }
  set.seed(seed)

  message('\nCalculating homotypic nearest neighbour distances...')

  hnndist <- rep(NA, nrow(spotcalldf))
  for(i in 1:nrow(codebook)){
    message(paste0(i, ' of ', nrow(codebook), '...'), appendLF = F)
    subcoords <- coords[g==rownames(codebook)[i]]
    if(length(subcoords) <=1 ){ next }
    if(length(subcoords) <=9 ){
      subdists <- rep(0, length(subcoords))
      for(j in 1:length(subcoords)){
        d <- Mod(subcoords - subcoords[j])
        d[d==0] <- NA
        subdists[j] <- min(d, na.rm = T)
      }
    }else{
      dt <- suppressWarnings(tripack::tri.mesh(Re(subcoords), Im(subcoords)))

      if( inherits(dt, 'try-error') ){
        # Cannot build DT, calculate NN the old-fashioned way
        subdists <- rep(0, length(subcoords))
        for(j in 1:length(subcoords)){
          d <- Mod(subcoords - subcoords[j])
          d[d==0] <- NA
          subdists[j] <- min(d, na.rm = T)
        }
      }else{
        nns <- tripack::neighbours(dt)
        subdists <- sapply(1:length(nns), function(j){
          min(Mod(subcoords[nns[[j]]] - subcoords[j]))
        })
      }

    }
    hnndist[g==rownames(codebook)[i]] <- subdists
  }
  return( hnndist )
}

##

spatialIsAdjacent <- function(
    spotcalldf,
    queryLocations,
    gridDistance = 1
){
  ## Check params
  coords <- quer_coords <- NULL
  if(is.data.frame(spotcalldf)){
    coords <- spotcalldf$WX + 1i * spotcalldf$WY
    if( is.null(coords) ){
      stop('Unable to determine which columns to extract coordinate data from!')
    }
  }
  if(is.complex(spotcalldf)){
    coords <- spotcalldf
  }
  if( is.character(queryLocations) & is.data.frame(spotcalldf) ){
    if( !(queryLocations %in% colnames(spotcalldf) & length(queryLocations)==1) ){
      stop('Invalid queryLocations column!')
    }
    bool <- spotcalldf[,queryLocations]
    bool <- as.logical(bool)
    if( any(is.na(bool)) ){
      stop('queryLocations column specified does not return a boolean vector!')
    }
    quer_coords <- (spotcalldf$WX + 1i * spotcalldf$WY)[bool]
  }
  if( is.numeric(queryLocations) & is.data.frame(spotcalldf) ){
    if( !(queryLocations %in% c(1:ncol(spotcalldf)) & length(queryLocations)==1) ){
      stop('Invalid queryLocations column!')
    }
    bool <- spotcalldf[,queryLocations]
    bool <- as.logical(bool)
    if( any(is.na(bool)) ){
      stop('queryLocations column specified does not return a boolean vector!')
    }
    quer_coords <- (spotcalldf$WX + 1i * spotcalldf$WY)[bool]
  }
  if( is.complex(queryLocations) ){
    quer_coords <- queryLocations
  }
  if(is.null(coords) | is.null(quer_coords)){
    stop("Unable to parse coordinate values!")
  }
  if( round(gridDistance) != gridDistance ){
    warning('Rounding gridDistance up to nearest integer')
  }
  gridDistance <- ceiling(gridDistance)
  if(gridDistance < 1){
    stop('gridDistance too small!')
  }

  ## Build grid
  grid <- getRasterCoords( matrix(0, nrow = 2 * gridDistance + 1,  2 * gridDistance + 1) )
  centre <- grid[ceiling(nrow(grid)/2), ceiling(ncol(grid)/2)]
  grid <- as.vector( grid - centre )
  quer_coords <- sapply(quer_coords, function(x){
    return( x + grid )
    })
  maxx <- pmax( max(Re(coords)), max(Re(quer_coords)) )

  quer_idx <- Re(quer_coords) + maxx * (Im(quer_coords) - 1)
  ref_idx <- Re(coords) + maxx * (Im(coords) - 1)

  return(ref_idx %in% quer_idx)
}

##

spatialClusterLeiden <- function(
    spotcalldf,
    leidenResolution = 0.01,
    minNeighbours = 3,
    maxInterSpotDistancePixels = 5,
    seed = 12345,
    distanceMetric = 'COS',
    params = get('params', envir = globalenv())
){
  if( !(distanceMetric %in% colnames(spotcalldf)) ){
    stop('Invalid distanceMetric specified!')
  }
  codebook <- params$ordered_codebook
  if(exists('seed', envir=params)){
    seed <- params$seed
  }
  set.seed(seed)
  rownames(spotcalldf) <- 1:nrow(spotcalldf)


  message('\nClustering homotypic pixels and finding cluster centroids...')


  filtout_idx <- c() #Rownames to drop
  cluster_idx <- rep(0, nrow(spotcalldf)) #Cluster identity per pixel
  for(i in 1:nrow(codebook)){
    message(paste0(i, ' of ', nrow(codebook), '...'), appendLF = F)

    if( sum(spotcalldf$g==rownames(codebook)[i])<=1 ){
      filtout_idx <- c(filtout_idx, rownames(spotcalldf)[spotcalldf$g==rownames(codebook)[i]])
      next
    }
    if( sum(spotcalldf$g==rownames(codebook)[i])<=minNeighbours ){
      filtout_idx <- c(filtout_idx, rownames(spotcalldf)[spotcalldf$g==rownames(codebook)[i]])
      next
    }
    subdf <- spotcalldf[spotcalldf$g==rownames(codebook)[i],]

    dt <- suppressWarnings(tripack::tri.mesh(subdf$WX, subdf$WY))

    if(inherits(dt, 'try-error')){
      ## Unable to build mesh, select the best performer. No cluster ID, no clusterSize
      subdf <- subdf[order(subdf[,distanceMetric], decreasing = F),]
      filtout_idx <- c(filtout_idx, rownames(subdf)[-1])
      next
    }

    ## Build nearest neighbour graph
    nns <- tripack::neighbours(dt)
    nndf <- do.call(rbind, lapply(1:length(nns), function(j){
      return( data.frame('from' = j, 'to' = nns[[j]]) )
    }))
    nndf[,c('from_x', 'from_y')] <- subdf[nndf$from,c('WX', 'WY')]
    nndf[,c('to_x', 'to_y')] <- subdf[nndf$to,c('WX', 'WY')]
    nndf$dist <- Mod( (nndf$to_x + 1i * nndf$to_y) - (nndf$from_x + 1i*nndf$from_y) )
    nndf$weight <- 1/nndf$dist

    ## Only consider links that are close together
    nndf <- nndf[nndf$dist <= maxInterSpotDistancePixels,]
    if( nrow(nndf)==0 ){
      subdf <- subdf[order(subdf[,distanceMetric], decreasing = F),]
      filtout_idx <- c(filtout_idx, rownames(subdf)[-1])
      next
    }
    gra <- igraph::graph_from_data_frame(
      nndf[,c('from', 'to', 'weight')],
      directed = F)
    clusts <- igraph::cluster_leiden(gra, resolution_parameter = leidenResolution)
    clusts <- data.frame('idx' = clusts$names, 'cluster' = clusts$membership)

    ## Pick best performing spot as cluster centroid
    subdf$cluster <- clusts[match(1:nrow(subdf), clusts$idx),'cluster']
    subdf <- subdf[order(subdf[,distanceMetric], decreasing = F),]
    loners <- which(tabulate(subdf$cluster) <= minNeighbours)

    # Save cluster identity
    cluster_idxi <- rep(0, nrow(spotcalldf))
    cluster_idxi[match( rownames(subdf), rownames(spotcalldf) )] <- subdf$cluster
    cluster_idxi[is.na(cluster_idxi)] <- 0
    cluster_idx[cluster_idxi!=0] <- cluster_idxi[cluster_idxi!=0] + max(cluster_idx)

    ## Note down non centroids
    filtout_idx <- c(filtout_idx, rownames(subdf)[
      duplicated(subdf$cluster) |
        is.na(subdf$cluster) | subdf$cluster %in% loners
    ])
  }

  clusterDF <- data.frame(
    'CLUSTER' = cluster_idx,
    'CLUSTER_SIZE' = tabulate(cluster_idx + 1)[cluster_idx + 1],
    'CENTROID' = !(rownames(spotcalldf) %in% filtout_idx)
  )
  clusterDF$CLUSTER_SIZE[clusterDF$CLUSTER==0] <- 1

  return(clusterDF)
}
