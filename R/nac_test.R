# CC Test functions -------------------------------------------------------
#' @export
CCTest <- function(A, K, type = 2, z = NULL, y = NULL, ...) {

  if (is.null(z)) z <- spec_clust(A, K, ...)

  if (!is.null(y)) type <- 2

  if(type == 1){
    y <- z
  }else if(type == 2){
    if (is.null(y)) y <- spec_clust(A, K+1, ...)
  }else{
    return(print("wrong test stats type!"))
  }

  Y <- label_vec2mat(y)

  # each row for A%*%Y is multinomial given degree
  if(K == 1 & type == 1) # this only should happen for type 1 test
    result <- list(total=Inf, each=Inf)
  else
    result <- CondChisqTest(X = A%*% Y, z = z, K=K)

  return(c(result, list(z=z,y=y)))
}

CondChisqTest <- function(X, z, K=NULL, d=NULL) {
  #if (is.null(z))   z <- rep(1,dim(X)[1])
  if (is.null(d))   d <- Matrix::rowSums(X) #degree for each node

  # Tstat <- rep(0,K)
  if (is.null(K)){
    Kvec <- sort(unique(z))
    K <- length(Kvec)
  } else {
    Kvec <- 1:K
  }
  n = nrow(X)
  L = ncol(X)
  #if (K==1) return( list(total=Inf, each=Inf) ) # this only should happen for type 1 test

  Tstat <- sapply(Kvec, function(k) {
    idx_k <- z==k;
    if (!any(idx_k)) return(0); # if the community is empty, then the test stats is 0
    # cat('dim(X[idx,] = ',str(X[idx,]), 'dim(X) = ',dim(X),'\n')
    chisq_test(X[idx_k,], d[idx_k])
  } )

  #list(total= sum(Tstat)/sqrt(K), each=Tstat) #separate normalization
  list(total= (sum(Tstat) - n*(L-1))/sqrt(2*n*(L-1)), each=Tstat) #normalize for once
}

chisq_test <- function(X, d) {
  if (is.null(dim(X))) return(0) # assume that X only has one row, hence became a vector
  n <- dim(X)[1]
  L <- dim(X)[2]

  # estimate rho, hence E
  D <- sum(d)
  #cat('D = ',D,', n = ',n,'\n')
  if (D == 0 || n == 1) return(0)
  rho <- Matrix::colSums(X) / D  # assuming all the rows of X have the same prob. vector
  # E <- d %*% t(rho)
  E <- matrix(d,ncol=1) %*% matrix(rho,nrow=1)

  # compute test stat.
  E <- as.vector(E)
  X <- as.vector(X)
  idx <- E > 0;
  if (sum(idx) == 0 || n == 1) return(0)
  Tstat <- sum( (X[idx] - E[idx])^2 / E[idx] ) #normalize for once

  # #normalize each term
  # gamma <-  sqrt((n-1)*(L-1))
  #
  # ((Tstat/gamma) - gamma) / sqrt(2)
}


sampleEveryComm <- function(z, K, ratio=0.5) {
  return(
    as.vector( unlist(sapply(1:K, function(k) {
      index_k <- which(z == k)
      sample(index_k, size = round(length(index_k)*ratio))
    })))
  )
}


# CCsub
# first
#' @export
CCsub <- function(A, K, z=NULL, ratio = 0.5, fromEachCommunity=TRUE,...){
  if (is.null(z)) z <- spec_clust(A, K, ...)

  n <- length(z)
  if (fromEachCommunity) {
    # split nodes clustered into halves
    index1 <- sampleEveryComm(z, K, ratio)
    # index1 <- unlist(sapply(1:K, function(k) {
    #   index_k <- which(z == k)
    #   sample(index_k, size = round(length(index_k)*ratio))
    # }))
    # index1 <- as.vector(index1)
  } else {
    index1 <- sample(n, round(n*ratio))
  }

  index2 <- (1:n)[-index1]
  # y1 <- z[index1]
  # y1[which(y1 == 1)[1:round(length(which(y1 == 1))/2)]] <- K+1

  # perform K+1 clustering on one half
  # if(sum(A[index1, index1]) == 0){
  #   return(list(total = 0))
  # }else{
  y1 <- spec_clust(A[index1, index1], K+1, tau=0.25) # working
  # }
  #y1 <- spec_clust(A[index1, index1], K, tau=0.25) # not working
  #y1 <- z[-index2] # not working

  z2 <- z[index2]
  # CCTest(A[index2, index1], K, z= z2, y = y1)$total
  return( CCTest(A[index2, index1], K, z = z2, y = y1) )
}

# # split nodes first and then cluster with K+1 and K
# CCsub2 <- function(A, K, ratio = 0.5){
#   n <- nrow(A)
#   index1 <- sample(n, size = round(n*ratio))
#   index2 <- (1:n)[-index1]
#   y1 <- spec_clust(A[index1, index1], K+1, tau=0.25) # not working
#   z1 <- spec_clust(A[index2, index2], K, tau=0.25)
#   # z1 <- z[index2]
#   # y1 <- z[index1]
#   CCTest(A[index2, index1], K, z= z1, y = y1)$total
# }

chisq_test_known <- function(X,E) {
  n <- dim(X)[1]
  L <- dim(X)[2]

  E <- as.vector(E)
  X <- as.vector(X)

  idx <- E > 0;

  if (sum(idx) == 0 || n == 1) {
    return(0)
  }

  Tstat <- sum( (X[idx] - E[idx])^2 / E[idx] )
  # gamma <-  sqrt((n-1)*(L-1))
  gamma <-  sqrt(n*(L-1))
  ((Tstat/gamma) - gamma) / sqrt(2)
}
