# NAC Test functions -------------------------------------------------------

#' NAC test
#'
#' @description  The NAC test to measure the goodness-of-fit of the DCSBM to network data.
#'
#' The function computes the NAC+ or NAC statistics in the paper below.
#'
#' @section References: [Adjusted chi-square test for degree-corrected block models](https://arxiv.org/abs/2012.15047),
#' Linfan Zhang, Arash A. Amini, arXiv preprint arXiv:2012.15047, 2020.
#'
#' @seealso  [snac_test]
#' @param A adjacency matrix
#' @param K number of communities
#' @param z label vector for rows of adjacency matrix. If not given, will be calculated by
#' the spectral clustering
#' @param y label vector for columns of adjacency matrix.
#' @param plus whether or not use column label vector with (`K`+1) communities, default is TRUE
#' @param ... parameters for spectral clustering.
#' @return A list of result
#' \item{stat}{NAC or NAC+ test statistic.}
#' \item{z}{row label vector.}
#' \item{y}{column label vector.}
#' @examples
#' A <- sample_dcpp(500, 10, 4, 0.1)
#' nac_test(A, K = 4)$stat
#' @export
nac_test <- function(A, K, z = NULL, y = NULL, plus = T, ...) {

  if (is.null(z)) z <- spec_clust(A, K, ...)

  if (!is.null(y)) plus <- T

  if (plus){
    if (is.null(y)) y <- spec_clust(A, K+1, ...)
  }else{
    y <- z
  }

  Y <- label_vec2mat(y)

  # each row for A%*%Y is multinomial given degree
  if(K == 1 & (!plus)) # this only should happen for nac test
    result <- list(total=Inf, each=Inf)
  else
    result <- CondChisqTest(X = A%*% Y, z = z, K=K)

  return(list(stat = result, z=z, y=y))
}

#' SNAC test
#'
#' @description The SNAC test to measure the goodness-of-fit of the DCSBM to network data.
#'
#' The function computes the SNAC+ or SNAC statistics in the paper below.
#' The row label vector of the adjacency matrix could be given through `z` otherwise will
#' be estimated by [spec_clust]. One can specify the ratio of nodes used to estimate column
#' label vector. If `plus = TRUE`, the column labels will be estimated by [spec_clust] with
#' (`K`+1) clusters, i.e. performing SNAC+ test, otherwise with `K` clusters SNAC test.
#' One can also get multiple test statistics with repeated random subsampling on nodes.
#'
#' @section References: [Adjusted chi-square test for degree-corrected block models](https://arxiv.org/abs/2012.15047),
#' Linfan Zhang, Arash A. Amini, arXiv preprint arXiv:2012.15047, 2020.
#'
#' @param A adjacency matrix
#' @param K desired number of communities
#' @param z label vector for rows of adjacency matrix. If not provided, will be estimated by
#' [spec_clust]
#' @param ratio ratio of subsampled nodes from the network
#' @param fromEachCommunity whether subsample from each estimated community or the full network,
#' default is TRUE
#' @param plus whether or not use column label vector with (K+1) communities, default is TRUE
#' @param niter number of times the statisitcs are computed
#' @param ... parameters for spectral clustering.
#' @return A list of result
#' \item{stat}{SNAC or SNAC+ test statistic.}
#' \item{z}{row label vector.}
#'
#' @examples
#' A <- sample_dcpp(500, 10, 4, 0.1)
#' snac_test(A, K = 4)$stat
#'
#' @export
snac_test <- function(A, K, z=NULL, ratio = 0.5,
                      fromEachCommunity=TRUE,
                      plus = TRUE,
                      niter = 1, ...){
  if (is.null(z)) z <- spec_clust(A, K, ...)

  n <- length(z)
  stat <- c()
  for (i in 1:niter){
    if (fromEachCommunity) {
      # split nodes clustered into halves
      index1 <- sampleEveryComm(z, K, ratio)
    } else {
      # split the whole network into halves
      index1 <- sample(n, round(n*ratio))
    }

    index2 <- (1:n)[-index1]

    # perform K+1 clustering on one half
    # if(sum(A[index1, index1]) == 0){
    #   return(list(total = 0))
    # }else{
    if (plus){
      y1 <- spec_clust(A[index1, index1], K+1, ...)
    }else{
      y1 <- spec_clust(A[index1, index1], K, ...)
    }
    # }

    z2 <- z[index2]

    stat[i] <- nac_test(A[index2, index1], K, z = z2, y = y1)$stat
  }
  return(list(stat = stat, z = z))

}

#' Estimate community number with SNAC+
#'
#' Applying SNAC+ test sequentially to estimate community number of a network
#' fit to DCSBM
#'
#' @param A adjacency matrix
#' @param Kmin minimum candidate community number
#' @param Kmax maximum candidate community number
#' @param alpha significance level for rejecting the null hypothsis
#' @param labels a matrix with each column being a row label vector for a
#' candidate community number
#' @return estimated community number
#' @seealso [snac_test]
#' @keywords mod_sel
#' @export
snac_select <- function(A, Kmin=1, Kmax, alpha, labels =NULL) {
  if (is.null(labels)) labels <- sapply(Kmin:Kmax, function(k) nett::spec_clust(A, k))
  Tstat <- sapply(Kmin:Kmax, function(k) snac_test(A, k, z = labels[,k-Kmin+1])$stat )

  inrange <- Tstat < qnorm(1-alpha)
  K <- which(inrange)[1]
  # Kmax is not enough, we return Kmax for consistency
  if (is.na(K)) K <- length(Tstat)
  K <- K+Kmin-1

  return(K)
}

CondChisqTest <- function(X, z, K=NULL, d=NULL) {
  if (is.null(d))   d <- Matrix::rowSums(X) #degree for each node

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
  return((sum(Tstat) - n*(L-1))/sqrt(2*n*(L-1))) #normalize for once
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





# Previous functions ------------------------------------------------------

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
