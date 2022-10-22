
#' Spectral clustering (fast)
#'
#' Perform spectral clustering (with regularization) to estimate communities
#'
#' @param A Adjacency matrix (n x n)
#' @param K Number of communities
#' @param type ("lap" | "adj" | "adj2") Whether to use Laplacian or adjacency-based spectral clustering
#' @param tau Regularization parameter for the Laplacian
#' @param nstart argument from function 'kmeans'
#' @param niter argument from function 'kmeans'
#' @param ignore_first_col whether to ignore the first eigen vector when doing spectral clustering
#' @return A label vector of size n x 1 with elements in {1,2,...,K}
#' @keywords comm_detect
#' @export
spec_clust <- function(A, K, type="lap",
                       tau = 0.25, nstart = 20, niter = 10,
                       ignore_first_col = FALSE) {
  # # A should be a sparse matrix
  #
  # if (ignore_first_col && K > 1) {
  #   Uidx = 2:K
  # } else {
  #   Uidx = 1:K
  # }
  # if (type =="adj") {
  #   eig_res <- RSpectra::eigs_sym(A, K)
  #   U <- eig_res$vectors[ , Uidx]
  #
  # } else if (type == "adj2") {
  #   eig_res <- RSpectra::eigs_sym(A, K)
  #   # U <- eig_res$vectors[ , 1:K] %*% diag(abs(eig_res$values))
  #   U <- eig_res$vectors %*% diag(abs(eig_res$values))
  #   U <- U[ , Uidx]
  #
  # } else { # Laplacian-based
  #   n <- dim(A)[1]
  #   #D <- matrix(0, nrow = n[1], ncol = n[2])
  #   #At <- A + tau/n[1]
  #
  #   degs <- Matrix::rowSums(A)
  #   avgDeg <- mean(degs)
  #
  #   alpha0 <- tau*avgDeg
  #   degh <- degs + alpha0
  #
  #   idx <- degh > 0
  #   dtsqi <- degh
  #   dtsqi[idx] <- degh[idx]^(-0.5)
  #
  #   # Dtsqi <- Diagonal(x=dtsqi) # sparse diagonal matrix
  #   # # diag(D)<- (rowSums(At))^(-0.5)
  #   # # L <- D %*% A %*% D
  #   # Dtsqi_A <- Dtsqi %*% A
  #   Dtsqi_A <- dtsqi * A
  #
  #   Ax_fun <- function(x, args=NULL) {
  #     xt <- as.vector(x)*as.vector(dtsqi)
  #     as.numeric(Dtsqi_A %*% xt + alpha0*dtsqi*mean(xt))
  #   }
  #   # eig_res <- arpack(function(x, extra=NULL) {as.vector(as.matrix(L) %*% x)},
  #   #                   options=list(n=dim(A)[1], nev=K, ncv=K+3, which="LM", maxiter=10000),
  #   #                   sym=TRUE, complex=FALSE)
  #   # require(RSpectra)
  #   # eig_res <- igraph::arpack(Ax_fun, sym = TRUE, options = list(n=dim(A)[1], nev=K, ncv=K+3, which="LM", maxiter=10000))
  #   eig_res <- RSpectra::eigs_sym(Ax_fun, K, n = n)
  #   U <- as.matrix(eig_res$vectors[ , Uidx])
  #   # U <- as.matrix(eig_res$vectors)
  # }
  #
  # # if ( K > 1)
  # #   kclust <- kmeans( eig_res$vectors[,2:K], K , nstart = nstart)
  # # else
  # # kclust <- kmeans(U, K, nstart = nstart)
  # # kclust$cluster
  #
  # #return(ClusterR::KMeans_rcpp(U, K, num_init = nstart)$clusters)
  U = spec_repr(A, K, type = type, tau = tau, ignore_first_col = ignore_first_col)
  return(kmeans(U, K, nstart = nstart, iter.max = niter)$cluster)
}

#' Spectral Representation
#'
#' Provides a spectral representation of the network (with regularization)
#' based on the adjacency or Laplacian matrices
#'
#' @param A Adjacency matrix (n x n)
#' @param K Number of communities
#' @param type ("lap" | "adj" | "adj2") Whether to use Laplacian or
#'   adjacency-based spectral clustering
#' @param tau Regularization parameter for the Laplacian
#' @param ignore_first_col whether to ignore the first eigen vector
#' @return The n x K matrix resulting from a spectral embedding of the network into
#'   R^K
#' @keywords net_repr
#' @export
spec_repr = function(A, K, type="lap",
                     tau = 0.25,  ignore_first_col = FALSE) {
  # A should be a sparse matrix

  if (ignore_first_col && K > 1) {
    Uidx = 2:K
  } else {
    Uidx = 1:K
  }
  if (type =="adj") {
    eig_res <- RSpectra::eigs_sym(A, K)
    U <- eig_res$vectors[ , Uidx]

  } else if (type == "adj2") {
    eig_res <- RSpectra::eigs_sym(A, K)
    U <- eig_res$vectors %*% diag(abs(eig_res$values))
    U <- U[ , Uidx]

  } else { # Laplacian-based
    n <- dim(A)[1]
    degs <- Matrix::rowSums(A)
    avgDeg <- mean(degs)

    alpha0 <- tau*avgDeg
    degh <- degs + alpha0

    idx <- degh > 0
    dtsqi <- degh
    dtsqi[idx] <- degh[idx]^(-0.5)
    Dtsqi_A <- dtsqi * A

    Ax_fun <- function(x, args=NULL) {
      xt <- as.vector(x)*as.vector(dtsqi)
      as.numeric(Dtsqi_A %*% xt + alpha0*dtsqi*mean(xt))
    }
    eig_res <- RSpectra::eigs_sym(Ax_fun, K, n = n)
    U <- as.matrix(eig_res$vectors[ , Uidx])
  }

  return(U)
}


# Spectral test based on Jing Lei's paper ---------------------------------

#' Adjusted spectral test
#'
#' @description The adjusted spectral goodness-of-fit test based on Poisson DCSBM.
#'
#' The test is a natural extension on Lei's work of testing goodness-of-fit
#' for SBM. The residual matrix \eqn{\tilde{A}} is computed from the DCSBM estimation
#' expectation of `A`. To speed up computation, the residual matrix uses Poisson variance instead.
#' Specifically,
#'
#' \deqn{
#'  	\tilde{A}_{ij} = (A_{ij} - \hat P_{ij}) / ( n \hat P_{ij})^{1/2}, \quad
#'  	\hat P_{ij} = \hat \theta_i \hat \theta_j \hat B_{\hat{z}_i, \hat{z}_j} \cdot 1\{i \neq j\}
#' }
#'
#' where \eqn{\hat{\theta}} and \eqn{\hat{B}} are computed using [estim_dcsbm] if not provided.
#'
#' @description  Adjusted spectral test
#' @param A adjacency matrix.
#' @param K number of communities.
#' @param z label vector for rows of adjacency matrix. If not given, will be calculated by
#' the spectral clustering.
#' @param DC whether or not include degree correction in the parameter estimation.
#' @param theta give the propensity parameter directly.
#' @param B give the connectivity matrix directly.
#' @param cluster_fct community detection function to get `z` , by default using [spec_clust].
#' @param ... additional arguments for `cluster_fct`.
#' @return Adjusted spectral test statistics.
#'
#' @section References:
#' Details of modification can be seen at [Adjusted chi-square test for degree-corrected block models](https://arxiv.org/abs/2012.15047),
#' Linfan Zhang, Arash A. Amini, arXiv preprint arXiv:2012.15047, 2020.
#'
#' The original spectral test is from [A goodness-of-fit test for stochastic block models](https://projecteuclid.org/euclid.aos/1452004791)
#' Lei, Jing, Ann. Statist. 44 (2016), no. 1, 401--424. doi:10.1214/15-AOS1370.
#'
#' @export
adj_spec_test <- function(A, K, z = NULL, DC = TRUE, theta = NULL, B = NULL,
                          cluster_fct = spec_clust, ...) {
  if (is.null(z)) {
    z <- cluster_fct(A, K, ...)
  }

  if(is.null(theta) | is.null(B)){
    if(DC){
      out <-  estimDCSBM(A,z)
      if(is.null(theta))
        theta <- out$theta
      if(is.null(B))
        B <- out$B
    }else{
      theta <- rep(1, nrow(A))
      B <- estimSBM(A, z)
    }
  }

  idx <- Matrix::rowSums(A) != 0 # remove nodes with zero degree
  A <- A[idx, idx]
  z <- z[idx]
  theta <- theta[idx]

  n <- nrow(A)
  summaryA <- Matrix::summary(A)
  ii <- summaryA$i
  jj <- summaryA$j
  xx <- summaryA$x
  m <- length(ii)
  # compute A / sqrt(P)
  # The calculation can be simplified by only working above diagonal (not implemented yet).
  pp <- theta[ii]*theta[jj]*B[cbind(z[ii],z[jj])]
  xxnew <- xx / sqrt(pp)
  # xxnew <- sapply(1:m, function(r) {
  #   p <- theta[ii[r]]*theta[jj[r]]*B[z[ii[r]],z[jj[r]]]
  #   xx[r] / sqrt(n*p)
  # })
  A_over_stdA <- Matrix::sparseMatrix(c(ii, n), c(jj, n), x = c(xxnew, 0))

  # compute Lx = (A/sqrt(P))*x - sqrt(P)*x
  # note that sqrt(P) = sqrt(\Theta) Z sqrt(B) Z^TRUE sqrt(\Theta)
  # if (version == 1) {
  Z <- label_vec2mat(z)
  Bstd <- sqrt(B)
  sqrt_theta <- sqrt(theta)
  ThetaSqrtZ <- Matrix::Diagonal(x=sqrt_theta) %*% Z
  ThetaSqrtZBstd <- ThetaSqrtZ %*% Bstd
  Lx_fun <- function(x, args=NULL) {
    # as.numeric( A_over_stdA %*% x -  ThetaSqrtZ %*% (Bstd %*% (t(ThetaSqrtZ) %*% x)) )
    # as.numeric( A_over_stdA %*% x -  ThetaSqrtZBstd %*% (t(ThetaSqrtZ) %*% x) )
    as.numeric( A_over_stdA %*% x - ThetaSqrtZBstd %*% t(t(sqrt_theta*x) %*% Z) )
  }

  lam1 <- RSpectra::eigs_sym(Lx_fun, 1, n = n, which="LM")$values
  tstat <- n^(2/3)*(abs(lam1)/sqrt(n) - 2)
  return(tstat)
}
