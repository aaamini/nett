# Spectral clustering -----------------------------------------------------
# Spectral clustering with perturbation
#' @export
spec_clust <- function(A, K, type="lap",
                       tau = 0.25, nstart = 20, niter = 10,
                       ignore_first_col = F) {
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
    # U <- eig_res$vectors[ , 1:K] %*% diag(abs(eig_res$values))
    U <- eig_res$vectors %*% diag(abs(eig_res$values))
    U <- U[ , Uidx]

  } else { # Laplacian-based
    n <- dim(A)[1]
    #D <- matrix(0, nrow = n[1], ncol = n[2])
    #At <- A + tau/n[1]

    degs <- Matrix::rowSums(A)
    avgDeg <- mean(degs)

    alpha0 <- tau*avgDeg
    degh <- degs + alpha0

    idx <- degh > 0
    dtsqi <- degh
    dtsqi[idx] <- degh[idx]^(-0.5)

    # Dtsqi <- Diagonal(x=dtsqi) # sparse diagonal matrix
    # # diag(D)<- (rowSums(At))^(-0.5)
    # # L <- D %*% A %*% D
    # Dtsqi_A <- Dtsqi %*% A
    Dtsqi_A <- dtsqi * A

    Ax_fun <- function(x, args=NULL) {
      xt <- as.vector(x)*as.vector(dtsqi)
      as.numeric(Dtsqi_A %*% xt + alpha0*dtsqi*mean(xt))
    }
    # eig_res <- arpack(function(x, extra=NULL) {as.vector(as.matrix(L) %*% x)},
    #                   options=list(n=dim(A)[1], nev=K, ncv=K+3, which="LM", maxiter=10000),
    #                   sym=T, complex=F)
    # require(RSpectra)
    # eig_res <- igraph::arpack(Ax_fun, sym = T, options = list(n=dim(A)[1], nev=K, ncv=K+3, which="LM", maxiter=10000))
    eig_res <- RSpectra::eigs_sym(Ax_fun, K, n = n)
    U <- as.matrix(eig_res$vectors[ , Uidx])
    # U <- as.matrix(eig_res$vectors)
  }

  # if ( K > 1)
  #   kclust <- kmeans( eig_res$vectors[,2:K], K , nstart = nstart)
  # else
  # kclust <- kmeans(U, K, nstart = nstart)
  # kclust$cluster

  #return(ClusterR::KMeans_rcpp(U, K, num_init = nstart)$clusters)
  return(kmeans(U, K, nstart = nstart, iter.max = niter)$cluster)
}

# Spectral test based on Jing Lei's paper ---------------------------------
# new fast spectral test
#' @export
spec_test <- function(A, K, z = NULL, theta = NULL, B = NULL, DCSBM = T) {
  if (is.null(z)) {
    z <- spec_clust(A, K)
  }

  if(is.null(theta) | is.null(B)){
    if(DCSBM){
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
  # note that sqrt(P) = sqrt(\Theta) Z sqrt(B) Z^T sqrt(\Theta)
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
  # ThetaZ <- diag(theta) %*% Z
  # P  <- ThetaZ %*% (B %*% t(ThetaZ))
  # L <- (A-P)/sqrt(P)

  # AoverSqrtP <- A/sqrt(P)
  # mean(abs(A_over_stdA -AoverSqrtP))
  # mean(abs(ThetaSqrtZ %*% (Bstd %*% (t(ThetaSqrtZ))) - sqrt(P)))
  # vecloss(L %*%x, Lx_fun(x))

  # } else {
  #   #Bstd <- sqrt(B)
  #   sqrt_theta <- sqrt(theta)
  #   sqrtThetaB <- sqrt(theta*B[z,])
  #   #sqrtThetaB <- sqrt_theta * Bstd[z,]
  #   Lx_fun <- function(x, args=NULL) {
  #     # as.numeric( A_over_stdA %*% x - (sqrt_theta * Bstd[z,]) %*% aggregate(x*sqrt_theta, list(label=z), sum)$x )
  #     as.numeric( A_over_stdA %*% x - sqrtThetaB %*% aggregate(x*sqrt_theta, list(label=z), sum)$x )
  #   }
  #
  # }

  # need to be careful determining the variance of A
  # it is likely to get estimated P greater than 1
  # which can be a problem when variance is p*(1-p)
  # L <- (A-P)/sqrt(n*(1-P)*P)

  # # require(RSpectra)
  # # don't forget to divide sqrt(n)
  # lam1 <- eigs_sym(Lx_fun, 1, n = n, which="LA")$values/sqrt(n)
  # lamn <- eigs_sym(Lx_fun, 1, n = n, which="SA")$values/sqrt(n)
  # sig <- max(abs(lam1), abs(lamn))
  # return(c(lam1, lamn, n^(2/3) *(sig - 2)))

  # lam1 <- eigs_sym(L, 1, which="LM")$values
  lam1 <- RSpectra::eigs_sym(Lx_fun, 1, n = n, which="LM")$values
  tstat <- n^(2/3)*(abs(lam1)/sqrt(n) - 2)
  return(tstat)
}
