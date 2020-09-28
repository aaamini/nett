
# AA: Functions are mostly provided by Aiyou Chen (most likely!)

#' @export
pp_conn <- function(n, p, lambda, csizes) {
  # create a planted partition connectivity matrix
  # the P matrix : block density
  # lambda is the desired avg degree, n is size, p is cout/cint, csizes is pi
  if (sum(csizes) != 1) {
    csizes = csizes/sum(csizes)
  }
  csizes = csizes * n
  Kc = length(csizes)
  P0 = p + diag(rep(1 - p, Kc))
  # cat('cout/cin = ', p, '\n')
  P0 = round(pmax(P0, t(P0))/max(P0), 2)
  scale = sum(diag(csizes) %*% P0 %*% diag(csizes)/sum(csizes))
  Pmat = pmin(P0 * lambda/scale, 1)
  return(Pmat)
}

#' @export
quickDCSBM <- function(n, lambda, K, oir,  theta=NULL) {
  pr <- (1:K)/K
  B <- pp_conn(n, oir, lambda, pr)
  if (is.null(theta)) theta <- EnvStats::rpareto(n, 2/3, 3)
  z <- sample(K, n, replace=T, prob=pr)

  list(adj=fastDCSBM(z, B, theta=theta), labels=z, B=B)
}

#' @export
fastDCSBM <- function(z, Pmat, theta=1) {
  csizes = tabulate(z)
  n = length(z)
  if (length(theta) == 1) theta = rep(theta,n)
  A <- fastDCSBM.internal(csizes, Pmat, theta)
  #A = fastSBM.internal(csizes, Pmat)  # E[A_{ij}] = B_{z[sig[i]], z[sig[j]]}
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig) # or order(sig)
  return(A[siginv,siginv]) # E[A_{siginv[i], siginv[j]}] = B_{z[i], z[j]}
}

#' @export
fastSBM <- function(z, Pmat){
  csizes = tabulate(z)
  A <- fastSBM.internal(csizes, Pmat)
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig)
  return(A[siginv,siginv])
}


# generate a K-block random graph
fastSBM.internal <- function(csizes, Pmat) {
  # csizes, K by 1, contains the size of each block
  # Pmat is K by K connectivity density matrix
  # store the graph on a sparse matrix
  # generate each block

  sample_block <- function(nr, p, nc = nr, sym=T) {
    # nr is number of rows
    # nc is number of columns
    if (p == 0) return(cbind(NULL,NULL))
    if (sym) {
      m = nr*(nr-1)/2
    } else {
      m = nr*nc
    }

    nnz_count <- rbinom(1, m, p)
    if (nnz_count == 0) return(cbind(NULL,NULL))
    index= sort(MySampleFun(m, nnz_count))

    return(lin2sub(index, nr, sym=sym))
  }

  # K blocks
  subs = NULL
  cumcsizes = cumsum(csizes)
  K = length(csizes)
  if (K == 1) {
    if (length(Pmat) != 1) stop("Pmat should be scalar for ER graphs.")
    subs = sample_block(csizes, Pmat, sym=T)
    n = csizes
  } else {
    if (nrow(Pmat) != K) stop("number of elements of csizes should be the same as nrow of Pmat.")
    n = sum(csizes)
    for (i in 1 : K)
      for (j in i : K) {
        if (i == j) {
          tmp = sample_block(csizes[i], Pmat[i, i], sym=T)
        } else {
          tmp = sample_block(csizes[i], Pmat[i, j], csizes[j], sym=F)
        }
        if (i > 1) {
          tmp[, 1] = tmp[, 1] + cumcsizes[i - 1]
        }
        if (j > 1) {
          tmp[, 2] = tmp[, 2] + cumcsizes[j - 1]
        }
        subs = rbind(subs, tmp)
      }
    # subs = cbind(c(subs[, 1], subs[, 2]),
    #                  c(subs[, 2], subs[, 1]))
    # # indexmat = unique(indexmat)
    # ii = c(subs[, 1], n)
    # jj = c(subs[, 2], n)
    # xx = c(rep(1, nrow(subs)), 0)
    # return(sparseMatrix(ii, jj, x = xx))
  }
  subs = rbind(subs, subs[,2:1])

  ii = c(subs[, 1], n)
  jj = c(subs[, 2], n)
  xx = c(rep(1, nrow(subs)), 0)
  Adj = Matrix::sparseMatrix(ii, jj, x = xx)
  Adj

  # return(sparseMatrix(subs[,1], subs[,2], x = 1, dims=c(n,n))) # this causes speed loss
}

fastDCSBM.internal <- function(csizes, Pmat, theta) {
  # generate a degree-corrected block model
  # P(Aij = 1|ci=k,cj=l,theta) = theta_i * theta_j * P_kl
  # input: csizes - size of each block, Pmat - density across blocks, theta - adjustment of degrees
  # first generate a KblockGraph, and then flip 1s to 0st with probability 1-thetai*thetaj

  A <- fastSBM.internal(csizes, Pmat)
  n <- nrow(A)
  summaryA <- Matrix::summary(A)
  ii <- summaryA$i
  jj <- summaryA$j
  xx <- summaryA$x
  upOnes <- (ii > jj)
  ii <- ii[upOnes]
  jj <- jj[upOnes]
  xx <- xx[upOnes]
  flip1s <- as.numeric(runif(length(ii)) < theta[ii] * theta[jj])
  xx <- xx * flip1s # TODO: xx could still contain 0s
  B <- Matrix::sparseMatrix(c(ii, jj, n),
                    c(jj, ii, n),
                    x = c(xx, xx, 0))
  return(B)
}


## brute force Poisson DCSBM
# consider using c++ to speed up if needed
# Poi_DCSBM <- function(z, B, theta){
#   n = length(z)
#   Theta = diag(theta)
#   Z <- label_vec2mat(z)
#   P <- Theta %*% Z %*% B %*% t(Z) %*% Theta
#   A <- matrix(nrow = n, ncol = n)
#   for (i in 1:n) {
#     for (j in i:n) {
#       A[i,j] <- rpois(1, P[i,j])
#       A[j,i] <- A[i,j]
#     }
#   }
#   return(as(A, "sparseMatrix"))
# }
Poi_DCSBM <- function(z, B, theta){
  n = length(z)
  csizes <- as.vector(table(z))

  sample_block <- function(p, theta1, theta2 = theta1, sym) {
    # nr is number of rows
    # nc is number of columns
    # p is an nr*nc matrix containing parameters
    P <- p*(theta1 %*% t(theta2))
    nr <- length(theta1)
    nc <- length(theta2)
    if (p == 0) return(cbind(NULL,NULL))
    if (sym) {
      A <- matrix(0, nrow = nr, ncol = nr)
      realP <- P[upper.tri(P)]
      sumP <- sum(realP)
      S <- rpois(1, sumP)
      if(S > 0){
        A[upper.tri(A)] <- rmultinom(1, S, realP/sumP)
      }
      return(A)
    } else {
      sumP = sum(P)
      S <- rpois(1, sumP)
      if(S > 0){
        return(matrix(rmultinom(1, S, as.vector(P)/sumP), nrow = nr))
      }else{
        return(matrix(0, nrow = nr, ncol = nc))
      }
    }
  }

  K = length(csizes)
  cumcsizes = cumsum(csizes)
  indices = cbind(c(1,cumcsizes[1:(K-1)]+1), cumcsizes)

  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1 : K){
    for (j in i : K) {
      if (i == j) {
        A[indices[i,1]:indices[i,2], indices[i,1]:indices[i,2]] =
          sample_block(B[i, i], theta[indices[i,1]:indices[i,2]], sym=T)
      } else {
        A[indices[i,1]:indices[i,2], indices[j,1]:indices[j,2]] =
          sample_block(B[i, j], theta[indices[i,1]:indices[i,2]], theta[indices[j,1]:indices[j,2]], sym=F)
      }
    }
  }
  A <- A+t(A)
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig)
  A <- A[siginv,siginv]
  return(as(A, "sparseMatrix"))
}

Poi_SBM <- function(z, B){
  n = length(z)
  csizes <- as.vector(table(z))

  sample_block <- function(p, nr, nc = nr, sym) {
    # nr is number of rows
    # nc is number of columns
    # p is connectivity probability

    if (p == 0) return(cbind(NULL,NULL))
    if (sym) {
      m = nr*(nr-1)/2
      A <- matrix(0, nrow = nr, ncol = nr)
      S <- rpois(1, m*p)
      if(S > 0){
        A[upper.tri(A)] <- rmultinom(1, S, rep(p, m))
      }
      return(A)
    } else {
      m = nr *nc
      S <- rpois(1, m*p)
      if(S > 0){
        return(matrix(rmultinom(1, S, rep(p,m)), nrow = nr))
      }else{
        return(matrix(0, nrow = nr, ncol = nc))
      }
    }
  }

  K = length(csizes)
  cumcsizes = cumsum(csizes)
  indices = cbind(c(1,cumcsizes[1:(K-1)]+1), cumcsizes)

  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1 : K){
    for (j in i : K) {
      if (i == j) {
        A[indices[i,1]:indices[i,2], indices[i,1]:indices[i,2]] =
          sample_block(B[i, i], nr = csizes[i], sym=T)
      } else {
        A[indices[i,1]:indices[i,2], indices[j,1]:indices[j,2]] =
          sample_block(B[i, j], nr = csizes[i], nc = csizes[j], sym=F)
      }
    }
  }
  A <- A+t(A)
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig)
  A <- A[siginv,siginv]
  return(as(A, "sparseMatrix"))
}

