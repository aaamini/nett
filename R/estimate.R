
# Parameter estimation ----------------------------------------------------

# The idea behind the following function is due to Aiyou Chen
# require A to be a sparse matrix
# This is faster than looping over [K] x [K].

#' Block sum of an adjacency matrix
#'
#' Compute the block sum of an adjacency matrix given a label vector.
#' @param A adjacency matrix.
#' @param z label vector.
#' @return A K x L matrix with (k,l)-th element as \eqn{sum_{i,j} A_{i,j} 1{z_i = k, z_j = l}}
#' @keywords estimation
#' @export
compute_block_sums <- function(A, z) {
  # A: a sparse adjacency matrix
  # z: a label vector
  # Outputs the matrix B[k,l] = sum_{i,j} A[i,j] 1{z_i = k, z_j = l}
  sA = Matrix::summary(A)
  return(as.matrix(Matrix::sparseMatrix(i = z[sA$i], j = z[sA$j], x = sA$x)))
}


# estimate parameters of the block model -- can be removed since estim_dcsbm is enough
# in fact estim_sbm = function(A,z) estim_dcsbm(A,z)$B
estim_sbm <- function(A,z) estim_dcsbm(A,z)$B
# estim_sbm <- function(A, z) {
#   ns <- as.vector(table(z))
#   if (length(ns) == 1) { # number of clusters == 1
#     return(as.matrix(max(sum(A),1)/(ns*(ns-1))))
#   }
#   # pmax(computeBlockSums(A, z),1) / (outer(ns,ns) - diag(ns))
#   # computeBlockSums(A, z) / (outer(ns,ns) - diag(ns))
#   Bsum = computeBlockSums(A, z)
#   Bsum[Bsum == 0] = 1 # avoid estimating 0
#   Bsum / (ns %*% t(ns) - diag(ns))
# }

# estimate parameters of the degree-corrected block model
# fast implementation

#' Estimate model parameters of a DCSBM
#'
#' Compute the block sum of an adjacency matrix given a label vector.
#' @param A adjacency matrix.
#' @param z label vector.
#' @return A list of result
#' \item{B}{estimated connectivity matrix.}
#' \item{theta}{estimated node propensity parameter.}
#' @details
#' \deqn{\hat B_{k\ell} = \frac{N_{k\ell}(\hat z)}{m_{k\ell} (\hat z)}, \quad \hat \theta_i =  \frac{n_{\hat z_i}(\hat z) d_i}{\sum_{j : \hat z_j = \hat z_i} d_i}}
#' where \eqn{N_{k\ell}(\hat{z})} is the sum of the elements of \code{A} in block \eqn{(k,\ell)}
#' specified by labels \eqn{\hat z}, \eqn{n_k(\hat z)} is the number of nodes in community \eqn{k}
#' according to \eqn{\hat z} and \eqn{m_{k\ell}(\hat z) = n_k(\hat z) (n_\ell(\hat z) - 1\{k = \ell\})}
#' @keywords estimation
#' @export
estim_dcsbm <- function(A,z) {
  ns <- as.vector(table(z)) # nk = tabulate(z)
  degs = Matrix::rowSums(A)

  if (length(ns) == 1) { # number of clusters == 1
    Bsum = sum(A)
    B = as.matrix(max(Bsum,1)/max(ns*(ns-1)),1)
    theta = degs*ns/Bsum

  } else {  # number of clusters > 1
    Bsum = compute_block_sums(A, z)
    Bsum[Bsum == 0] = 1
    B = Bsum / pmax(ns %*% t(ns) - diag(ns), 1)

    total_clust_degs = rowSums(Bsum)  # sum of degrees within each cluster
    theta = degs*ns[z]/total_clust_degs[z]
  }

  return(list(B=B, theta=theta))
}

# TODO: check dimensions
# TODO: extend to non-consequential labels, and automatically detect max. label
# Likelihood computations -------------------------------------------------

#' Log likelihood of a DCSBM (fast with poi = TRUE)
#'
#' Compute the log likelihood of a DCSBM, using estimated parameters
#' B, theta based on the given label vector
#'
#' The log likelihood is calculated by
#' \deqn{\ell(\hat B,\hat \theta, \hat \pi, \hat z \mid A) =
#' \sum_i \log \hat \pi_{z_i} + \sum_{i < j} \phi(A_{ij};\hat \theta_i \hat \theta_j \hat B_{\hat{z}_i \hat{z}_j} )}
#' where \eqn{\hat B}, \eqn{\hat \theta} is calculated by [estim_dcsbm],
#' \eqn{\hat{\pi}_k} is the proportion of nodes in community k.
#'
#' @param A adjacency matrix
#' @param z label vector
#' @param poi whether to use Poisson version of likelihood
#' @param eps truncation threshold for the Bernoulli likelihood,
#' used when parameter phat is close to 1 or 0.
#' @return log likelihood of a DCSBM
#' @seealso [eval_dcsbm_loglr], [eval_dcsbm_bic]
#' @keywords estimation
#' @export
eval_dcsbm_like <- function(A, z, poi = TRUE, eps = 1e-6) {
  Bsum = compute_block_sums(A,z)
  degs = Matrix::rowSums(A)
  total_clust_degs = Matrix::rowSums(Bsum)
  theta = degs/total_clust_degs[z] # unit normalization for theta

  Z = label_vec2mat(z)
  Zth = Z * as.numeric(theta)

  sA = Matrix::summary(A)
  ii = sA$i
  jj = sA$j
  xx = sA$x # xx could still contain 0s
  nz_idx = (xx > 0) & (jj > ii) # pick the upper triangular part
  ix = ii[nz_idx]
  jx = jj[nz_idx]
  aa = xx[nz_idx]

  zi = z[ix]
  zj = z[jx]
  pp =  truncate_to_ab(theta[ix] * theta[jx] * Bsum[(zi-1)*ncol(Bsum)+zj], eps, 1-eps)

  if (!poi) { # Bernoulli, slow computation, high mem
    warning('Bernoulli likelihood computation is slow. Try "poi = TRUE" option.')
    term1 = sum(aa * log(pp/(1-pp)))
    mm = truncate_to_ab(Zth %*% Bsum %*% t(Zth), eps, 1-eps) # mean matrix
    term2 = sum( log(1-mm[which(upper.tri(mm))]) )

  } else { # Poisson
    term1 = sum(aa * log(pp))
    # The next line computes \sum_{i < j} Phat_{ij} where Phat_{ij} = \theta_i \theta_j Bsum_{z_i, z_j}
    # assuming unit theta-normalization
    # that is, \theta \sum_{z_i = k} \theta_i = 1
    #term2 = sum( (1-colSums(Zth^2))*diag(Bsum)/2 ) + (sum(Bsum) - sum(diag(Bsum)))/2  # fast
    # so that n is hided
    term2 = sum(Bsum)/2 - sum(colSums(Zth^2)*diag(Bsum)/2)
    term2 = -term2

    # # slow, high mem approach
    # mm = truncate_to_ab(Zth %*% Bsum %*% t(Zth), eps, 1-eps) # mean matrix
    # term2 = - sum( mm[which(upper.tri(mm))] )
  }

  n = nrow(A)
  ns = as.vector(table(z))
  ns[ns == 0] = 1
  term3 = sum(ns * log(ns/n))
  return(term1 + term2 + term3)
}


#' Log-likelihood ratio of two DCSBMs (fast with poi = TRUE)
#'
#' Computes the log-likelihood ratio of one DCSBM relative to another, using
#' estimated parameters `B` and `theta` based on the given label vectors.
#'
#' The log-likehood ratio is computed between two DCSBMs specified by the columns
#' of `labels`. The function computes the log-likelihood ratio of the model with
#' `labels[ , 2]` w.r.t. the model with `labels[ , 1]`. This is often used with two
#' label vectors fitted using different number of communities (say `K` and `K+1`).
#'
#' When `poi` is set to `TRUE`, the function uses fast sparse matrix computations
#' and is scalable to large sparse networks.
#'
#' @param A adjacency matrix
#' @param labels a matrix with two columns representing two different label
#'   vectors
#' @param poi whether to use Poisson version of likelihood (instead of Bernoulli)
#' @param eps truncation threshold for the Bernoulli likelihood, used when
#'   parameter phat is close to 1 or 0.
#' @return log-likelihood ratio
#' @seealso [eval_dcsbm_like], [eval_dcsbm_bic]
#' @keywords mod_sel
#' @export
eval_dcsbm_loglr = function(A, labels, poi = TRUE, eps = 1e-6) {
  eval_dcsbm_like(A, labels[ , 2], poi = poi, eps = eps) - eval_dcsbm_like(A, labels[ , 1], poi = poi, eps = eps)
}

#' Compute BIC score
#' @description compute BIC score when fitting a DCSBM to network data
#' @param A adjacency matrix
#' @param z label vector
#' @param K number of community in \code{z}
#' @param poi whether to use Poisson version of likelihood
#' @return BIC score
#' @details the BIC score is calculated by -2*log likelihood minus \eqn{K\times(K + 1)\times log(n)}
#' @section References:
#' BIC score is originally proposed in
#' [Likelihood-based model selection for stochastic block models](https://projecteuclid.org/euclid.aos/1494921948)
#' Wang, YX Rachel, Peter J. Bickel, The Annals of Statistics 45, no. 2 (2017): 500-528.
#'
#' The details of modified implementation can be found in
#' [Adjusted chi-square test for degree-corrected block models](https://arxiv.org/abs/2012.15047),
#' Linfan Zhang, Arash A. Amini, arXiv preprint arXiv:2012.15047, 2020.
#'
#' @seealso [eval_dcsbm_like], [eval_dcsbm_loglr]
#' @keywords mod_sel
#' @export
eval_dcsbm_bic = function(A, z, K, poi) {
  n = length(labels)
  # eval_dcsbm_like(A, z = z, poi = poi) - K*(K + 1)*log(n)/2
  -2*eval_dcsbm_like(A, z = z, poi = poi) + K*(K + 1)*log(n)
}

# Old functions -----------------------------------------------------------
# This is only for testing
computeBlockSums2 <- function(A, z) {
  K <- length(unique(z))
  Bsum <- matrix(0, K, K)
  idx <- sapply(1:K, function(k) z==k) # potentially probelmatic if K >= n

  for (i in 1:K) {
    Bsum[i,i] <- sum(A[idx[,i],idx[,i]])
  }
  if (K > 1){
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        Bsum[j,i] <- Bsum[i,j] <- sum(A[idx[,i],idx[,j]])
      }
    }
  }
  return(Bsum)
}


estimSBM <- function(A, z) {
  K <- length(unique(z))
  B <- matrix(0, K, K)
  idx <- sapply(1:K, function(k) z==k) # potentially probelmatic if K >= n
  # ns <- colSums(idx)
  ns <- as.vector(table(z))
  for (i in 1:K) {
    B[i,i] <- max( sum(A[idx[,i],idx[,i]]), 1) / (ns[i]*(ns[i]-1))
  }
  if (K > 1){
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        B[j,i] <- B[i,j] <- max( sum(A[idx[,i],idx[,j]]), 1) / (ns[i]*ns[j])
      }
    }
  }
  B[B>1] <- 1
  return(B)
}

# estimate parameters of the degree-corrected block model
#' @import stats
estimDCSBM <- function(A, z) {
  B = estimSBM(A, z)
  nk = tabulate(z)
  # tmp = B %*% nk
  degs = Matrix::rowSums(A)
  # theta = degs / tmp[z]
  #
  total_clust_degs = aggregate(degs, list(label=z), sum)$x
  # total_clust_degs = aggregate(degs, list(label=z), mean)$x
  # theta2 = degs / total_clust_degs[z]

  theta = degs*nk[z]/total_clust_degs[z]

  # total_clust_degs = aggregate(degs, list(label=z), mean)$x
  # # K = unique(z)
  # # total_clust_degs = sapply(1:K, function(k) mean(degs[z==k]))
  # theta2 = degs / total_clust_degs[z]


  return(list(B=B, theta=theta))
}
