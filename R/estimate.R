
# Parameter estimation ----------------------------------------------------

# The idea behind the following function is due to Aiyou Chen
# require A to be a sparse matrix
# This is faster than looping over [K] x [K].
#' @export
computeBlockSums <- function(A, z) {
  # A: a sparse adjacency matrix
  # z: a label vector
  # Outputs the matrix B[k,l] = sum_{i,j} A[i,j] 1{z_i = k, z_j = l}
  sA = Matrix::summary(A)
  return(as.matrix(Matrix::sparseMatrix(i = z[sA$i], j = z[sA$j], x = sA$x)))
}


# estimate parameters of the block model -- can be removed since estim_dcsbm is enough
# in fact estim_sbm = function(A,z) estim_dcsbm(A,z)$B
#' @export
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
#' @export
estim_dcsbm <- function(A,z) {
  ns <- as.vector(table(z)) # nk = tabulate(z)
  degs = Matrix::rowSums(A)

  if (length(ns) == 1) { # number of clusters == 1
    Bsum = sum(A)
    B = as.matrix(max(Bsum,1)/(ns*(ns-1)))
    theta = degs*ns/Bsum

  } else {  # number of clusters > 1
    Bsum = computeBlockSums(A, z)
    Bsum[Bsum == 0] = 1
    B = Bsum / (ns %*% t(ns) - diag(ns))

    total_clust_degs = rowSums(Bsum)  # sum of degrees within each cluster
    theta = degs*ns[z]/total_clust_degs[z]
  }

  return(list(B=B, theta=theta))
}

# TODO: check dimensions
# TODO: extend to non-consequential labels, and automatically detect max. label
# Likelihood computations -------------------------------------------------
#' @export
eval_dcsbm_like <- function(A, z, poi = T, eps = 1e-6) {
  Bsum = computeBlockSums(A,z)
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
    warning('Bernoulli likelihood computation is slow. Try "poi = T" option.')
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

# computes likelihood ratio of labels[ , 2]-model w.r.t. labels[ , 1]-model
#' @export
eval_dcsbm_loglr = function(A, labels, poi = T, eps = 1e-6) {
  eval_dcsbm_like(A, labels[ , 2], poi = poi, eps = eps) - eval_dcsbm_like(A, labels[ , 1], poi = poi, eps = eps)
}

# compute BIC score
#' @export
eval_dcsbm_bic = function(A, z, K, poi = T) {
  n = length(labels)
  eval_dcsbm_like(A, z = z, poi = poi) - K*(K + 1)*log(n)/2
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

#' @export
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
#' @export
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
