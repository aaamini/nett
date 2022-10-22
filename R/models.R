# Functions to sample from various random network models

#' Calculate the expected average degree of a DCSBM
#'
#' @param n number of nodes
#' @param pri distribution of node labels (K x 1)
#' @param B  connectivity matrix (K x K)
#' @param ex_theta expected value of theta
#' @return expected average degree of a DCSBM
#' @keywords models
#' @export
get_dcsbm_exav_deg <- function(n, pri, B, ex_theta = 1) {
  # Calculate the expected average degree of a DCSBM
  # pri: class prior
  # ex_theta: expected theta
  pri = as.vector(pri)
  as.numeric( (n-1) * t(pri) %*% B %*% pri * (ex_theta^2) )
}

#' Generate planted partition (PP) connectivity matrix
#'
#' Create a degree-corrected planted partition connectivity matrix with a given
#' average expected degree.
#'
#' @param n the number of nodes
#' @param oir out-in-ratio
#' @param lambda the expected average degree
#' @param pri the prior on community labels
#' @param theta node connection propensity parameter of DCSBM
#' @param normalize_theta whether to normalize theta so that max(theta) == 1
#' @param d diagonal of the connectivity matrix. An all-one vector by default.
#' @return The connectivity matrix B of the desired DCSBM.
#' @keywords models
#' @export
pp_conn <- function(n, oir, lambda, pri, theta = rep(1,n), normalize_theta = FALSE, d = rep(1, length(pri))) {
  K = length(pri)
  if (sum(pri) != 1)  pri = pri/sum(pri)
  #if (length(theta) == 1) theta = rep(theta,n)
  if (K > 1) {
    B0 = oir + diag(d-oir)  # diag(scalar) returns an empty matrix if scalar is < 1
  } else {
    B0 = as.matrix(1)
  }

  if (normalize_theta) theta = theta / max(theta)
  # scale = get_sbm_exav_deg(n, theta %*% label_vec2mat(z) / n, B0)
  scale = get_dcsbm_exav_deg(n, pri, B0, mean(theta))
  B = B0*lambda/scale
  if (max(B) > 1) {
    warning("Probabilities truncated at 1.")
    B = pmin(B,1)
  }
  list(B=B, theta=theta)
}

#' Generate randomly permuted connectivity matrix
#'
#' Creates a randomly permuted DCPP connectivity matrix with a given average
#' expected degree
#'
#' The connectivity matrix is a convex combination of a random symmetric
#' permutation matrix and the matrix of all ones, with weights gamm and 1-gamma.
#'
#' @param n number of nodes
#' @param K number of communities
#' @param lambda expected average degree
#' @param gamma a measure of out-in-ratio (convex combination parameter)
#' @param pri the prior on community labels
#' @param theta node connection propensity parameter of DCSBM, by default
#'   E(theta) = 1
#' @return connectivity matrix B of the desired DCSBM.
#' @keywords models
#' @export
gen_rand_conn = function(n, K, lambda, gamma = 0.3, pri = rep(1,K)/K, theta = rep(1,n)) {

  if (sum(pri) != 1)  pri = pri/sum(pri)

  B = matrix(runif(K^2),K)
  B = (B + t(B))/2
  # main structure
  rand_perm_mat = rsymperm(K) # diag(K)[, sample(K)]
  B = (1-gamma)*rand_perm_mat + gamma*B
  scale = get_dcsbm_exav_deg(n, pri, B, mean(theta))
  B*lambda/scale
}


#' Sample from a DCPP
#'
#' Sample from a degree-corrected planted partition model
#'
#' @param n number of nodes
#' @param lambda average degree
#' @param K number of communities
#' @param oir out-in ratio
#' @param theta propensity parameter, if not given will be samples from a Pareto
#' distribution with scale parameter 2/3 and shape parameter 3
#' @param pri prior distribution of node labels
#' @param normalize_theta whether to normalize theta so that max(theta) == 1
#' @return an adjacency matrix following a degree-corrected planted parition model
#' @seealso [sample_dcsbm], [sample_tdcsbm]
#' @export
sample_dcpp <- function(n, lambda, K, oir, theta = NULL,
                       pri = rep(1,K)/K, normalize_theta = FALSE) {
  if (is.null(theta)) theta = EnvStats::rpareto(n, 2/3, 3)
  out = pp_conn(n, oir, lambda, pri, theta, normalize_theta = normalize_theta)
  z = sample(K, n, replace=TRUE, prob=pri)

  list(adj = sample_dcsbm(z, out$B, theta = out$theta), labels = z, B = out$B, theta=out$theta)
}

gkern = function(x,y) exp(-sum((x-y)^2))

#' Sample from a DCLVM
#'
#' A DCLVM with \eqn{K} clusters has edges generated as
#' \deqn{
#'  	 E[\,A_{ij} \mid x, \theta\,] \;\propto\; \theta_i \theta_j e^{- \|x_i - x_j\|^2}
#' }
#' where \eqn{x_i = 2 e_{z_i} + w_i}, \eqn{e_k} is the \eqn{k}th basis vector of \eqn{R^d}, \eqn{w_i \sim N(0, I_d)},
#' and \eqn{\{z_i\} \subset [K]^n}. The proportionality constant is chosen such
#' that the overall network has expected average degree \eqn{\lambda}.
#' To calculate the scaling constant, we approximate \eqn{E[e^{- \|x_i - x_j\|^2}]}
#' for \eqn{i \neq j} by generating random `npairs` \eqn{\{z_i, z_j\}} and average over them.
#'
#' Sample form a degree-corrected latent variable model with Gaussian kernel
#' @param z a vector of cluster labels
#' @param lambda desired average degree of the network
#' @param theta degree parameter
#' @param npairs number of pairs of \eqn{\{z_i, z_j\}}
#' @return Adjacency matrix of DCLVM
#' @export
sample_dclvm = function(z, lambda, theta, npairs = NULL) {
  n = nrow(z)
  if (is.null(npairs)) npairs = n*log10(n)/2

  ex_theta = mean(theta)

  # Approximate computation of E[K(z_i,z_j)] for i \neq j assuming {z_i} are i.i.d.
  # Generate random npairs {z_i, z_j} and average over them
  pidx = unique(matrix(sample(n, 2*npairs, TRUE), ncol=2))
  # pidx[apply(pidx, 1, function(x) x[1] != x[2]),]
  ex_kern = mean( sapply(1:nrow(pidx), function(i) gkern(z[pidx[i,1],], z[pidx[i,2],])) )

  scale = lambda/((n-1)*(ex_theta)^2*ex_kern)

  return( sample_dclvm_cpp(z, scale, theta) )
}

#' Sample from a DCSBM
#'
#' Sample an adjacency matrix from a degree-corrected block model (DCSBM)
#'
#' @param z Node labels (\eqn{n * 1})
#' @param B Connectivity matrix (\eqn{K * K})
#' @param theta Node connectivity propensity vector (\eqn{n * 1})
#' @return An adjacency matrix following DCSBM
#' @seealso [sample_dcpp], [fast_sbm], [sample_tdcsbm]
#' @examples
#' B = pp_conn(n = 10^3, oir = 0.1, lambda = 7, pri = rep(1,3))$B
#' head(sample_dcsbm(sample(1:3, 10^3, replace = TRUE), B, theta = rexp(10^3)))
#' @export
sample_dcsbm = function(z, B, theta=1) {
  n = length(z)
  if (length(theta) == 1) theta = rep(theta,n)
  sample_dcsbm_cpp(z-1, B, theta)
}

#' Sample from a DCER
#'
#' Sample an adjacency matrix from a degree-corrected Erdős–Rényi model (DCER).
#'
#' @param theta Node connectivity propensity vector (\eqn{n * 1})
#' @return An adjacency matrix following DCSBM
#' @export
sample_dcer = function(theta) {
  sample_dcer_cpp(theta)
}

#' Sample from a SBM (fast)
#'
#' Samples an adjacency matrix from a stochastic block model (SBM)
#'
#' The function implements a fast algorithm for sampling sparse SBMs, by only
#' sampling the necessary nonzero entries. This function is adapted almost
#' verbatim from the original code by Aiyou Chen.
#'
#' @param z Node labels (\eqn{n * 1})
#' @param B Connectivity matrix (\eqn{K * K})
#' @return An adjacency matrix following SBM
#'
#' @examples
#' B = pp_conn(n = 10^4, oir = 0.1, lambda = 7, pri = rep(1,3))$B
#' head(fast_sbm(sample(1:3, 10^4, replace = TRUE), B))
#'
#' @keywords models
#' @export
fast_sbm <- function(z, B){
  if (!isSymmetric(B)) stop("B has to be a symmetric matrix.")
  csizes = tabulate(z, nrow(B))
  A <- fast_sbm.internal(csizes, B)
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig)
  return(A[siginv,siginv])
}

# generate a K-block random graph
fast_sbm.internal <- function(csizes, Pmat) {
  # csizes, K by 1, contains the size of each block
  # Pmat is K by K connectivity density matrix
  # store the graph on a sparse matrix
  # generate each block

  sample_block <- function(nr, p, nc = nr, sym=TRUE) {
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
    subs = sample_block(csizes, Pmat, sym=TRUE)
    n = csizes
  } else {
    if (nrow(Pmat) != K) stop("number of elements of csizes should be the same as nrow of Pmat.")
    n = sum(csizes)
    for (i in 1 : K)
      for (j in i : K) {
        if (i == j) {
          tmp = sample_block(csizes[i], Pmat[i, i], sym=TRUE)
        } else {
          tmp = sample_block(csizes[i], Pmat[i, j], csizes[j], sym=FALSE)
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

# Previous model simulation code ---------------------------------------------------------


quickDCSBM <- function(n, lambda, K, oir, theta = NULL,
                       pri = rep(1,K)/K, normalize_theta = FALSE) {
  # pri <- (1:K)/K
  # B0 = oir + diag(rep(1 - oir, K))
  # # B <- pp_conn(n, oir, lambda, pr)
  if (is.null(theta)) theta = EnvStats::rpareto(n, 2/3, 3)
  out = pp_conn(n, oir, lambda, pri, theta, normalize_theta = normalize_theta)
  # if (normalize_theta) theta = theta / max(theta)
  # # scale = get_sbm_exav_deg(n, theta %*% label_vec2mat(z) / n, B0)
  # # scale = get_dcsbm_exav_deg(n, pri, B0, mean(theta))
  # scale = get_dcsbm_exav_deg(n, pri, B0, mean(theta))
  # B = B0*lambda/scale
  # if (max(B) > 1) {
  #   warning("Probabilities truncated at 1.")
  #   B = pmin(B,1)
  # }
  z <- sample(K, n, replace=TRUE, prob=pri)

  list(adj=sample_tdcsbm(z, out$B, theta = out$theta), labels = z, B = out$B, theta=out$theta)
}

#' Sample truncated DCSBM (fast)
#'
#' Sample an adjacency matrix from a truncated degree-corrected block model
#' (DCSBM) using a fast algorithm.
#'
#' The function  samples an adjacency matrix from a truncated DCSBM, with
#' entries having Bernoulli distributions with mean \deqn{ E[A_{ij} | z] =
#' B_{z_i, z_j} \min(1, \theta_i \theta_j).} The approach uses the masking idea
#' of Aiyou Chen, leading to fast sampling for sparse networks. The masking,
#' however, truncates \eqn{\theta_i \theta_j} to at most 1, hence
#' we refer to it as the truncated DCSBM.
#'
#' @param z Node labels (\eqn{n * 1})
#' @param B Connectivity matrix (\eqn{K * K})
#' @param theta Node connectivity propensity vector (\eqn{n * 1})
#' @return An adjacency matrix following DCSBM
#' @examples
#' B = pp_conn(n = 10^4, oir = 0.1, lambda = 7, pri = rep(1,3))$B
#' head(sample_tdcsbm(sample(1:3, 10^4, replace = TRUE), B, theta = rexp(10^4)))
#' @export
sample_tdcsbm <- function(z, B, theta=1) {
  csizes = tabulate(z)
  n = length(z)
  if (length(theta) == 1) theta = rep(theta,n)
  A <- sample_tdcsbm.internal(csizes, B, theta)

  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig) # or order(sig)
  return(A[siginv,siginv]) # E[A_{siginv[i], siginv[j]}] = B_{z[i], z[j]}
}


sample_tdcsbm.internal <- function(csizes, Pmat, theta) {
  # generate a degree-corrected block model
  # P(Aij = 1|ci=k,cj=l,theta) = theta_i * theta_j * P_kl
  # input: csizes - size of each block, Pmat - density across blocks, theta - adjustment of degrees
  # first generate a KblockGraph, and then flip 1s to 0st with probability 1-thetai*thetaj

  A <- fast_sbm.internal(csizes, Pmat)
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

#' @import stats
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
          sample_block(B[i, i], theta[indices[i,1]:indices[i,2]], sym=TRUE)
      } else {
        A[indices[i,1]:indices[i,2], indices[j,1]:indices[j,2]] =
          sample_block(B[i, j], theta[indices[i,1]:indices[i,2]], theta[indices[j,1]:indices[j,2]], sym=FALSE)
      }
    }
  }
  A <- A+t(A)
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig)
  A <- A[siginv,siginv]
  return(as(A, "sparseMatrix"))
}

#' @import stats
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
          sample_block(B[i, i], nr = csizes[i], sym=TRUE)
      } else {
        A[indices[i,1]:indices[i,2], indices[j,1]:indices[j,2]] =
          sample_block(B[i, j], nr = csizes[i], nc = csizes[j], sym=FALSE)
      }
    }
  }
  A <- A+t(A)
  sig = order(z) # z[sig] = 1 1 1 ... 2 2 2 ... 3 3 3 3 ...
  siginv = Matrix::invPerm(sig)
  A <- A[siginv,siginv]
  return(as(A, "sparseMatrix"))
}

