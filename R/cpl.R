#' CPL algorithm for community detection (fast)
#'
#' The Conditional Pseudo-Likelihood (CPL) algorithm for fitting
#' degree-corrected block models
#'
#' The function implements the CPL algorithm as described in the paper below. It
#' relies on the `mixtools` package for fitting a mixture of multinomials to a
#' block compression of the adjacency matrix based on the estimated labels and
#' then reiterates.
#'
#' Technically, `fast_cpl` fits a stochastic block model (SBM) conditional on
#' the observed node degrees,to account for the degree heterogeneity within
#' communities that is not modeled well in SBM. CPL can also be used to
#' effectively estimate the parameters of the degree-corrected block model
#' (DCSBM).
#'
#' The code is an adaptation of the original R code by Aiyou Chen with slight
#' simplifications.
#'
#' @section References: For more details, see [Pseudo-likelihood methods for
#'   community detection in large sparse
#'   networks](https://projecteuclid.org/euclid.aos/1382547514), A. A. Amini, A.
#'   Chen, P. J. Bickel, E. Levina, Annals of Statistics 2013, Vol. 41 (4),
#'   2097—2122.
#'
#'
#' @param Amat adjacency matrix of the network
#' @param K desired number of communities
#' @param ilabels initial label vector (if not provided, initial labels are
#'   estimated using [spec_clust])
#' @param niter number of iterations
#' @return Estimated community label vector.
#' @examples
#' head(fast_cpl(igraph::as_adj(polblogs), 2), 50)
#' @keywords comm_detect
#'
#' @export
fast_cpl <- function(Amat, K, ilabels = NULL, niter = 10) {
  check_pkg_and_stop("mixtools")
  # initial labeling: random compression + spectral clustering
  if (K == 1) return(rep(1, nrow(A)))

  if (is.null(ilabels)) {
    ilabels = spec_clust(Amat, K = K)
  }

  # define temporary variables
  elabels = ilabels
  tlabels = elabels
  # Ohat = as.matrix(BmatFun(Amat, elabels))
  Ohat = compute_block_sums(Amat, elabels)
  # print(Ohat)
  theta = Ohat/apply(Ohat, 1, sum)

  for (niter in 1:niter) {
    # colBmat = as.matrix(colBmatFun(Amat, elabels))
    colBmat = colBmatFun(Amat, elabels)
    tb = as.numeric(table(elabels))
    if (niter > 1) {
      theta = fit$theta
    }
    fit = mixtools::multmixEM(colBmat, lambda = tb/sum(tb), theta = theta, maxit = 300, epsilon = 1e-4)
    maxp = do.call(pmax, data.frame(fit$posterior))
    for (i in 1 : K) {
      tlabels[fit$posterior[, i] == maxp] = i
    }
    # print mis-matching matrix
    # cat('niteration', niter, '\n')
    if (length(unique(tlabels)) == K) {
      elabels = tlabels
      # print(table(clabels, elabels))
      # print(fit$lambda)
    } else {
      break
    }
  }

  # # estimating parameters for the block model
  # L = sum(colBmat)
  # # colBmat = as.matrix(colBmatFun(Amat, elabels))
  # colBmat = colBmatFun(Amat, elabels)
  # fit = mixtools::multmixEM(colBmat, lambda = fit$lambda, theta = fit$theta, maxit = 300, epsilon = 1e-4)
  # # Ohat =  as.matrix(BmatFun(Amat, elabels))
  # Ohat =  computeBlockSums(Amat, elabels)
  #
  # pihat = fit$lambda
  # invTheta = solve(fit$theta)
  # tauhat = c(diag(1/pihat) %*% t(invTheta) %*% rowSums(Ohat))/L
  # Shat = diag(tauhat) %*% fit$theta %*% solve(Ohat/L) %*% t(fit$theta) %*% diag(tauhat)
  # Shat[Shat<0] = 0

  return(elabels)
  # return(list(labels = elabels, pi = pihat, S = Shat))
}

# compress columns only
colBmatFun <- function(Amat, ilabels, transpose = FALSE) {
  # Amat must be a sparse matrix
  sAmat = data.frame(Matrix::summary(Amat))
  if (transpose) {
    ii = sAmat$j
    jj = ilabels[sAmat$i]
  } else {
    ii = sAmat$i
    jj = ilabels[sAmat$j]
  }
  colBmat = Matrix::sparseMatrix(ii, jj, x = sAmat$x)

  # return(colBmat)
  as.matrix(colBmat)
}
