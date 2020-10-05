# CPL
#' @export
fastCPL <- function(Amat, Kcommunities, ilabels = NULL, Iter = 10) {
  # initial labeling: random compression + spectral clustering
  if (Kcommunities == 1) return(rep(1,nrow(A)))

  if (is.null(ilabels)) {
    ilabels = spec_clust(Amat, K = Kcommunities)
  }

  # define temporary variables
  elabels = ilabels
  tlabels = elabels
  # Ohat = as.matrix(BmatFun(Amat, elabels))
  Ohat = computeBlockSums(Amat, elabels)
  # print(Ohat)
  theta = Ohat/apply(Ohat, 1, sum)

  for (iter in 1:Iter) {
    # colBmat = as.matrix(colBmatFun(Amat, elabels))
    colBmat = colBmatFun(Amat, elabels)
    tb = as.numeric(table(elabels))
    if (iter > 1) {
      theta = fit$theta
    }
    fit = mixtools::multmixEM(colBmat, lambda = tb/sum(tb), theta = theta, maxit = 300, epsilon = 1e-4)
    maxp = do.call(pmax, data.frame(fit$posterior))
    for (i in 1 : Kcommunities) {
      tlabels[fit$posterior[, i] == maxp] = i
    }
    # print mis-matching matrix
    # cat('iteration', iter, '\n')
    if (length(unique(tlabels)) == Kcommunities) {
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
