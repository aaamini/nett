# label_vec2mat <- function(z, K=NULL){
#   if (is.null(K)) {
#     K <- max(unique(z))
#   }
#   # diag(K)[z,]
#   Diagonal(K)[z,]
# }

#' Convert label matrix to vector
#' @param Z  a cluster assignment matrix
#' @return A label vector that follows the assignment matrix
#' @export
#' @keywords utils
label_mat2vec <- function(Z){
  max.col(Z)
}

#' Compute confusion matrix
#' @param z a label vector
#' @param y a label vector
#' @param K number of labels in both \code{z} and \code{y}
#' @return A `K`x`K` confusion matrix between `z` and `y`
#' @export
#' @keywords evaluation
compute_confusion_matrix <- function (z, y, K=NULL) {
  # Compute the confusion matrix between labels "y" and "z"
  # z,y Two sets of labels
  # K   number of labels in both "c" and "e"

  if (is.null(K)) K = max(c(z,y))
  t(label_vec2mat(z,K)) %*% label_vec2mat(y,K)
}

# M = label_vec2mat(c,K)'*label_vec2mat(e,K);

#' Convert label vector to matrix
#' @param z a label vector
#' @param K number of labels in \code{z}
#' @param sparse whether the output should be sparse matrix
#' @return A cluster assignment matrix that follows from the label vector `z`
#' @export
#' @keywords utils
label_vec2mat <- function(z, K=NULL, sparse=FALSE) {
  if (is.null(K)) K <- max(z)

  if (K==1)
    return( as.matrix( rep(1,length(z)) , ncol=1) )
  else {
    if (sparse) {
      return( Matrix::Diagonal(K)[z,] )
    } else {
      return( diag(K)[z,] )
    }
  }
}

#' Compute normalized mutual information (NMI)
#'
#' Compute the NMI between two label vectors with the same cluster number
#' @param z a label vector
#' @param y a label vector
#' @return NMI between `z` and `y`
#' @export
#' @keywords evaluation
compute_mutual_info  <- function(z,y) {
  # normMUI Computes the normalized mutual information between two clusters
  #  Labels should be either vectors or n x k matrices

  # c = turn_into_column_ifvec(c);
  # e = turn_into_column_ifvec(e);

  # if ( !is.null(dim(z)) ) z = label_mat2vec(z)
  # if ( !is.null(dim(y)) ) z = label_mat2vec(y)

  CM = compute_confusion_matrix(z,y)
  normCM = CM / sum(CM); # normalized confusion matrix
  IDX = CM == 0 # index of zero elements of CM so that we can avoid them

  jointEnt = - sum( (normCM[!IDX])*log(normCM[!IDX]) )
  indpt = matrix(rowSums(normCM),ncol=1) %*% matrix(colSums(normCM),nrow=1)
  muInfo = sum(normCM[!IDX] * log(normCM[!IDX] / indpt[!IDX]) )

  muInfo / jointEnt
}


## Written by Aiyou Chen
lin2sub <- function(index, nr, sym=FALSE) {
  # traslate linear index to row and column indices
  # nr is the number of rows
  # index is the linear index assuming columns are concatenated into a vector
  if (sym) {
    kk = 1 : (nr-1)
    ss = c(.5, kk * (kk + 1)/2 + .5)
    # icut = as.numeric(cut(index, ss))
    # cc = icut + 1
    cc = ceiling(sqrt(1+8*index)/2 + 0.5)
    rr = index - (ss[cc - 1] - .5)
  } else {
    cc = ceiling(index/nr)
    rr = index - (cc-1) * nr
  }
  return(cbind(rr,cc))
}
###
truncate_to_ab = function(x, a, b) {
  pmin(pmax(x, a), b)
}

#' The usual "printf" function
#' @param ... printing object
#' @return the value of the printing object
#' @keywords utils
#' @export
printf <- function(...) cat(sprintf(...))


customSampleFun <- function(m, n) {
  # sample n numbers from 1 : m
  # sample() is too slow when m is large
  while (1) {
    x = ceiling(runif(n*2) * m)
    y = unique(x)
    if (length(y) >= n) {
      break
    }
  }
  return(sample(y, n))
}

# MySampleFun <- customSampleFun
# MySampleFun <- function(m,n) sample.int(m, n, useHash = TRUE) # gives error when n > m/2
MySampleFun <- function(m,n) igraph::sample_seq(1,m,n)



# check MySampleFun
# MySampleFun <- sample

check_pkg_and_stop = function(pkg, func_name = NULL) {
  if (is.null(func_name)) {
    func_name = "this function"
  } else {
    func_name = paste0(func_name, "()")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package \"%s\" needed for %s. Please install it.", pkg, func_name),
      call. = FALSE)
  }
}

#' Generate random symmetric permutation matrix
#'
#' Generate a random symmetric permutation matrix (recursively)
#'
#' @param K size of the matrix
#' @returns A random `K` x `K` symmetric permutation matrix
#' @keywords utils
#' @export
rsymperm = function(K) {
  if (K == 1) return(1)
  B = matrix(0, K, K)
  r = sample(1:K, 1)
  B[1,r] = B[r,1] = 1
  idx = setdiff(1:K, c(1,r))
  nidx = length(idx)
  if (nidx > 0) {
    B[idx, idx] = rsymperm(nidx)
  }
  B
}
