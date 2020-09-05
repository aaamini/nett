# label_vec2mat <- function(z, K=NULL){
#   if (is.null(K)) {
#     K <- max(unique(z))
#   }
#   # diag(K)[z,]
#   Diagonal(K)[z,]
# }

#' @export
label_mat2vec <- function(Z){
  max.col(Z)
}

#' @export
compute_confusion_matrix <- function (z, y, K=NULL) {
  # Compute the confusion matrix between labels "y" and "z"
  # z,y Two sets of labels
  # K   number of labels in both "c" and "e"

  if (is.null(K)) K = max(c(z,y))
  t(label_vec2mat(z,K)) %*% label_vec2mat(y,K)
}

# M = label_vec2mat(c,K)'*label_vec2mat(e,K);

#' @export
label_vec2mat <- function(z, K=NULL, sparse=F) {
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

#' @export
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
#' @export
lin2sub <- function(index, nr, sym=F) {
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

truncate_to_ab = function(x, a, b) {
  pmin(pmax(x, a), b)
}

