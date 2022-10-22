
l2_norm = function(x) sqrt(sum(x^2))

#' Sinkhorn--Knopp matrix scaling
#'
#' Implements the Sinkhorn--Knopp algorithm for transforming a square matrix
#' with positive entries to a stochastic matrix with given common row and column
#' sums (e.g., a doubly stochastic matrix).
#'
#' Computes diagonal matrices D1 and D2 to make D1*A*D2 into a matrix with
#' given row/column sums. For a symmetric matrix `A`, one can set `sym = TRUE` to
#' compute a symmetric scaling D*A*D.
#'
#' @param A input matrix
#' @param sums desired row/column sums
#' @param niter number of iterations
#' @param tol convergence tolerance
#' @param sym whether to compute symmetric scaling D A D
#' @param verb whether to print the current change
#' @return Diagonal matrices D1 and D2 to make D1*A*D2 into a matrix with
#' given row/column sums.
#' @export
#' @keywords utils
sinkhorn_knopp = function(A, sums = rep(1, nrow(A)),
                          niter = 100, tol = 1e-8, sym = FALSE, verb = FALSE) {


  delta = Inf
  r = c = rep(1, nrow(A))
  converged = FALSE
  t = 1
  while( t <= niter && !converged) {
    r = sums / (A %*% c)
    cnew = sums / (t(A) %*% r)

    # Symmetric Sinkhorn--Knopp algorithm could oscillate between two sequences,
    # need to bring the two sequences together (See for example "Algorithms For
    # The Equilibration Of Matrices And Their Application To Limited-memory
    # Quasi-newton Methods")
    if (sym) cnew = (cnew + r)/2

    delta = l2_norm(cnew-c)
    if (delta < tol) converged = TRUE
    if (verb) nett::printf("err = %3.5e\n", delta)
    c = cnew
    t = t+1
  }
  list(r = as.numeric(r), c = as.numeric(c))
}


# row_normalize = function(mat) {
#   scales = 1/rowSums(mat)
#   list(mat = mat * scales[row(mat)], d = scales)
# }

# col_normalize = function(mat) {
#   scales = 1/colSums(mat)
#   list(mat = mat * scales[col(mat)], d = scales)
# }


# # test
set.seed(123)
n = 5
B = matrix(runif(n^2), nrow=n)
A = (B + t(B))/2
h = c(2,5,4,1,1.4)
out = sinkhorn_knopp( A , sums=h, sym=TRUE)
out
r = out$r
c = out$c
rowSums(diag(r) %*% A %*% diag(c))
rowSums(diag(r) %*% A %*% diag(r))

d = r / h
diag(d) %*% A %*% diag(d) %*% h

