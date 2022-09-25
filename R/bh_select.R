#' Beth-Hessian model selection
#'
#' Estimate the number of communities under block models by using the spectral
#' properties of network Beth-Hessian matrix with moment correction.
#'
#' @section References:
#' [Estimating the number of communities in networks by spectral methods](https://arxiv.org/abs/1507.00827),
#' Can M. Le, Elizaveta Levina, arXiv preprint arXiv:1507.00827, 2015
#'
#' @param A adjacency matrix.
#' @param Kmax the maximum number of communities to check.
#' @return A list of result
#' \item{K}{estimated the number of communities}
#' \item{rho}{eigenvalues of the Beth-Hessian matrix}
#' @importFrom methods as
#' @keywords mod_sel
#'
#' @export
bethe_hessian_select <- function(A, Kmax){
    if (Kmax <= 2)
      Kmax <- 2

    d <- Matrix::colSums(A)
    n <- nrow(A)
    I <- as(diag(rep(1, n)), "dgCMatrix")
    D <- as(diag(d), "dgCMatrix")
    r <- sqrt(mean(d))
    BH <- (r^2 - 1) * I - r * A + D
    rho <- sort(RSpectra::eigs_sym(BH, Kmax, which = "SA")$values)
    diff <- rho[2:Kmax] - 5 * rho[1:(Kmax - 1)]
    return(list(K = max(which(diff > 0)), values = rho))
}
