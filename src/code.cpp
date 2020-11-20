#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;

// arma::sp_mat sample_dcsbm(arma::uvec z, arma::mat Pmat, arma::vec theta) {
// [[Rcpp::export]]
arma::sp_mat sample_dcsbm_cpp(arma::uvec z, arma::mat Pmat, arma::vec theta) {
  int n = theta.n_elem;
  int locs_len = 2*n;
  int loc_idx = 0;

  // theta.print();
  arma::umat locs(2, locs_len);

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      // Rcout << locs_len << "\n";
      // Rcout << z(i) << z(j) << "\n";
      // Rcout << loc_idx << " " << locs_len << "\n";
      double p = Pmat(z(i),z(j))*theta(i)*theta(j);
      // Rcout << i+1 << ", " << j+1 << ": " << p;

      if (R::runif(0,1) < p) {
        //Rcout << "+";
        if (loc_idx >= locs_len) {
          locs_len = round(locs_len*1.5);
          // Rcout << locs_len << "\n";
          locs.resize(2, locs_len);
        }
        locs(0, loc_idx) = i;
        locs(1, loc_idx) = j;
        //Rcout << "\n";
        //locs.t().print();
        loc_idx++;
      }
      // Rcout << "\n";
    }
  }
  locs.resize(2, loc_idx);

  locs = arma::join_horiz(locs, arma::join_vert(locs.row(1),locs.row(0)));
  //locs.print();
  arma::sp_mat A(locs, arma::ones<arma::vec>(locs.n_cols), n, n);
  return A;
}

// [[Rcpp::export]]
arma::sp_mat sample_r1(arma::vec theta, double p = 1) {
  // vec tempu(theta.n_elem, fill::randu);

  int n = theta.n_elem;
  arma::umat locs(2,0);
  for (int i = 0; i < n; i++) {
    arma::vec tempu = arma::randu(n);
    arma::uvec idx = find(tempu < p * theta(i) * theta);

    // Rcout << join_horiz(tempu, theta(i) * theta);
    int m = idx.n_elem;
    // Rcout << join_horiz(i*ones<uvec>(m),idx).t();
    locs = arma::join_horiz(locs, arma::join_horiz(i*arma::ones<arma::uvec>(m),idx).t());
  }

  // Rcout << locs;
  arma::sp_mat A(locs, arma::ones<arma::vec>(locs.n_cols));

  return A;
}

/*** R
sample_dcsbm_cpp(c(0,0,0,1,1,1), matrix(c(1,0,0,1),nrow=2), theta = c(1,1,1,1,1,1))
# sample_r1(c(0.1,0.5,0.9,0.9,0.7))
*/
