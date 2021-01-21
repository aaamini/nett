#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma;

// [[Rcpp::export]]
double pair_dist2(arma::vec x, arma::vec y) {
  return( arma::sum(arma::square(x-y)) );
}

// [[Rcpp::export]]
arma::mat pair_dist2_mat(arma::mat z) {
  int n = z.n_rows;
  arma::mat D(n, n, arma::fill::zeros);

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      D(i,j) = pair_dist2(z.row(i).t(), z.row(j).t());
      D(j,i) = D(i,j);
    }
  }

  return(D);
}


//
// [[Rcpp::export]]
arma::sp_mat sample_dclvm_cpp(arma::mat z, double scale, arma::vec theta) {
  int n = z.n_rows;
  int locs_len = 2*n;
  int loc_idx = 0;

  // theta.print();
  arma::umat locs(2, locs_len);

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      // Rcout << locs_len << "\n";
      // Rcout << z(i) << z(j) << "\n";
      // Rcout << loc_idx << " " << locs_len << "\n";
      double kern_sim = exp(-pair_dist2(z.row(i).t(), z.row(j).t()));
      double p = scale*kern_sim*theta(i)*theta(j);
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
arma::sp_mat sample_dcer_cpp(arma::vec theta) {
  int n = theta.n_elem;
  int locs_len = 2*n;
  int loc_idx = 0;

  // theta.print();
  arma::umat locs(2, locs_len);

  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      double p = theta(i)*theta(j);
      if (R::runif(0,1) < p) {
        if (loc_idx >= locs_len) {
          locs_len = round(locs_len*1.5);
          locs.resize(2, locs_len);
        }
        locs(0, loc_idx) = i;
        locs(1, loc_idx) = j;
        loc_idx++;
      }
    }
  }
  locs.resize(2, loc_idx);

  locs = arma::join_horiz(locs, arma::join_vert(locs.row(1),locs.row(0)));

  arma::sp_mat A(locs, arma::ones<arma::vec>(locs.n_cols), n, n);
  return A;
}



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
