#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat center_variables(const arma::mat &kM,
                           const arma::vec &kw,
                           const arma::ivec &kT) {
  // Auxiliary variables
  const int knt = kM.n_rows;
  const int kp = kM.n_cols;
  const int kn = kT.n_rows;
  
  // Generate matrix of centered variables
  arma::mat M(knt, kp);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute weighted group mean
    arma::rowvec fac(kp, arma::fill::zeros);
    double denom = 0.0;
    for (int t = 0 ; t < kT(i) ; ++t) {
      fac += kM.row(j + t) * kw(j + t);
      denom += std::pow(kw(j + t), 2);
    }
    
    // Center variables
    fac /= denom;
    for (int t = 0 ; t < kT(i) ; ++t) {
      M.row(j + t) = kM.row(j + t) - kw(j + t) * fac;
    }
    
    // Increase counter
    j += kT(i);
  }
  
  // Return matrix of centered variables
  return M;
}


// [[Rcpp::export]]
arma::vec update_alpha(const arma::vec &ku,
                       const arma::vec &kw,
                       const arma::ivec &kT) {
  // Auxiliary variable
  const int kn = kT.n_rows;
  
  // Update incidental parameters given update of structural parameters 
  arma::vec alpha(kn);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute weighted individual mean
    double num = 0.0;
    double denom = 0.0;
    for (int t = 0 ; t < kT(i) ; ++t) {
      num += ku(j + t) * kw(j + t);
      denom += std::pow(kw(j + t), 2);
    }
    
    // Store update
    alpha(i) = num / denom;
    
    // Increase counter
    j += kT(i);
  }
  
  // Return matrix of centered variables
  return alpha;
}


// [[Rcpp::export]]
arma::vec variance_alpha(const arma::mat &kV,
                         const arma::mat &kX,
                         const arma::vec &kw,
                         const arma::ivec &kT) {
  // Auxiliary variable
  const int kp = kX.n_cols;
  const int kn = kT.n_rows;
  
  // Update incidental parameters given update of structural parameters 
  arma::vec var_alpha(kn);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute weighted individual mean
    arma::rowvec num(kp, arma::fill::zeros);
    double denom = 0.0;
    for (int t = 0 ; t < kT(i) ; ++t) {
      num += kX.row(j + t) * kw(j + t);
      denom += kw(j + t);
    }
    
    // Compute standard error
    const arma::rowvec kj = num / denom;
    var_alpha(i) = 1.0 / denom + arma::as_scalar(kj * kV * kj.t());
    
    // Increase counter
    j += kT(i);
  }
  
  // Return matrix of centered variables
  return var_alpha;
}
