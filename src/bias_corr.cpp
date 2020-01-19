#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec group_sums(const arma::vec &kv,
                     const arma::ivec &kT) {
  // Auxiliary variables
  const int kn = kT.n_rows;
  
  // Compute group sums
  arma::vec sum(kn, arma::fill::zeros);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute individual sum
    for (int t = 0 ; t < kT(i) ; ++t) {
      sum(i) += kv(j + t);
    }
    
    // Increase counter
    j += kT(i);
  }
  
  // Return vector
  return sum;
}


// [[Rcpp::export]]
arma::vec group_sums_bias(const arma::mat &kM,
                          const arma::vec &kw,
                          const arma::ivec &kT) {
  // Auxiliary variables
  const int kp = kM.n_cols;
  const int kn = kT.n_rows;
  
  // Compute sum of weighted group sums
  arma::vec wsum(kp, arma::fill::zeros);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute numerator and denominator
    arma::rowvec num(kp, arma::fill::zeros);
    double denom = 0.0;
    for (int t = 0 ; t < kT(i) ; ++t) {
      num += kM.row(j + t);
      denom += kw(j + t);
    }
    
    // Add weighted group sum
    wsum += (num / denom).t();
    
    // Increase counter
    j += kT(i);
  }
  
  // Return vector
  return wsum;
}


// [[Rcpp::export]]
arma::mat group_sums_cov(const arma::mat &kM1,
                         const arma::mat &kM2,
                         const arma::ivec &kT) {
  // Auxiliary variables
  const int kp = kM1.n_cols;
  const int kn = kT.n_rows;
  
  // Compute covariance matrix
  arma::mat V(kp, kp, arma::fill::zeros);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Add group specific covariance
    arma::rowvec sum(kp, arma::fill::zeros);
    for (int t = 0 ; t < kT(i) ; ++t) {
      for (int s = t + 1 ; s < kT(i) ; ++s) {
        V += kM1.row(j + t).t() * kM2.row(j + s);
      }
    }
    
    // Increase counter
    j += kT(i);
  }
  
  // Return matrix
  return V;
}


// [[Rcpp::export]]
arma::vec group_sums_spectral(const arma::mat &kM,
                              const arma::vec &kv,
                              const arma::vec &kw,
                              const int kL,
                              const arma::ivec &kT) {
  // Auxiliary variables
  const int kp = kM.n_cols;
  const int kn = kT.n_rows;
  
  // Compute sum of weighted group sums
  arma::vec wsum(kp, arma::fill::zeros);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute numerator
    arma::rowvec num(kp, arma::fill::zeros);
    for (int l = 1 ; l <= kL ; ++l) {
      for (int t = l ; t < kT(i) ; ++t) {
        num += kM.row(j + t) * kv(j + t - l) * kT(i) / (kT(i) - l);
      }
    }
    
    // Compute denominator
    double denom = 0.0;
    for (int t = 0 ; t < kT(i) ; ++t) {
      denom += kw(j + t);
    }
    
    // Add weighted group sum
    wsum += (num / denom).t();
    
    // Increase counter
    j += kT(i);
  }
  
  // Return vector
  return wsum;
}


// [[Rcpp::export]]
arma::mat group_sums_var(const arma::mat &kM,
                         const arma::ivec &kT) {
  // Auxiliary variables
  const int kp = kM.n_cols;
  const int kn = kT.n_rows;
  
  // Compute covariance matrix
  arma::mat V(kp, kp, arma::fill::zeros);
  int j = 0;
  for (int i = 0 ; i < kn ; ++i) {
    // Compute group sum
    arma::rowvec sum(kp, arma::fill::zeros);
    for (int t = 0 ; t < kT(i) ; ++t) {
      sum += kM.row(j + t);
    }
    
    // Add group specific covariance
    V += sum.t() * sum;
    
    // Increase counter
    j += kT(i);
  }
  
  // Return matrix
  return V;
}
