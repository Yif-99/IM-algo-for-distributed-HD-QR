// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat generateData( int N, int p, double rho, double B) {
  arma::mat X(N, p+1, arma::fill::zeros);
  X.col(0) = arma::randn<arma::vec>(N);
  
  for (int i = 1; i < p; ++i) {
    X.col(i) = rho * X.col(i - 1) + std::sqrt(1 - std::pow(rho, 2)) * arma::randn<arma::vec>(N);
  }
  X.col(0) = normcdf(X.col(0));
  X = clamp(X, -B, B);
  X.col(p) = arma::randn<arma::vec>(N);
  return X;
}
