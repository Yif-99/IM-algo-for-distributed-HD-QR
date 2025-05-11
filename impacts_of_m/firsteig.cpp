#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
using namespace arma;
using namespace RcppParallel;



// [[Rcpp::export]]
vec first_eigenvalue(mat A,int K,int p_g){
  vec eigen(K);
  vec storage;
  for(int i = 0; i < K; i++){
    storage = svd(A.cols(i*p_g,(i+1)*p_g-1));
    eigen(i) = std::pow(storage(0),2);
    storage.reset();
  }
  return eigen;
}