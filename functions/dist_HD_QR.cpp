#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
using namespace arma;
using namespace RcppParallel;


vec Shrink(vec z, const vec& rho){
  int n = z.n_elem;
  for(int i=0; i<n; i++){
    z(i) = sign(z(i))*std::max(abs(z(i))-rho(i),0.0);
  }
  return z;
}

vec pen_deriv_scad(vec beta,const double a,const double lambda,const int p){
  for(int i=0; i < p; i++){
    if(beta(i)<=lambda){
      beta(i) = lambda;
    }
    else{
      beta(i) = std::max(0.0,(a*lambda-beta(i))/(a-1));
    }
  }
  return beta;
}

inline double K_bar_uniform(double input){
  double output;
  if(input < (-1)){
    output = 0.0;
  }else if(input > 1){
    output = 1.0;
  }else{
    output = (input+1)/(2.0);
  }
  return output;
}

vec z_uni_update(vec b,double tau,double phi,int n,double h){
  for(int i =0;i<n;i++){
    if(b(i) < ((tau-1)/(n*phi)-h) ){
      b(i) += (1-tau)/(n*phi);
    }else if(b(i) > (tau/(n*phi)+h) ){
      b(i) -= tau/(n*phi);
    }else{
      b(i) = (2*n*phi*h*b(i) + h - 2*tau*h)/(1+2*n*phi*h);
    }
  }
  return b;
}


struct GradParallel : public Worker {
  const arma::mat& X;
  const arma::vec& beta;
  const arma::vec& y;
  const double h;
  const double tau;
  const int N;
  arma::vec value;
  GradParallel(const arma::mat& X, const arma::vec& beta, const arma::vec& y, double h, double tau, int N)
    : X(X), beta(beta), y(y), h(h), tau(tau), N(N), value(zeros<vec>(beta.n_elem)) {}
  GradParallel(const GradParallel& gradParallel, Split)
    : X(gradParallel.X), beta(gradParallel.beta), y(gradParallel.y), h(gradParallel.h), tau(gradParallel.tau)
    , N(gradParallel.N), value(zeros<vec>(gradParallel.beta.n_elem)) {}
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      arma::rowvec xi = X.row(i);
      value += (K_bar_uniform(as_scalar(xi*beta-y[i])/h) - tau) * xi.t();
    }
  }
  void join(const GradParallel& rhs) {
    value += rhs.value;
  }
};
 
arma::vec parallelGrad(const arma::mat& X, const arma::vec& beta, const arma::vec& y, double h, double tau) {
  int N = X.n_rows;
  GradParallel gradParallel(X, beta, y, h, tau, N);
  parallelReduce(0, N, gradParallel);
  gradParallel.value /= N;
  return gradParallel.value;
}


vec SQR_local_L1(const mat& X,const vec& y,vec beta,const double tau,const double h
                   ,double pen_lambda,const vec& gamma,const vec& eta,const int G,double phi
                   ,const int itermax, bool verbose = false){
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int p_g = p/G;
  const vec constzeros(n,fill::zeros);
  double epsa = 0.0005;  double epsr = 0.0005;
  double pri_res=1;  double dual_res=1;

  vec u(n,fill::zeros); 
  vec w = u;
  vec xbeta = X*beta/G;
  w = ( y - z_uni_update(y-G*xbeta-G*u,tau,phi/G,n,h) )/G; 
  u += xbeta - w; 
  vec w_old = w;
  vec lambda_weight(p,fill::ones);
  lambda_weight *= pen_lambda;
  
  int iter1 = 1;
  while(iter1 <= itermax && (pri_res > epsa || dual_res > epsr) ){
    w_old = w;
    for(int i = 0; i < G; i++){
      beta.rows(p_g*i,p_g*(i+1)-1) = Shrink(beta.rows(p_g*i,p_g*(i+1)-1)- ((1.0)/eta(i))*phi
                          *X.cols(p_g*i,p_g*(i+1)-1).t()*(xbeta-w+u)-gamma.rows(p_g*i,p_g*(i+1)-1)/eta(i),((1.0)/eta(i))*lambda_weight.rows(p_g*i,p_g*(i+1)-1));
    }
    xbeta = X*beta/G;
    w = ( y - z_uni_update(y-G*xbeta-G*u,tau,phi/G,n,h) ) / G;
    u += xbeta - w;
    
    dual_res = phi*norm((w_old-w),2);
    pri_res = sqrt(G)*norm(xbeta - w,2);
    iter1++;
  }
  return beta;
}



// [[Rcpp::export]]
Rcpp::List dist_SQR_uniform_L1(const mat& X,const vec& y,const double tau,vec beta
                               ,const int M, double pen_lambda,const vec& eta, const int G
                               ,int s=10 ,double h = -1, double b =-1,double c0=0.8
                             ,double phi=1.0,const int itermax=100, int iterT=10, bool verbose=false){
  
  const int N = X.n_rows;
  const int p = X.n_cols;
  const int n = N/M;
  if(M == 1){
    // non-distributed setting
    iterT = 1;
    h = std::max(0.05,0.5*std::pow(log(p)/N,0.25));
    b = h;
  }else{
    if(s>15 || h<0){
      s = std::min(s,15);
      h = c0*std::pow(s*log(p)/N,0.25);
    }
    if(s>15 || b<0){
      b = c0*std::pow(s*s*log(p)/n,0.25);
    }
  }

  mat X_local = X.rows(0,n-1);
  vec y_local = y.rows(0,n-1);
  vec grad_global(p,fill::zeros);
  vec grad_local(p,fill::zeros);
  int t = 0;
  while(t < iterT){
    if(M > 1){
      grad_global = parallelGrad(X,beta,y,h,tau);
      grad_local = parallelGrad(X_local,beta,y_local,b,tau);
    }
    beta = SQR_local_L1(X_local,y_local,beta,tau,b,pen_lambda
                          ,grad_global-grad_local,eta,G,phi,itermax,verbose);
    t++;
  }
  return  Rcpp::List::create(Rcpp::Named("coeff") = beta
                               ,Rcpp::Named("iter_T") = t
                               ,Rcpp::Named("pen_lambda") = pen_lambda
                               ,Rcpp::Named("global_bandwidth") = h
                               ,Rcpp::Named("local_bandwidth") = b
  );
}



vec SQR_local_scad(const mat& X,const vec& y,const double tau,vec beta,const double h
                     ,const vec& lambda_weight,const vec& gamma,const vec& eta,const int G
                       ,double phi,const int itermax, bool verbose = false){
  // initialize parameter and data
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int p_g = p/G;
  const vec constzeros(n,fill::zeros);
  double epsa = 0.0005;  double epsr = 0.0005;
  double pri_res=1;  double dual_res=1;

  vec u(n,fill::zeros); 
  vec w = u;
  vec w_old = w;
  vec xbeta = X*beta/G;
  w = ( y - z_uni_update(y-G*xbeta-G*u,tau,phi/G,n,h) ) / G;
  u += xbeta - w;
  int iter2 = 1;
  while(iter2 <= itermax && (pri_res > epsa || dual_res > epsr) ){
    w_old = w;
    for(int i = 0; i < G; i++){
        beta.rows(p_g*i,p_g*(i+1)-1) = Shrink(beta.rows(p_g*i,p_g*(i+1)-1)- ((1.0)/eta(i))*phi
                         *X.cols(p_g*i,p_g*(i+1)-1).t()*(xbeta-w+u)-gamma.rows(p_g*i,p_g*(i+1)-1)/eta(i)
                               ,((1.0)/eta(i))*lambda_weight.rows(p_g*i,p_g*(i+1)-1));
    }
    xbeta = X*beta/G;
    w = ( y - z_uni_update(y-G*xbeta-G*u,tau,phi/G,n,h) ) / G;
    u += xbeta - w;
    //==========residual computation=========
    dual_res = phi*norm((w_old-w),2);
    pri_res = sqrt(G)*norm(xbeta - w,2);
    iter2++;
  }
  return beta;
}




vec SQR_local_ora(const mat& X,const vec& y ,const double tau,vec beta,const double h
                  ,const vec& gamma,double phi, int itermax, bool verbose=false){
  // initialize parameter and data
  const int n = X.n_rows;
  double epsa = 0.0005;  double epsr = 0.0005;
  double pri_res=1;  double dual_res=1;
  
  vec u(n,fill::zeros); 
  vec w = u;
  vec w_old = w;
  vec xbeta = X*beta;
  w = ( y - z_uni_update(y-1*xbeta-1*u,tau,phi/1.0,n,h) ) ;
  u += xbeta - w;
  
  int iter2 = 1;
  while(iter2 <= itermax && (pri_res > epsa || dual_res > epsr) ){
    w_old = w;
    beta =  beta - solve(X.t()*X, X.t()*(xbeta-w+u) + gamma/phi);
    xbeta = X*beta;
    w = ( y - z_uni_update(y-1*xbeta-1*u,tau,phi/1.0,n,h) ) ;
    u += xbeta - w;
    //==========residual computation=========
    dual_res = phi*norm((w_old-w),2);
    pri_res = norm(xbeta - w,2);
    iter2++;
  }
  return beta;
}



// [[Rcpp::export]]
Rcpp::List dist_SQR_uniform_scad(const mat& X,const vec& y,const double tau, vec beta
                                   ,const int M,double pen_lambda,vec& eta,const int G
                                ,int s=10 ,double h = -1, double b =-1, double c0=0.8
                                   ,bool ifora=true,const uvec ora_supp = arma::zeros<uvec>(1)
                                ,double phi=1, int itermax=100,int iterT=10,bool verbose=false){
  const int N = X.n_rows;
  const int p = X.n_cols;
  int n = N/M;
  if(M == 1){
    // non-distributed setting
    h = std::max(0.05,0.5*std::pow(log(p)/N,0.25));
    b = h;
  }else{
    if(s > 15 || h < 0){
      s = std::min(s,15);
      h = c0*std::max(std::pow(s,1.5)*log(p)/N,std::pow(s*s*log(s)/N,0.5) );  
    }
    if(s > 15 || b < 0){
      b = c0*std::max(std::pow(s*log(p)/N,0.5), std::max(std::pow(s,1.5)*log(p)/n,std::pow(s*s*log(s)/n,0.5) ) );
    } 
  }

  vec grad_global(p,fill::zeros);
  vec grad_local(p,fill::zeros);
  mat X_local = X.rows(0,n-1);
  vec y_local = y.rows(0,n-1);
  mat X_local_ora; vec beta_ora;
  if(ifora){
    X_local_ora = X_local.cols(arma::find(ora_supp == 1));
    beta_ora = beta.rows(arma::find(ora_supp == 1));
  }
  vec gamma(p,fill::zeros);
  vec lambda_weight(p,fill::ones);
  
  int t = 0;
  while(t < iterT){
    if(M > 1){
      grad_global = parallelGrad(X,beta,y,h,tau);
      grad_local = parallelGrad(X_local,beta,y,b,tau);
      gamma = grad_global - grad_local;
      if(ifora & (t == (iterT -1))){
        beta_ora = SQR_local_ora(X_local_ora,y_local,tau,beta.rows(arma::find(ora_supp == 1)),b
                                  ,gamma.rows(arma::find(ora_supp == 1)),phi=phi,itermax=500,verbose=verbose);
      }
    }
    lambda_weight = pen_deriv_scad(abs(beta),3.7,pen_lambda,p);
    beta = SQR_local_scad(X_local,y_local,tau,beta,b,lambda_weight
                            ,gamma,eta,G,phi,itermax,verbose=verbose);
    t++;
  }

  if(ifora){
    return  Rcpp::List::create(Rcpp::Named("coeff") = beta
                                 ,Rcpp::Named("coeff_ora") = beta_ora
                                 ,Rcpp::Named("iter_T") = t
                                 ,Rcpp::Named("pen_lambda") = pen_lambda
                                 ,Rcpp::Named("global_bandwidth") = h
                                 ,Rcpp::Named("local_bandwidth") = b
    );
  }else{
    return  Rcpp::List::create(Rcpp::Named("coeff") = beta
                                 ,Rcpp::Named("iter_T") = t
                                 ,Rcpp::Named("pen_lambda") = pen_lambda
                                 ,Rcpp::Named("global_bandwidth") = h
                                 ,Rcpp::Named("local_bandwidth") = b
    );
  }
}

