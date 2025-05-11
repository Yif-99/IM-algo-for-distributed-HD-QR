#========functions=============
check <- function(x,tau){
  x*(tau - (x<0))
}
#========
check_error <- function(beta,X,y){
  N = nrow(X)
  X_a <- cbind(rep(1,N), X)
  return( sum(check(y-X_a%*%beta,tau=tau))/N )
} 
#========
HBIC_conquer_ini <- function(lambda,beta_pre,X,y,iterT=4){
  n = length(y)
  X <- cbind(rep(1,n), X)
  result = dist_SQR_uniform_scad(X,y,test_v,tau,beta=beta_pre,M=1,lambda,eta,G
                                  ,h=-1,b=-1, ifora=FALSE,ora_supp=1*(beta_pre != 0)
                                 ,phi=1.0/p,itermax=3000,iterT=iterT,verbose=FALSE)
  checkloss = sum(check(y-X%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff !=0)*log(log(n))/n*1*log(p)  
  return(list(hbic = hbic, model=result))
}
ini_est <- function(lam_seq,iterT=4){
  clusterExport(cl7, list("beta_pre"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_conquer_ini(lambda,beta_pre,X_master,y_master,iterT)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result = ret_conquer_seq[2,min_idx_conquer]$model
  lambda_ini = result$pen_lambda  
  beta_ini = result$coeff
  return(list( coef=result$coeff,lambda=lambda_ini,size= sum( beta_ini!=0) ))
}
#========
HBIC_distlasso <- function(c0=0.5,lambda,betaini,iterT,verbose=FALSE){
  h = c0 * (hat_s_ini*log(p)/N)^(0.25);
  b = c0 * (hat_s_ini*log(p)/n)^(0.25);
  result = dist_SQR_uniform_L1(X_a_train,y_train,test_v,tau,beta=betaini,m,lambda,eta,G
                               ,s=hat_s_ini,h=h,b=b,c0=c0,phi=1.0/p, itermax=3000,iterT=iterT,verbose)
  checkloss = sum(check(y_train-X_a_train%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1*log(p)
  return(list(hbic = hbic, model=result))
}
#
dist_l1_est <- function(c0=0.5,lam_seq,iterT,beta_pre){
  clusterExport(cl7, list("beta_pre","hat_s_ini"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_distlasso(c0,lambda,betaini=beta_pre,iterT)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist_lasso = ret_conquer_seq[2,min_idx_conquer]$model
  return(result_dist_lasso)
}
#================
HBIC_distscad <- function(c0=0.5,lambda,betaini,iterT,verbose=FALSE){
  h = c0 * max( c(hat_s_ini^(1.5) * log(p)/N , (hat_s_ini^(2)*log(hat_s_ini)/N)^(0.5)) )
  b =  min(0.12,c0 * max( c((hat_s_ini*log(p)/N)^(0.5),  hat_s_ini^1.5 *log(p)/n,  (hat_s_ini^2*log(hat_s_ini)/n)^(0.5))  ) )
  result = dist_SQR_uniform_scad(X_a_train,y_train,test_v,tau,beta=betaini,m,lambda,eta,G
                                  ,h=h,b=b, ifora=FALSE,ora_supp=1
                                 ,phi=1.0/p,itermax=3000,iterT=iterT,verbose=verbose)
  checkloss = sum(check(y_train-X_a_train%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1*log(p)
  return(list(hbic = hbic, model=result))
}
#
dist_scad_est <- function(c0=0.5,lam_seq,iterT,beta_pre){
  clusterExport(cl7, list("beta_pre","hat_s_ini"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lam) {
    HBIC_distscad(c0,lambda=lam,betaini=beta_pre,iterT)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist_scad = ret_conquer_seq[2,min_idx_conquer]$model
  return(result_dist_scad)
}
###########################
sigma_dist = function(X,y,beta,k,s){
  H = matrix(0,nrow =s ,ncol = s )
  e = y -X%*%beta;
  for(ii in 1:n){
    H = H + dnorm(e[ii]/k)*  X[ii, ] %*% t(X[ii, ]) / (n*k)
  }
  if(qr(H)$rank < 1*s){
    count = 0; sigma2 = 0
  }else{
    count = 1; sigma2 = sigmahat(X,y ,beta, v_test, tau , k=k, h=h ) 
  }
  return(list(count=count,sigma2=sigma2))
}
sigma_est <- function(X,y,beta,k,s,m){
  clusterExport(cl7, list('ck','test_v','tau','hat_s_ini','result_dist_scad', 'v_test','size_dist_scad5','support',"h",'inside'))
  seq = NULL
  seq <- do.call(cbind, parLapply(cl7, 1:m, function(site) {
    sigma_dist(X[(n*(site-1)+1):(n*site),],y[(n*(site-1)+1):(n*site)],beta,k,s)
  }))
  count = sum(as.numeric(seq[1,]))
  sigma2=0
  for(site in 1:m){
    sigma2 = sigma2 +   (seq[2,][[site]])/count
  }
  return(list(sigma2=sigma2,count=count))
}
###########################
library(Rcpp)
library(RcppArmadillo)
library("conquer")
library("parallel")
library(foreach)
library(doParallel)
library(data.table)
sourceCpp('../functions/dist_HD_QR_inference.cpp'); sourceCpp('firsteig.cpp')
