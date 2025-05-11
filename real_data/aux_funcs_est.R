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
  result = dist_SQR_uniform_scad(X,y,tau,beta=beta_pre,M=1,lambda,eta,G
                                 ,s=-1,h=-1,b=-1,c0=c0,ifora=FALSE,ora_supp=1*(beta_pre != 0)
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
#========
HBIC_distlasso <- function(c0=0.8,lambda,betaini,iterT,verbose=FALSE){
  h = c0 * (hat_s_ini*log(p)/N)^(0.25);
  b = c0 * (hat_s_ini*log(p)/n)^(0.25);
  result = dist_SQR_uniform_L1(X_a_train,y_train,tau,beta=betaini,m,lambda,eta,G
                               ,s=hat_s_ini,h=h,b=b,c0=c0,phi=1.0/p, itermax=3000,iterT=iterT,verbose)
  checkloss = sum(check(y_train-X_a_train%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1*log(p)
  return(list(hbic = hbic, model=result))
}
dist_l1_est <- function(c0=0.8,lam_seq,iterT,beta_pre){
  clusterExport(cl7, list("beta_pre","hat_s_ini" ))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_distlasso(c0,lambda,betaini=beta_pre,iterT)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist_lasso = ret_conquer_seq[2,min_idx_conquer]$model
  return(result_dist_lasso)
}
#========
#========
HBIC_distscad <- function(c0=0.8,lambda,betaini,iterT,verbose=FALSE){
  h =  c0 * max( c(hat_s_ini^(1.5) * log(p)/N , (hat_s_ini^(2)*log(hat_s_ini)/N)^(0.5)) )
  b =  min(0.12,c0 * max( c((hat_s_ini*log(p)/N)^(0.5),  hat_s_ini^1.5 *log(p)/n,  (hat_s_ini^2*log(hat_s_ini)/n)^(0.5))  ) )
  result = dist_SQR_uniform_scad(X_a_train,y_train,tau,beta=betaini,m,lambda,eta,G
                                 ,s=hat_s_ini,h=h,b=b,c0=c0,ifora=FALSE,ora_supp=1*(betaini != 0)
                                 ,phi=1.0/p,itermax=30000,iterT=iterT,verbose=verbose)
  checkloss = sum(check(y_train-X_a_train%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1*log(p)
  return(list(hbic = hbic, model=result))
}

dist_scad_est <- function(c0=0.8,lam_seq,iterT,beta_pre){
  clusterExport(cl7, list("beta_pre","hat_s_ini"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lam) {
    HBIC_distscad(c0,lambda=lam,betaini=beta_pre,iterT=iterT)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist_scad = ret_conquer_seq[2,min_idx_conquer]$model
  return(result_dist_scad)
}
#========
HBIC_global <- function(lambda,betaini,pen="Lasso",verbose=FALSE){
  X_a_train = cbind(rep(1,nrow(X_train)),X_train)
  if(pen == "Lasso"){
    result = dist_SQR_uniform_L1(X_a_train,y_train,tau,beta=betaini,M=1,lambda
                                 ,etag,G,phi=1.0/nrow(X_train),itermax=3000,iterT=1,verbose=verbose)
  }else if(pen=="SCAD"){
    result = dist_SQR_uniform_scad(X_a_train,y_train,tau,beta=betaini,M=1,lambda,etag,G,ifora=FALSE,ora_supp=1*(betaini!= 0)
                                   ,phi=1/nrow(X_train),itermax=3000,iterT=4,verbose=verbose)
  }
  checkloss = sum(check(y_train-X_a_train%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1*log(p)
  return(list(hbic = hbic, model=result))
}

global_est  <- function(lam_seq,beta_pre, pen="Lasso",verbose=FALSE){
  clusterExport(cl7, list("beta_pre"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_global(lambda,beta_pre,pen,verbose)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  return(ret_conquer_seq[2,min_idx_conquer]$model)
}


###########################
library(Rcpp)
library(RcppArmadillo)
library("conquer")
library("parallel")
library(foreach)
library(doParallel)
library(data.table)
sourceCpp('../functions/dist_HD_QR.cpp')
sourceCpp('firsteig.cpp')
