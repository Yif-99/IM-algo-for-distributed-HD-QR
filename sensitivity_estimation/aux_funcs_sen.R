#========functions=============
check <- function(x,tau){
  x*(tau - (x<0))
}
#========
HBIC_conquer_ini <- function(lambda,pen,X,y){
  n = length(y)
  result = conquer.reg(X,y,lambda=lambda,tau=tau,kernel='uniform',penalty=pen
                       ,epsilon=0.001,iteTight=3)
  checkloss = sum(check(y-X%*%result$coeff[-1],tau=tau))
  hbic = log(checkloss) + sum(result$coeff[-1]!=0)*log(log(n))/n*1.5*log(p)
  return(list(hbic = hbic, model=result))
}
#
ini_est <- function(lam_seq){
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_conquer_ini(lambda,"scad",X_master,y_master)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist = ret_conquer_seq[2,min_idx_conquer]$model
  lambda_ini = result_dist$lambda  
  beta_ini = result_dist$coeff[-1]
  return(beta_ini)
}
#========
HBIC_global <- function(lambda,betaini,pen="Lasso",verbose=FALSE){
  if(pen == "Lasso"){
    result = dist_SQR_uniform_L1(X,y,tau,beta=betaini,M=1,lambda
                                 ,etag,G,phi=1.0,itermax=300,iterT=1,verbose=verbose)
  }else if(pen=="SCAD"){
    result = dist_SQR_uniform_scad(X,y,tau,beta=betaini,M=1,lambda,etag,G,ifora=FALSE,ora_supp=1*(beta0 != 0)
                                   ,phi=1,itermax=300,iterT=4,verbose=verbose)
  }
  checkloss = sum(check(y-X%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1.5*log(p)
  return(list(hbic = hbic, model=result))
}
#
global_est  <- function(lam_seq,beta_pre, pen="Lasso",verbose=FALSE){
  clusterExport(cl7, list("beta_pre"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_global(lambda,beta_pre,pen,verbose)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  return(ret_conquer_seq[2,min_idx_conquer]$model)
}
#========
HBIC_distlasso <- function(c0=0.5,lambda,betaini,T,verbose=FALSE){
  h = c0 * (hat_s_ini*log(p)/N)^(0.25);
  b = c0 * (hat_s_ini*log(p)/n)^(0.25);
  result = dist_SQR_uniform_L1(X,y,tau,beta=betaini,m,lambda,eta,G
                               ,s=hat_s_ini,h=h,b=b,c0=c0,phi=1.0, itermax=200,iterT=T,verbose)
  checkloss = sum(check(y-X%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1.5*log(p)
  return(list(hbic = hbic, model=result))
}
#
dist_l1_est <- function(c0=0.5,lam_seq,T,beta_pre){
  clusterExport(cl7, list("beta_pre","hat_s_ini"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lambda) {
    HBIC_distlasso(c0,lambda,betaini=beta_pre,T)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist_lasso = ret_conquer_seq[2,min_idx_conquer]$model
  return(result_dist_lasso)
}
#================
HBIC_distscad <- function(c0=0.5,lambda,betaini,T,verbose=FALSE){
  h = c0 * max( c(hat_s_ini^(1.5) * log(p)/N , (hat_s_ini^(2)*log(hat_s_ini)/N)^(0.5)) )
  b = c0 * max( c((hat_s_ini*log(p)/N)^(0.5),  hat_s_ini^1.5 *log(p)/n,  (hat_s_ini^2*log(hat_s_ini)/n)^(0.5))  )
  result = dist_SQR_uniform_scad(X,y,tau,beta=betaini,m,lambda,eta,G
                                 ,s=hat_s_ini,h=h,b=b,c0=c0,ifora=FALSE,ora_supp=1*(beta0 != 0)
                                 ,phi=1.0,itermax=300,iterT=T,verbose=verbose)
  checkloss = sum(check(y-X%*%result$coeff,tau=tau))
  hbic = log(checkloss) + sum(result$coeff!=0)*log(log(N))/N*1.5*log(p)
  return(list(hbic = hbic, model=result))
}
#
dist_scad_est <- function(c0=0.5,lam_seq,T,beta_pre){
  clusterExport(cl7, list("beta_pre","hat_s_ini"))
  ret_conquer_seq = NULL
  ret_conquer_seq <- do.call(cbind, parLapply(cl7, lam_seq, function(lam) {
    HBIC_distscad(c0,lambda=lam,betaini=beta_pre,T=T)
  }))
  min_idx_conquer = which.min(ret_conquer_seq[1,])
  result_dist_scad = ret_conquer_seq[2,min_idx_conquer]$model
  return(result_dist_scad)
}

###########################################
library(Rcpp);library(RcppArmadillo)
library("conquer");library("parallel")
library(foreach);library(doParallel)
sourceCpp('../functions/dist_HD_QR.cpp'); sourceCpp('../functions/generate_data.cpp')
####################################################
G = 20; eta = rep(7.5,G) ; etag = rep(4,G) 
############## PARALLEL COMPUTATION #################
cl7 <- makeCluster(2); clusterEvalQ(cl7, library("conquer"))
vars_to_export <- list("HBIC_distlasso","HBIC_distscad","HBIC_global", "HBIC_conquer_ini"
                       ,'check','beta0',"G", "eta","etag", "p",  'check')
clusterExport(cl7, vars_to_export); clusterEvalQ(cl7, Rcpp::sourceCpp('../functions/dist_HD_QR.cpp'))
clusterExport(cl7, list('tau','n',"m", "N"))
####################################################
