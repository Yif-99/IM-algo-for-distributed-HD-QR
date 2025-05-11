source("aux_funcs_est.R")
X <- fread("suvtrain.csv", header = TRUE) # full model
# X <- fread("suvtrain_reduced.csv", header = TRUE)  # reduced model

X = as.matrix(X) ; y=log(X[,1]); X = X[,-1]
X[ ,1:15] = scale(X[ ,1:15])

set.seed(42); random_order = sample(nrow(X)); X = X[random_order,]; y = y[random_order]
N = 300000   
p = ncol(X)
X_train = X[1:N,];X_pred = X[-(1:N),];
y_train = y[1:N]; y_pred = y[-(1:N)]; N_pred = nrow(X_pred)
X_a <- cbind(rep(1,nrow(X)), X); X_a_train = X_a[1:N,]; X_a_pred = X_a[-(1:N),];

m = 80; n = N/m
X_master = X_train[1:n,]; y_master = y_train[1:n]

G = 5; eta = first_eigenvalue(X_a_train[1:n,],G,(p+1)/G)/p  +  1
etag = first_eigenvalue(X_a_train[ ,],G,(p+1)/G)/N + 1
c0 = 0.8
ifglobal = 1; n_cols = 5; source('store_matrix.R')

############## PARALLEL COMPUTATION #################
cl7 <- makeCluster(2); clusterEvalQ(cl7, library("conquer"))
clusterExport(cl7, list( "X_master","y_master" ,"X_train","X_a_train","y_train" ))
clusterExport(cl7, list("eta" , 'check','c0', "m", "G", "p", "n", "N"
                        , "HBIC_distlasso","HBIC_distscad",'HBIC_global',"HBIC_conquer_ini"))
if(ifglobal){
  clusterExport(cl7, list("etag"))
}
clusterEvalQ(cl7, Rcpp::sourceCpp('../functions/dist_HD_QR.cpp'))
####################################################



for(tau in c(0.1,0.3,0.5,0.7,0.9)){  
  
  if(tau==0.5){inside=3}else if(tau==0.7){inside=4}else if(tau==0.3){inside=2}else if(tau==0.1){inside=1}else if(tau==0.9){inside=5}
  clusterExport(cl7, list("tau", 'inside') ); 
  ### Local initialization
  beta_pre =  c(quantile(y_train,tau),rep(0,p))
  if(m<=30){
    lam_seq = seq(0.02, 0.001, -0.001) 
  }else{
    lam_seq = seq(0.01, 0.0001, -0.0005) 
  }
  result_ini = ini_est(lam_seq,iterT = 6)
  beta_ini=result_ini$coef; lambda_ini[inside] = result_ini$lambda
  names(beta_ini) = c("Intercept",names(X_train[1,]))
  coef_ini_record[,inside] = beta_ini
  size_ini[ inside] = result_ini$size;  hat_s_ini = min(c( size_ini[ inside],10))
  check_pred_ini[inside] = check_error(beta_ini,X_pred,y_pred)
  ###global
  lam_seq = seq(0.0008, 0.0001, -0.0001) 
  beta_pre = beta_ini;
  result_global_scad  = global_est(lam_seq,beta_pre,verbose=FALSE)
  lambda_global_lasso[ inside] = result_global_scad$pen_lambda
  names(result_global_scad$coeff) = c("Intercept",names(X_train[1,]))
  coef_global_record_lasso[,inside] = result_global_scad$coeff
  size_global_lasso[inside] = sum(result_global_scad$coeff!=0)
  check_pred_global_lasso[inside] = check_error(result_global_scad$coeff,X_pred,y_pred)
  beta_pre =  coef_global_record_lasso[,inside];
  result_global_scad  = global_est(lambda_global_lasso[ inside],beta_pre,pen = 'SCAD',verbose=FALSE)
  names(result_global_scad$coeff) = c("Intercept",names(X_train[1,]))
  coef_global_record [,inside] = result_global_scad$coeff
  lambda_global_scad[ inside] = result_global_scad$pen_lambda
  size_global_scad[inside] = sum(result_global_scad$coeff!=0)
  check_pred_global[inside] = check_error(result_global_scad$coeff,X_pred,y_pred)
  ### Contraction Stage (iterT=1)
  if(m<=50){
    lam_seq =  seq(0.001, 0.0001, -0.0001)
  }else{
    lam_seq =  seq(0.002, 0.0001, -0.0001)
  }
  beta_pre = beta_ini; result_dist_lasso = dist_l1_est(c0,lam_seq,iterT=1,beta_pre)
  lambda_lasso1[inside] = result_dist_lasso$pen_lambda
  size_dist_lasso[inside] = sum(result_dist_lasso$coeff!=0)
  names(result_dist_lasso$coeff) = c("Intercept",names(X[1,]))
  check_pred_lasso1[inside] = check_error(result_dist_lasso$coeff,X_pred,y_pred)
  k_count = 1
  while (k_count < 6 ) {
    k_count = k_count + 1;
    result_old = result_dist_lasso
    beta_pre = result_old$coeff
    result_dist_lasso =  dist_l1_est(c0,lam_seq,iterT=1,beta_pre)
  }
  lambda_lasso5[inside] = result_dist_lasso$pen_lambda
  size_dist_lasso[inside] = sum(result_dist_lasso$coeff!=0)
  coeflasso_record[,inside] = result_dist_lasso$coeff
  names(result_dist_lasso$coeff) = c("Intercept",names(X[1,]))
  check_pred_lasso[inside] = check_error(result_dist_lasso$coeff,X_pred,y_pred)
  ### Tightening Stage (iterT>=2)
  beta_pre = coeflasso_record[,inside];
  result_dist_scad = dist_scad_est(c0,lam_seq,iterT=1,beta_pre )
  lambda_dist1[inside] = result_dist_scad$pen_lambda
  size_dist_scad1[inside] = sum(result_dist_scad$coeff!=0)
  names(result_dist_scad$coeff) = c("Intercept",names(X[1, ]))
  t_count = 2;
  while (t_count <= 5) {
    t_count = t_count + 1;
    result_old = result_dist_scad
    beta_pre = result_old$coeff
    result_dist_scad =  dist_scad_est(c0,lam_seq,iterT=1,beta_pre)
  }
  size_dist_scad5[inside] = sum(result_dist_scad$coeff!=0)
  names(result_dist_scad$coeff) = c("Intercept",names(X[1, ]))
  coef_record[,inside] = result_dist_scad$coeff
  check_pred_scad[inside] = check_error(result_dist_scad$coeff,X_pred,y_pred)
}

check_pred_null = check_pred_ini
for(tau in c(0.1,0.3,0.5,0.7,0.9)){  
  if(tau==0.5){inside=3}else if(tau==0.7){inside=4}else if(tau==0.3){inside=2}else if(tau==0.1){inside=1}else if(tau==0.9){inside=5}
  beta_pre =  c(quantile(y_train,tau),rep(0,p))
  check_pred_null[inside] = check_error(beta_pre,X_pred,y_pred)
  
}
  

 
stopCluster(cl7)

rm("X","y","X_train","y_train","X_pred","y_pred","X_a",'X_a_pred',"X_a_train", "cl7","X_master","y_master", "beta_ini" 
   ,"result_dist_scad","result_dist_lasso" ,"beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove" )
save.image(file=paste0("Results/real_est_n_",n,"_m_",m,"_p_",p,".RData"))



