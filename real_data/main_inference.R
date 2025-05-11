source("aux_funcs_infer.R")
X <- fread("suvtrain.csv", header = TRUE) # full model
# X <- fread("suvtrain_reduced.csv", header = TRUE)  # reduced model

X = as.matrix(X) ; y=log(X[,1]); X = X[,-1]
X[ ,1:15] = scale(X[ ,1:15])

set.seed(42); random_order = sample(nrow(X)); X = X[random_order,]; y = y[random_order]
N =  nrow(X)-1;
p = ncol(X)
X_train = X[1:N,];X_pred = X[-(1:N),];
y_train = y[1:N]; y_pred = y[-(1:N)]; N_pred = nrow(X_pred)
X_a <- cbind(rep(1,nrow(X)), X); X_a_train = X_a[1:N,]; X_a_pred = X_a[-(1:N),];

m = 65; n = N/m  #5 x 13 x 43 x 149
X_master = X_train[1:n,]; y_master = y_train[1:n]

G = 5 ; eta = first_eigenvalue(X_a_train[1:n,],G,(p+1)/G)/p  +  1
c0 = 0.8; ifglobal = 0; n_cols = 3; 
source('store_matrix.R')

test_v = c(5,7) # c(3,9)
# test_v = 13 # 4


pvalue = (rep(0,n_cols) )
T_IM = (rep(0,n_cols) )
sigma_IM = (array(0,c(length(test_v),length(test_v),n_cols)) )


############## PARALLEL COMPUTATION #################
cl7 <- makeCluster(2); clusterEvalQ(cl7, library("conquer"))
clusterExport(cl7, list( "X_master","y_master" ,"X_train","X_a_train","y_train" ))
clusterExport(cl7, list("eta", 'check','c0', "m", "G", "p", "n", "N",'sigma_dist'
                        , "HBIC_distlasso","HBIC_distscad", "HBIC_conquer_ini"))
clusterEvalQ(cl7, Rcpp::sourceCpp('../functions/dist_HD_QR_inference.cpp'))
####################################################


#=======main=======
inside = 0
for(tau in c(0.3,0.5,0.7)){  
  inside = inside + 1
  clusterExport(cl7, list("tau", 'inside','test_v') ); 
  ##########################
  ### Local initialization
  beta_pre =  c(quantile(y_train,tau),rep(0,p))
  if(m<=30){
    lam_seq = seq(0.02, 0.001, -0.001)  
  }else{
    lam_seq = seq(0.02, 0.001, -0.001)  
  }
  result_ini = ini_est(lam_seq,iterT = 6)
  beta_ini=result_ini$coef; lambda_ini[inside] = result_ini$lambda
  names(beta_ini) = c("Intercept",names(X_train[1,]))
  coef_ini_record[,inside] = beta_ini
  size_ini[ inside] = result_ini$size;  hat_s_ini = min(c( size_ini[ inside],10))
  # ### Contraction Stage (iterT=1)
  if(m<=50){
    lam_seq =  seq(0.001, 0.0001, -0.0001) 
  }else{
    lam_seq =  seq(0.006, 0.0001, -0.0005)  
  }
  beta_pre = beta_ini; result_dist_lasso = dist_l1_est(c0,lam_seq,iterT=1,beta_pre)
  lambda_lasso1[inside] = result_dist_lasso$pen_lambda
  size_dist_lasso[inside] = sum(result_dist_lasso$coeff!=0)
  k_count = 1
  while (k_count < 10) {
    k_count = k_count + 1;
    result_old = result_dist_lasso
    beta_pre = result_old$coeff
    result_dist_lasso =  dist_l1_est(c0,lam_seq,iterT=1,beta_pre)
  }
  lambda_lasso5[inside] = result_dist_lasso$pen_lambda
  size_dist_lasso[inside] = sum(result_dist_lasso$coeff!=0)
  coeflasso_record[,inside] = result_dist_lasso$coeff
  ###########################  
  ## Tightening Stage (iterT>=2)
  beta_pre = coeflasso_record[,inside]; 
  result_dist_scad = dist_scad_est(c0,lam_seq,iterT=1,beta_pre )
  lambda_dist1[inside] = result_dist_scad$pen_lambda
  size_dist_scad1[inside] = sum(result_dist_scad$coeff!=0)
  t_count = 2;
  while (t_count <= 15) {
    t_count = t_count + 1;
    result_old = result_dist_scad
    beta_pre = result_old$coeff
    result_dist_scad =  dist_scad_est(c0,lam_seq,iterT=1,beta_pre)
  }
  size_dist_scad5[inside] = sum(result_dist_scad$coeff!=0)
  coef_record[,inside] = result_dist_scad$coeff 
  
  h = result_dist_scad$global_bandwidth
  ck =  (1.5*dnorm(qnorm(tau,sd=1),sd=1)^2/(2*qnorm(tau,sd=1)^2+1))^(1/3)  
  support = (result_dist_scad$coeff!=0);  support[test_v] = 1;   
  v_test = test_v
  for(j in 1:length(test_v)){
    v_test[j] = sum(support[1:test_v[j]])
  }
  sigma_hat = sigma_est(X_a_train[, support==1],y_train,result_dist_scad$coeff[support==1],k=ck*(size_dist_scad5[inside]/(n*m))^(1/3),s = sum(support),m=m)
  sigma_IM[,, inside] = sigma_hat$sigma2
  T_IM[ inside] = N*t(result_dist_scad$coeff[test_v]-0) %*% MASS::ginv(sigma_IM[,,inside]+0*diag(1,length(test_v))) %*% (result_dist_scad$coeff[test_v]-0)
  pvalue[inside] = pchisq(T_IM[inside], length(test_v), lower.tail = FALSE)
  
}

T_IM
pvalue

stopCluster(cl7)

rm("X","y","X_train","y_train","X_pred","y_pred", "X_a_train", "cl7","X_master","y_master", "beta_ini" 
   ,"result_dist_scad","result_dist_lasso" ,"beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove" )
save.image(file=paste0("Results/real_inference_n_",n,"_m_",m,'_tau_',tau,"_p_",p,".RData"))