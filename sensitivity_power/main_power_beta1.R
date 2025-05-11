rm(list = ls())
###############################################
##### Empirical Power Analysis for beta_1=1.5 #####
test_v = 1; c_v = 1.5;
signal = 1.5;  idx = 1:(5); alpha = 0.05
##################################
sim = 500
######################
p = 2000
N = 2000; m = 10; 
n = N/m; 
tau = 0.5
if(tau == 0.3){Tmax=8}else if(tau==0.5){ Tmax=8}else{ Tmax=8}
######################
Dset=c(0.5,0.51,0.52,0.53,0.54,0.55)
#######################################################
source("store_matrix_power.R")
source("aux_funcs_power.R")
#=======main=======

ck =  0.52 # 0.52,0.57,0.62,0.67,0.72

et1 = proc.time()[3]
for(ii in 1:sim){  
  set.seed(100+ii)  
  X = generateData(N,p,rho,B)
  e = X[,p+1];  X = X[,1:p];  X_master = X[1:n,]
  beta0 <- rep(0,p); beta0[idx] <- signal; beta0[test_v] <- signal
  ora_support = (beta0 != 0);  ora_support[test_v] = 1
  y <- X%*%beta0 + 0.7*e*X[,1];     y_master = y[1:n]
  clusterExport(cl7, list("X_master", "X" ) )
  inside = 0;
  for(tau in Dset){ 
    inside = inside + 1
    beta0 <- rep(0,p); beta0[idx] <- signal;  beta0[1] <- beta0[1] + 0.7*qnorm(tau)
    clusterExport(cl7, list("ii",'inside' ,"beta0","y","y_master","tau") )
    cat("\n##############################\n##############################\n p: "
        ,p,", tau: ",tau,", sim: ", ii,"\n")
    ##########################################
    ### Global test
    global_result = conquer(X[,idx],y,tau,kernel="Gaussian",ci="both",alpha=0.05)
    global_result =  SQR_ora(X[,idx], y, tau, global_result$coeff[-1] ,phi = 1/200, itermax=1000, alpha=0.05,k=-1,eps=1e-4, verbose = FALSE)
    T_global0[ii,inside] = 1*( c_v<global_result$lower[test_v] | c_v>global_result$upper[test_v]    )
    local_result = conquer(X[1:n,idx],y[1:n],tau,kernel="Gaussian",ci="both",alpha=0.05)
    local_result =  SQR_ora(X[1:n,idx], y[1:n], tau, local_result$coeff[-1] ,phi = 1/200, itermax=1000, alpha=0.05,k=-1,eps=1e-4, verbose = FALSE)
    T_local0[ii,inside] = 1*( c_v<local_result$lower[test_v] | c_v>local_result$upper[test_v] )
    ##########################################
    ### Local initialization
    beta_ini = ini_est(seq(0.5, 0.07, -0.01))
    ini_beta[, ii, inside] = beta_ini;  hat_s_ini = sum(beta_ini != 0)
    size_ini[ii, inside] = hat_s_ini;    AE_ini[ii, inside] = sum(abs(beta_ini - beta0))
    ##########################################
    ### Contraction Stage (T=1)
    beta_pre = ini_beta[,ii,inside];
    lam_seq =  exp(seq(-5, -.8, 0.065))
    k_count = 0;
    while (k_count < 6 ) {
      k_count = k_count + 1;
      if(k_count > 1){
        result_old = result_dist_lasso
        beta_pre = result_old$coeff
      }
      result_dist_lasso =  dist_l1_est(c0,lam_seq,T=1,beta_pre)
    }
    # ##########################################
    ## Tightening Stage (T>=2)
    hat_s_ini = min(c( size_ini[ii,inside],15))
    lam_seq =  exp(seq(-4, -0.6, 0.05))
    t_count = 1; ifora=FALSE;
    while ( t_count < Tmax) {
      t_count = t_count + 1;
      if(t_count > 2){
        beta_pre = result_dist_scad$coeff;
      }
      if(t_count == Tmax){
        ifora = 1
      }
      result_dist_scad =  dist_scad_est(c0,lam_seq,T=1,beta_pre)
    }
    
    lambda_dist[ii,inside] = result_dist_scad$pen_lambda
    size_dist[ii,inside] = sum(result_dist_scad$coeff!=0)
    AE_dist[ii,inside] = sum( abs(result_dist_scad$coeff[]-beta0[]) )
    dist_beta_record[,ii,inside] = result_dist_scad$coeff
    iter_num[ii,inside] = t_count
    b = result_dist_scad$local_bandwidth
    h = result_dist_scad$global_bandwidth
    support = (result_dist_scad$coeff!=0);  support[test_v] = 1;
    sigma2 = 0
    for(site in 1:m){
      sigma2 = sigma2 + sigmahat(X[(n*(site-1)+1):(n*site),(support==1)],y[(n*(site-1)+1):(n*site)],result_dist_scad$coeff[(support==1)],
                                 test_v, tau,  k = ck* ((5+log(N))/N)^(1/3)  ,h=h)
    }
    sigma2 = sigma2/m
    T_IM[ii,inside] = N*t(result_dist_scad$coeff[test_v]-c_v) %*% solve(sigma2)  %*% (result_dist_scad$coeff[test_v]-c_v)
    sigma_IM[,,ii,inside] = sigma2
    T_ora[ii,inside] = N*t(result_dist_scad$coeff_ora[test_v]-c_v) %*% solve(sigma2)  %*% (result_dist_scad$coeff_ora[test_v]-c_v)
    sigma_ora[,,ii,inside] = sigma2
    # Test
    test_IM[ii,inside] = (T_IM[ii,inside] > qchisq(1-alpha, df = length(test_v))) * 1.0
    test_ora[ii,inside] = (T_ora[ii,inside] > qchisq(1-alpha, df = length(test_v))) * 1.0
    ########################################
  }
 
}

for(kkk in 1:length(Dset)){
  store_matrix_table[1, ] =  Dset
  store_matrix_table[2,((kkk-1)*1+1):(kkk*1)] =  c(mean(test_IM[1:sim,kkk]))
  store_matrix_table[3,((kkk-1)*1+1):(kkk*1)] =  c(mean(T_global0[1:sim,kkk]))
  store_matrix_table[4,((kkk-1)*1+1):(kkk*1)] =  c(mean(T_local0[1:sim,kkk]))
}
et2=proc.time()[3]
write.csv(store_matrix_table,paste0("Results/sen_power_beta1_n_",n,"_p_",p,"_m_",m,'_ck_',ck,".csv"),row.names=FALSE)
 

stopCluster(cl7)

rm("X","y","beta0","e","cl7","X_master","y_master","idx","beta_ini","beta_ora"
   ,"result_dist_scad","result_dist_lasso", "beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove","vars_to_export")

save.image(file=paste0("Results/sen_power_beta1_n_",n,"_p_",p,"_m_",m,'_ck_',ck,".RData"))
