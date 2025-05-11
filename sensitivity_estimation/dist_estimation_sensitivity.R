rm(list=ls())
############ check 1 #############
sim =  100
############ check 2 #############
p = 2000; # 2000 5000
tau = 0.3 # 0.3 0.5 0.7
N = 2000; m = 10; n = N/m;
truecoef = c(1.5,0,2,0,0,2,0,2,0,2,0,2,0,0,2,0,2,rep(0,p-17))
idx = (truecoef!=0); beta0 <- truecoef
source("store_matrix_sen.R")
source("aux_funcs_sen.R")

#=======main=======
et1 = proc.time()[3]
for(ii in  1:sim){  
  set.seed(100+50*ii)   
  X = generateData(N,p,0.8,30)
  e = X[,p+1]; X = X[,1:p]; X_master = X[1:n,]
  y <- X%*%beta0 + (e - qnorm(tau))*(0.2*X[,1] + 1);   y_master = y[1:n]
  clusterExport(cl7, list("tau","y","X","y_master","X_master"))
  inside = 0  
  ### Local initialization
  beta_pre =  beta0
  beta_ini = ini_est(seq(0.3, 0.07, -0.01))
  ini_beta[,ii,1] = beta_ini; hat_s_ini = sum(beta_ini!=0)
  size_ini[ii,1] = hat_s_ini; AE_ini[ii,1] = sum( abs(beta_ini-beta0) )
  hat_s_ini = min(c(hat_s_ini,15))
  ##########################################
  for( c0 in exp(c(-3,-2,-1,-0.5,0,0.5,1))  ){  
    inside = inside + 1
    clusterExport(cl7, list("ii",'inside','c0') )
    cat("##############################\n##############################\n log c0:"
          ,log(c0),", n: ",n,", m: ",m,", p: ",p,", sim: ", ii,"\n")
    ##########################################
    ### Contraction Stage (T=1)
    beta_pre = ini_beta[,ii,1];
    lam_seq =  exp(seq(-5, -.8, 0.065))
    k_count = 0;
    while (k_count < 10  ) {
      k_count = k_count + 1;
      if(k_count > 1){
        result_old = result_dist_lasso
        beta_pre = result_old$coeff
      }
      result_dist_lasso =  dist_l1_est(c0,lam_seq,T=1,beta_pre)
    }
    dist_beta_record_1[,ii,inside] = result_dist_lasso$coeff
    size_dist_lasso[ii,inside] = sum(result_dist_lasso$coeff!=0)
    AE_dist_lasso[ii,inside] = sum( abs(result_dist_lasso$coeff-beta0) )
    ### Tightening Stage (T>=2)
    lam_seq =   exp(seq(-4.5, -1.4, 0.08))
    beta_pre = dist_beta_record_1[,ii,inside]
    t_count = 1;
    while ( t_count < 15  ) {
      t_count = t_count + 1;
      if(t_count > 2){
        result_old = result_dist_scad
        beta_pre = result_old$coeff
      }
      result_dist_scad =  dist_scad_est(c0,lam_seq,T=1,beta_pre)
    }
    lambda_dist[ii,inside] = result_dist_scad$pen_lambda
    size_dist[ii,inside] = sum(result_dist_scad$coeff!=0)
    AE_dist[ii,inside] = sum( abs(result_dist_scad$coeff-beta0) )
    dist_beta_record[,ii,inside] = result_dist_scad$coeff
  }
}
et2=proc.time()[3]


for(kkk in 1:7){
    store_matrix_table[1,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_lasso[1:sim,kkk]),mean(size_dist_lasso[1:sim,kkk]))
    store_matrix_table[3,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist[1:sim,kkk]),mean(size_dist[1:sim,kkk]))
}
store_sd_table = store_matrix_table
for(kkk in 1:7){
  store_sd_table[1,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_lasso[1:sim,kkk]),sd(size_dist_lasso[1:sim,kkk]))
  store_sd_table[3,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist[1:sim,kkk]),sd(size_dist[1:sim,kkk]))
}
store_sd_table = store_sd_table / sqrt(sim)
store_matrix_table[,c(1,3,5,7,9,11,13)]

write.csv(store_sd_table,paste0("Results/sd_sensitivity_estimation_n_",n,"_p_",p,"_tau_",tau,".csv"),row.names=FALSE)
write.csv(store_matrix_table,paste0("Results/sensitivity_estimation_n_",n,"_p_",p,"_tau_",tau,".csv"),row.names=FALSE)


stopCluster(cl7)

rm("X","y","beta0","e","cl7","X_master","y_master","idx","beta_ini", 
   "result_dist_lasso","result_dist_scad","truecoef","beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove","vars_to_export")

save.image(file=paste0("Results/sensitivity_estimation_n_",n,"_p_",p,"_tau_",tau,".RData"))


