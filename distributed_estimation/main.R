rm(list=ls())
sim =  100
#################
p = 200 
n = 1000; m = 100; N = n*m 
truecoef = c(1.5,0,2,0,0,2,0,2,0,2,0,2,0,0,2,0,2,rep(0,p-17))
idx = (truecoef!=0); beta0 <- truecoef
#################
source("store_matrix.R")
source("aux_funcs.R")
#=======main=======
et1 = proc.time()[3]
for(ii in 1:sim){
  set.seed(100+50*ii)
  X = generateData(N,p,0.9,30)
  e = X[,p+1]; X = X[,1:p]; X_master = X[1:n,]
  clusterExport(cl7, list("X", "X_master")) 
  inside = 0   
  for( tau in c(0.3,0.5,0.7) ){
    inside = inside + 1
    y <- X%*%beta0 + (e - qnorm(tau))*(0.2*X[,1] + 1);   y_master = y[1:n]
    clusterExport(cl7, list("y_master","y","tau", "ii",'inside','c0') )
    cat("##############################\n##############################\n m:"
          ,m,', p: ',p,', n: ',n,", tau: ",tau,", sim: ", ii,"\n")
    ##########################################
    ### Local initialization
    beta_ini = ini_est(seq(0.3, 0.07, -0.01))
    ini_beta[,ii,inside] = beta_ini; hat_s_ini = sum(beta_ini!=0)
    size_ini[ii,inside] = hat_s_ini; AE_ini[ii,inside] = sum( abs(beta_ini-beta0) )
    hat_s_ini = min(c(hat_s_ini,15))
    ### Global L1
    beta_ini = ini_beta[,ii,inside];  beta_pre = beta_ini;
    lam_seq = seq(0.15, 0.01, -0.015)
    result_global  = global_est(lam_seq,beta_pre,verbose=FALSE)
    lambda_conquer[ii,inside] = result_global$pen_lambda
    size_global_lasso[ii,inside] = sum(result_global$coeff !=0)
    AE_global_lasso[ii,inside] = sum( abs(result_global$coeff -beta0) )
    ### Global SCAD
    beta_pre = result_global$coeff
    lam_seq = seq(lambda_conquer[ii,inside]+0.06, max(0.01,lambda_conquer[ii,inside]-0.02), -0.01)
    result_global = global_est(lam_seq,beta_pre,pen="SCAD",verbose=FALSE)
    lambda_conquer_scad[ii,inside] = result_global$pen_lambda
    size_global_scad[ii,inside] = sum(result_global$coeff !=0)
    AE_global_scad[ii,inside] = sum( abs(result_global$coeff -beta0) )
    ##########################################
    ### Contraction Stage (T=1)
    beta_pre = beta_ini;   lam_seq = seq(0.15, 0.055, -0.003)
    result_dist_lasso = dist_l1_est(c0,lam_seq,T=1,beta_pre)
    lambda_lasso1[ii,inside] = result_dist_lasso$pen_lambda
    lam_seq = seq(0.082, 0.037, -0.003);
    k_count = 1;
    while (k_count < 8 ) {
      k_count = k_count + 1;
      result_old = result_dist_lasso
      beta_pre = result_old$coeff
      result_dist_lasso =  dist_l1_est(c0,lam_seq,T=1,beta_pre)
    }
    dist_beta_record_1[,ii,inside] = result_dist_lasso$coeff
    size_dist_lasso[ii,inside] = sum(result_dist_lasso$coeff!=0)
    AE_dist_lasso[ii,inside] = sum( abs(result_dist_lasso$coeff-beta0) )
    ## Tightening Stage (T>=2)
    hat_s_ini = min(c( size_ini[ii,inside],15))
    if(p<=2000){
      lam_seq =  exp(seq(-4,-0.6,0.05))
    }else{
      lam_seq =  exp(seq(-3.5,-0.6,0.05))
    }
    beta_pre =  dist_beta_record_1[,ii,inside]; 
    result_dist_scad = dist_scad_est(c0,lam_seq,T=1,beta_pre )
    lambda_dist2[ii,inside] = result_dist_scad$pen_lambda
    size_dist_2[ii,inside] = sum(result_dist_scad$coeff!=0)
    AE_dist_2[ii,inside] = sum( abs(result_dist_scad$coeff-beta0) )
    size_dist_ora2[ii,inside] = sum(result_dist_scad$coeff_ora!=0)
    AE_dist_ora2[ii,inside] = sum( abs(result_dist_scad$coeff_ora-beta0[(beta0 != 0)]) )
    beta_ora = rep(0,p);    beta_ora[beta0 != 0]=result_dist_scad$coeff_ora
    dist_beta_record_2[,ii,inside] = result_dist_scad$coeff
    dist_beta_ora_record_2[,ii,inside] = result_dist_scad$coeff_ora
    if(p<=2000){
      lam_seq = exp(seq(-4,-0.6,0.05))
    }else{
      lam_seq = exp(seq(-3.5,-0.6,0.05)) 
    }
    result_old = result_dist_scad;
    t_count = 2;
    while (  t_count < 6  ) {
      t_count = t_count + 1;
      beta_pre = result_dist_scad$coeff;  
      result_dist_scad =  dist_scad_est(c0,lam_seq,T=1,beta_pre)
    }
    result_old = result_dist_scad
    lambda_dist6[ii,inside] = result_dist_scad$pen_lambda
    size_dist_6[ii,inside] = sum(result_dist_scad$coeff!=0)
    AE_dist_6[ii,inside] = sum( abs(result_dist_scad$coeff-beta0) )
    size_dist_ora6[ii,inside] = sum(result_dist_scad$coeff_ora!=0)
    AE_dist_ora6[ii,inside] = sum( abs(result_dist_scad$coeff_ora-beta0[(beta0 != 0)]) )
    beta_ora = rep(0,p);      beta_ora[beta0 != 0]=result_dist_scad$coeff_ora
    dist_beta_record_6[,ii,inside] = result_dist_scad$coeff
    dist_beta_ora_record_6[,ii,inside] = result_dist_scad$coeff_ora
    result_dist_scad$coeff =  dist_beta_record_6[,ii,inside]
    t_count = 6;
    result_old = result_dist_scad;
    while (  t_count < 10  ) {
      t_count = t_count + 1;
      beta_pre = result_dist_scad$coeff;  
      result_dist_scad =  dist_scad_est(c0,lam_seq,T=1,beta_pre)
    }
    lambda_dist10[ii,inside] = result_dist_scad$pen_lambda
    size_dist_10[ii,inside] = sum(result_dist_scad$coeff!=0)
    AE_dist_10[ii,inside] = sum( abs(result_dist_scad$coeff-beta0) )
    size_dist_ora10[ii,inside] = sum(result_dist_scad$coeff_ora!=0)
    AE_dist_ora10[ii,inside] = sum( abs(result_dist_scad$coeff_ora-beta0[(beta0 != 0)]) )
    beta_ora = rep(0,p);  beta_ora[beta0 != 0]=result_dist_scad$coeff_ora
    dist_beta_record_10[,ii,inside] = result_dist_scad$coeff
    dist_beta_ora_record_10[,ii,inside] = result_dist_scad$coeff_ora
  }
}
for(kkk in 1:3){
  store_matrix_table[1,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_lasso[1:sim,kkk]),mean(size_dist_lasso[1:sim,kkk]))
  store_matrix_table[2,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_2[1:sim,kkk]),mean(size_dist_2[1:sim,kkk]))
  store_matrix_table[3,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_6[1:sim,kkk]),mean(size_dist_6[1:sim,kkk]))
  store_matrix_table[4,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_10[1:sim,kkk]),mean(size_dist_10[1:sim,kkk]))
  store_matrix_table[5,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_ora2[1:sim,kkk]),mean(size_dist_ora2[1:sim,kkk]))
  store_matrix_table[6,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_ora6[1:sim,kkk]),mean(size_dist_ora6[1:sim,kkk]))
  store_matrix_table[7,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_ora10[1:sim,kkk]),mean(size_dist_ora10[1:sim,kkk]))
  store_matrix_table[14,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_global_lasso[1:sim,kkk]),mean(size_global_lasso[1:sim,kkk]))
  store_matrix_table[15,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_global_scad[1:sim,kkk]),mean(size_global_scad[1:sim,kkk]))
}
for(kkk in 1:3){
  store_sd_table[1,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_lasso[1:sim,kkk]),sd(size_dist_lasso[1:sim,kkk]))
  store_sd_table[2,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_2[1:sim,kkk]),sd(size_dist_2[1:sim,kkk]))
  store_sd_table[3,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_6[1:sim,kkk]),sd(size_dist_6[1:sim,kkk]))
  store_sd_table[4,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_10[1:sim,kkk]),sd(size_dist_10[1:sim,kkk]))
  store_sd_table[5,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_ora2[1:sim,kkk]),sd(size_dist_ora2[1:sim,kkk]))
  store_sd_table[6,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_ora6[1:sim,kkk]),sd(size_dist_ora6[1:sim,kkk]))
  store_sd_table[7,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_ora10[1:sim,kkk]),sd(size_dist_ora10[1:sim,kkk]))
  store_sd_table[14,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_global_lasso[1:sim,kkk]),sd(size_global_lasso[1:sim,kkk]))
  store_sd_table[15,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_global_scad[1:sim,kkk]),sd(size_global_scad[1:sim,kkk]))
}
store_sd_table = store_sd_table / sqrt(sim)
et2=proc.time()[3]


write.csv(store_sd_table,paste0("Results/sd_dist_est_n_",n,"_m_",m,"_p_",p,".csv"),row.names=FALSE)
write.csv(store_matrix_table,paste0("Results/dist_est_n_",n,"_m_",m,"_p_",p,".csv"),row.names=FALSE)

stopCluster(cl7)

rm("X","y","beta0","e",'result_global' ,"cl7","X_master","y_master","idx","beta_ini","beta_ora"
   ,"result_dist_scad","result_dist_lasso","truecoef","beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove","vars_to_export")
save.image(file=paste0("Results/dist_est_n_",n,"_m_",m,"_p_",p,".RData"))


