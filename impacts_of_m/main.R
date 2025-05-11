rm(list=ls())
sim =  100
#################
p = 2000 # choice 1
tau = 0.7 # choice 2
n = 200
truecoef = c(1.5,0,2,0,0,2,0,2,0,2,0,2,0,0,2,0,2,rep(0,p-17))
idx = (truecoef!=0); beta0 <- truecoef
source("store_matrix.R")
source("aux_funcs.R")
#=======main=======
inside = 0   
for(m in c(30,50,70,100,200) ){ 
  if(m==30){inside=1}else if(m==50){inside=2}else if(m==70){inside=3}else if(m==100){inside=4}else{inside=5}
  N = n*m  
  clusterExport(cl7, list("m", "N"))
  et1 = proc.time()[3]
  for(ii in 1:sim){
    set.seed(100+50*ii)
    X = generateData(N,p,0.8,30)
    e = X[,p+1]; X = X[,1:p]; X_master = X[1:n,]
    y <- X%*%beta0 + (e - qnorm(tau))*(0.2*X[,1] + 1);   y_master = y[1:n]
    clusterExport(cl7, list("X_master","y_master", "X","y","tau", "ii",'inside','eta') )
    cat("##############################\n##############################\n m:"
          ,m,", tau: ",tau,", sim: ", ii,"\n")
    ### Local initialization
    beta_pre =  beta0
    beta_ini = ini_est(seq(0.5, 0.07, -0.01))
    ini_beta[,ii,inside] = beta_ini; hat_s_ini = sum(beta_ini!=0)
    size_ini[ii,inside] = hat_s_ini; AE_ini[ii,inside] = sum( abs(beta_ini-beta0) )
    hat_s_ini = min(c(hat_s_ini,15))
    ## Global
    beta_ini = ini_beta[,ii,inside];  beta_pre = beta_ini;
    lam_seq =  exp(seq(-5, -2, 0.4))
    result_global  = global_est(lam_seq,beta_pre,verbose=FALSE)
    lambda_conquer[ii,inside] = result_global$pen_lambda
    size_global_lasso[ii,inside] = sum(result_global$coeff !=0)
    AE_global_lasso[ii,inside] = sum( abs(result_global$coeff -beta0) )
    beta_pre = result_global$coeff
    result_global = global_est(lam_seq,beta_pre,pen="SCAD",verbose=FALSE)
    lambda_conquer_scad[ii,inside] = result_global$pen_lambda
    size_global_scad[ii,inside] = sum(result_global$coeff !=0)
    AE_global_scad[ii,inside] = sum( abs(result_global$coeff -beta0) )
    ## Contraction Stage (T=1)
    beta_pre = ini_beta[,ii,inside];
    lam_seq = exp(seq(-5, -1.5, 0.072))
    k_count = 0;
    while (k_count < 6  ) {
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
    ## Tightening Stage (T>=2)
    hat_s_ini = min(c( size_ini[ii,inside],15))
    lam_seq =  exp(seq(-5.5, -1.5, 0.07))
    beta_pre = dist_beta_record_1[,ii,inside]
    t_count = 1;
    while ( t_count < 15 ) {
      t_count = t_count + 1
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

  for(kkk in 1:5){
    store_matrix_table[1,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist_lasso[1:sim,kkk]),mean(size_dist_lasso[1:sim,kkk]))
    store_matrix_table[3,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_dist[1:sim,kkk]),mean(size_dist[1:sim,kkk]))
    store_matrix_table[15,((kkk-1)*2+1):(kkk*2)] =  c(mean(AE_global_scad[1:sim,kkk]),mean(size_global_scad[1:sim,kkk]))
  }
  for(kkk in 1:5){
    store_sd_table[1,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist_lasso[1:sim,kkk]),sd(size_dist_lasso[1:sim,kkk]))
    store_sd_table[3,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_dist[1:sim,kkk]),sd(size_dist[1:sim,kkk]))
    store_sd_table[15,((kkk-1)*2+1):(kkk*2)] =  c(sd(AE_global_scad[1:sim,kkk]),sd(size_global_scad[1:sim,kkk]))
  }
  store_sd_table = store_sd_table / sqrt(sim)
  et2=proc.time()[3]
  write.csv(store_sd_table,paste0("Results/sd_effects_on_m_tau_",tau,"_m_",m,"_p",p,".csv"),row.names=FALSE)
  write.csv(store_matrix_table,paste0("Results/effects_on_m_tau_",tau,"_m_",m,"_p",p,".csv"),row.names=FALSE)

}

stopCluster(cl7)

rm("X","y","beta0","e",'result_global' ,"cl7","X_master","y_master","idx","beta_ini","beta_ora"
   ,"result_dist_scad","result_dist_lasso","truecoef","beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove","vars_to_export")
save.image(file=paste0("Results/effects_on_m_tau_",tau,"_m_",m,"_p",p,".RData"))


