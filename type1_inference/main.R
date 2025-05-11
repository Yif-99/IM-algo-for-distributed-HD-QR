rm(list = ls())
##### Empirical Type I Error for beta_2=0 #####
test_v = c(2); c_v = 0.0; true_signal = 0.0
signal = 1.5;  idx = 1:(5); alpha = 0.05
###############################################
sim = 500
######################
p = 2000
N = 2000; m = 10; 
n = N/m; 
####################################################
source("store_matrix.R")
source("aux_funcs.R")
#=======main=======
et1 = proc.time()[3]
for(ii in 1:sim){
  set.seed(100+ii^2)
  X = generateData(N ,p,0.5,30)
  e = X[,p+1];  X = X[,1:p];  X_master = X[1:n,]
  beta0 <- rep(0,p);  beta0[(idx)] = signal; beta0[test_v] = true_signal
  y <- X%*%beta0 + 0.7*e*X[,1];     y_master = y[1:n]
  clusterExport(cl7, list("X_master","y_master", "X","y" ) )
  for( tau in c(0.3,0.5,0.7) ){ 
    if(tau == 0.3){inside=1;Tmax=12}else if(tau==0.5){inside=2;Tmax=10}else{inside=3;Tmax=10}
    
    beta0 <- rep(0,p);  beta0[(idx)] = signal; beta0[test_v] = true_signal
    beta0[1] <- beta0[1] + 0.7*qnorm(tau)
    clusterExport(cl7, list('c0','tau', "ii",'inside' ,"beta0") )
    cat("\n##############################\n##############################\n m: "
        ,m,", tau: ",tau,", sim: ", ii,"\n")
    ##########################################
    ### Local initialization
    beta_pre = beta0
    beta_ini = ini_est(seq(0.5, 0.07, -0.01))
    ini_beta[,ii,inside] = beta_ini; hat_s_ini = sum(beta_ini!=0)
    size_ini[ii,inside] = hat_s_ini; AE_ini[ii,inside] = sum( abs(beta_ini-beta0) )
    hat_s_ini = min(c(hat_s_ini,15))
    ##########################################
    ### Contraction Stage (T=1)
    beta_pre = ini_beta[,ii,inside];
    lam_seq =  exp(seq(-5, -.8, 0.065))
    k_count = 0;
    while (k_count < 8 ) {
      k_count = k_count + 1;
      if(k_count > 1){
        result_old = result_dist_lasso
        beta_pre = result_old$coeff
      }
      result_dist_lasso =  dist_l1_est(c0,lam_seq,T=1,beta_pre)
    }
    dist_beta_record_1[,ii,inside] = result_dist_lasso$coeff
    size_dist_lasso[ii,inside] = sum(result_dist_lasso$coeff!=0)
    AE_dist_lasso[ii,inside] = sum( abs(result_dist_lasso$coeff - beta0 ) )
    ### Tightening Stage (T>=2)
    beta_pre = lasso_beta[,ii,inside];
    hat_s_ini = 4
    lam_seq =  exp(seq(-4, -0.6, 0.05))
    t_count = 1; ifora=FALSE;
    while ( t_count < Tmax) {  
      t_count = t_count + 1;
      if(t_count>2){
        beta_pre = result_dist_scad$coeff;  
      }
      if(t_count == Tmax){
        ifora = TRUE
      }
      result_dist_scad =  dist_scad_est(c0,lam_seq,T=1,beta_pre)
    }
    lambda_dist[ii,inside] = result_dist_scad$pen_lambda
    size_dist[ii,inside] = sum(result_dist_scad$coeff!=0)
    AE_dist[ii,inside] = sum( abs(result_dist_scad$coeff[]-beta0[]) )
    support = (beta0 != 0);    support[test_v] = 1;    
    beta_ora = rep(0,p);  beta_ora[(support==1)]=result_dist_scad$coeff_ora
    AE_dist_ora[ii,inside] = sum( abs(beta_ora - beta0 ) )
    dist_beta_record[,ii,inside] = result_dist_scad$coeff
    dist_beta_ora_record[,ii,inside] = beta_ora
    iter_num[ii,inside] = t_count
    b = result_dist_scad$local_bandwidth
    h = result_dist_scad$global_bandwidth
    ck =  (1.5*dnorm(qnorm(tau,sd=1),sd=1)^2/(2*qnorm(tau,sd=1)^2+1))^(1/3)  
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
    test_ora[ ii,inside] = (T_ora[ii,inside] > qchisq(1-alpha, df = length(test_v))) * 1.0
    ########################################
  }

   
}
et2=proc.time()[3]

for(kkk in 1:3){
  store_matrix_table[1, ] =  c(0.3,0.5,0.7)
  store_matrix_table[2,((kkk-1)*1+1):(kkk*1)] =  c(mean(test_IM[1:sim  ,kkk]))
  store_matrix_table[3,((kkk-1)*1+1):(kkk*1)] =  c(mean(test_ora[1:sim   ,kkk]))
}
write.csv(store_matrix_table,paste0("Results/type1_n_",n,"_p_",p,"_m_",m,".csv"),row.names=FALSE)

stopCluster(cl7)

rm("X","y","beta0","e","cl7","X_master","y_master","idx","beta_ini","beta_ora"
   ,"result_dist_scad","result_dist_lasso", "beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove","vars_to_export")

save.image(file=paste0("Results/type1_n_",n,"_p_",p,"_m_",m,".RData"))
