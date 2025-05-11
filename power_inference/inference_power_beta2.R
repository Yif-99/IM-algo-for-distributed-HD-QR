rm(list = ls())
###############################################
##### Empirical Power Analysis for beta_2=0 #####
test_v = 2; c_v = 0.0;
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
Dset=c(0,0.003,0.006,0.009,0.012,0.015) 
#######################################################
source("store_matrix_power.R")
source("aux_funcs_power.R")
#=======main=======

et1 = proc.time()[3]
for(ii in 1:sim){  
  set.seed(100+ii)  
  X = generateData(N,p,0.5,30)
  e = X[,p+1];  X = X[,1:p];   X_master = X[1:n,]
  clusterExport(cl7, list("X_master", "X",'ii' ) )
  inside = 0;
  
  results <- parLapply(cl7, 1:m, function(site){
    X_l <- X[(n*(site-1)+1):(n*site),];
    set.seed(120+ii)
    cv_result <- cv.glmnet(cbind(1,X_l[,-test_v]), X_l[,test_v]);
    w_result <- cv_result$glmnet.fit$beta[,cv_result$index[1]] / m
    return(list(w_result = w_result))
  })
  SDS_w_mean = rep(0,p);
  for (result in results) {
    SDS_w_mean <- SDS_w_mean + result$w_result
  }
  SDS_w_record[,ii] = SDS_w_mean
  
  for(delta in Dset){ 
    inside = inside + 1
    beta0 <- rep(0,p);  beta0[(idx)] = signal; beta0[test_v] = delta
    y <- X%*%beta0 + 0.7*e*X[,1];     y_master = y[1:n]
    beta0[1] <- beta0[1] + 0.7*qnorm(tau)
    lam_seq=seq(0.17, 0.04, -0.014)
    clusterExport(cl7, list("ii",'inside' ,"beta0","y","y_master","tau","lam_seq") )
    cat("\n##############################\n##############################\n m: "
        ,m,", log(c0):", log(c0),", (c0):",  (c0),", delta: ",delta,", sim: ", ii,"\n")
    ##########################################
    ### SDS
    SDS_beta_mean = rep(0,p+1);   
    results <- parLapply(cl7, 1:m, function(site){
      X_l <- X[(n*(site-1)+1):(n*site),];    y_l <- y[(n*(site-1)+1):(n*site)]
      set.seed(120+ii)
      beta_result <- conquer.cv.reg(X_l, y_l,tau=tau,lambdaSeq=lam_seq)$coeff.min / m   
      return(list(beta_result = beta_result))
    })
    for (result in results) {
      SDS_beta_mean <- SDS_beta_mean + result$beta_result
    }
    ###############SDS T=2,3,4#########################
    for(iter_SDS in 1:4){
      beta_test_SDS = nleqslv(SDS_beta_mean[test_v+1], func_2_solve)$x;  clusterExport(cl7, c("beta_test_SDS"))
      if(abs( SDS_beta_mean[test_v+1]- beta_test_SDS)<0.001){break;}
      SDS_beta_mean[test_v+1] = beta_test_SDS;  SDS_beta_mean[-test_v] = 0;
      results <- parLapply(cl7, 1:m, function(site){
        set.seed(120+ii)
        X_l <- X[(n*(site-1)+1):(n*site),]
        y_l <- y[(n*(site-1)+1):(n*site)]
        result <- conquer.cv.reg(X_l[,-test_v],y_l-X_l[,test_v]*beta_test_SDS,tau=tau,penalty="scad",lambdaSeq=lam_seq)$coeff.min 
        return(result)
      })
      SDS_beta_mean[-(test_v+1)] <- Reduce("+", results) / m
    }
    beta_test_SDS = nleqslv(SDS_beta_mean[test_v+1], func_2_solve)$x
    SDS_beta_mean[test_v+1] = beta_test_SDS;
    SDS_beta_record[,ii,inside]=SDS_beta_mean
    sigma = 0
    for(site in 1:m){
      X_l = X[(n*(site-1)+1):(n*site),]; y_l = y[(n*(site-1)+1):(n*site)]
      res_w = X_l[,test_v] - cbind(1,X_l[,-test_v])%*%SDS_w_mean;      res = y_l - cbind(1,X_l)%*%SDS_beta_mean
      band_k = 0.5 * ((1+log(n))/n)^0.5
      D_j = tau*(1-tau)*mean(res_w^2); Q_j = mean(dnorm(res/band_k)*res_w*X_l[,test_v])/band_k
      sigma = sigma + Q_j^(-2)*D_j/m
    }
    sigma_SDS[1,1,ii,inside] = sigma
    lower =  beta_test_SDS - qnorm(0.975)*sqrt(sigma/N); upper = beta_test_SDS + qnorm(0.975)*sqrt(sigma/N)
    T_SDS[ii,inside] = (c_v< lower || c_v>upper)
    ### Global test
    global_result = conquer(X[,idx],y,tau,tol = 1e-4,kernel="uniform",ci="asymptotic",alpha=0.05)
    global_result =  SQR_ora(cbind(1,X[,idx]), y, tau, global_result$coeff ,phi = 1/200, itermax=0, test_v+1,alpha=0.05,k=-1,eps=1e-8, verbose =0)
    if(length(test_v)==1){
      test_global[ii,inside] = 1*( c_v<global_result$lower[test_v] | c_v>global_result$upper[test_v]    )
    }else{
      T_global[ii,inside] = global_result$wald
      test_global[ii,inside] = (T_global[ii,inside] > qchisq(1-alpha, df = length(test_v))) * 1.0
    }
    local_result = conquer(X[1:n,idx],y[1:n],tau,kernel="uniform",ci="asymptotic",alpha=0.05)
    local_result =  SQR_ora(X[1:n,idx], y[1:n], tau, local_result$coeff[-1] ,phi = 1/200, itermax=0, test_v, alpha=0.05,k=-1,eps=1e-4, verbose =0)
    if(length(test_v)==1){
      test_local[ii,inside] = 1*( c_v<local_result$lower[test_v] | c_v>local_result$upper[test_v] )
    }else{
      T_local[ii,inside] = local_result$wald
      test_local[ii,inside] = (T_local[ii,inside] > qchisq(1-alpha, df = length(test_v))) * 1.0
    }
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
   
    ##########################################
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
    test_ora[ii,inside] = (T_ora[ii,inside] > qchisq(1-alpha, df = length(test_v))) * 1.0
    ########################################
  }

}
et2=proc.time()[3]

for(kkk in 1:length(Dset)){
  store_matrix_table[1, ] =  Dset
  store_matrix_table[2,((kkk-1)*1+1):(kkk*1)] =  c(mean(test_IM[1:sim,kkk]))
  store_matrix_table[3,((kkk-1)*1+1):(kkk*1)] =  c(mean(test_global[1:sim,kkk]))
  store_matrix_table[4,((kkk-1)*1+1):(kkk*1)] =  c(mean(test_local[1:sim,kkk]))
  store_matrix_table[5,((kkk-1)*1+1):(kkk*1)] =  c(mean(T_SDS[1:sim,kkk]))
}
write.csv(store_matrix_table,paste0("Results/power_beta2_n_",n,"_p_",p,"_m_",m,'_tau_',tau,".csv"),row.names=FALSE)

stopCluster(cl7)

rm("X","y","beta0","e","cl7","X_master","y_master","idx","beta_ini","beta_ora"
   ,"result_dist_scad","result_dist_lasso", "beta_pre","result_old")
all_objects <- ls();functions_to_remove <- sapply(all_objects, function(x) is.function(get(x)))
rm(list = all_objects[functions_to_remove]);rm("all_objects","functions_to_remove","vars_to_export")

save.image(file=paste0("Results/power_beta2_n_",n,"_p_",p,"_m_",m,'_tau_',tau,".RData"))

 
  