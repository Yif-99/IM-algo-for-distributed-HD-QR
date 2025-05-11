size_global_lasso = matrix(rep(0,sim*7),nrow =sim,ncol = 7);size_global_scad = matrix(rep(0,sim*7),nrow =sim,ncol = 7)
AE_global_lasso = matrix(rep(0,sim*7),nrow =sim,ncol = 7);AE_global_scad = matrix(rep(0,sim*7),nrow =sim,ncol = 7)

AE_ini = matrix(rep(0,sim*1),nrow =sim,ncol = 1)
AE_dist_lasso = matrix(rep(0,sim*7),nrow =sim,ncol = 7)
AE_dist = matrix(rep(0,sim*7),nrow =sim,ncol = 7);

size_ini = matrix(rep(0,sim*1),nrow =sim,ncol = 1)
size_dist_lasso = matrix(rep(0,sim*7),nrow =sim,ncol = 7)
size_dist = matrix(rep(0,sim*7),nrow =sim,ncol = 7);

lambda_conquer = matrix(rep(0,sim*7),nrow =sim,ncol = 7)
lambda_conquer_scad = lambda_conquer
lambda_lasso1 = matrix(rep(0,sim*7),nrow =sim,ncol = 7)
lambda_lasso5 = matrix(rep(0,sim*7),nrow =sim,ncol = 7)
lambda_dist = matrix(rep(0,sim*7),nrow =sim,ncol = 7)

dist_beta_record_1 = array(rep(0,sim*p*7),c(p,sim,7))
dist_beta_record = array(rep(0,sim*p*7),c(p,sim,7))
ini_beta = array(rep(0,sim*p*1),c(p,sim,1))


store_matrix_table = matrix(rep(0,15*14),nrow = 15,ncol = 14) #AE SIZE*7
store_sd_table = store_matrix_table

