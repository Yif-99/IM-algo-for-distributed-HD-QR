size_global_lasso = matrix(rep(0,sim*5),nrow =sim,ncol = 5);size_global_scad = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
AE_global_lasso = matrix(rep(0,sim*5),nrow =sim,ncol = 5);AE_global_scad = matrix(rep(0,sim*5),nrow =sim,ncol = 5)

AE_ini = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
AE_dist_lasso = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
AE_dist = matrix(rep(0,sim*5),nrow =sim,ncol = 5);

size_ini = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
size_dist_lasso = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
size_dist = matrix(rep(0,sim*5),nrow =sim,ncol = 5);

lambda_conquer = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
lambda_conquer_scad = lambda_conquer
lambda_lasso1 = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
lambda_lasso5 = matrix(rep(0,sim*5),nrow =sim,ncol = 5)
lambda_dist = matrix(rep(0,sim*5),nrow =sim,ncol = 5)

dist_beta_record_1 = array(rep(0,sim*p*5),c(p,sim,5))
dist_beta_record = array(rep(0,sim*p*5),c(p,sim,5))
ini_beta = array(rep(0,sim*p*5),c(p,sim,5))


store_matrix_table = matrix(rep(0,15*10),nrow = 15,ncol = 10) 
store_sd_table = store_matrix_table

