size_global_lasso = matrix(rep(0,sim*3),nrow =sim,ncol = 3);size_global_scad = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
AE_global_lasso = matrix(rep(0,sim*3),nrow =sim,ncol = 3);AE_global_scad = matrix(rep(0,sim*3),nrow =sim,ncol = 3)

AE_ini = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
AE_dist_lasso = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
AE_dist_2 = matrix(rep(0,sim*3),nrow =sim,ncol = 3); AE_dist_6 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);AE_dist_10 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
AE_dist_ora2 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);AE_dist_ora6 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);AE_dist_ora10 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)

size_ini = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
size_dist_lasso = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
size_dist_2 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);size_dist_6 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);
size_dist_10 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
size_dist_ora2 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);size_dist_ora6 = matrix(rep(0,sim*3),nrow =sim,ncol = 3);
size_dist_ora10 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)

lambda_conquer = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_conquer_scad = lambda_conquer
lambda_lasso1 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_lasso5 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_dist2 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_dist6 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_dist10 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)

dist_beta_record_1 = array(rep(0,sim*p*3),c(p,sim,3))
dist_beta_record_2 = array(rep(0,sim*p*3),c(p,sim,3))
dist_beta_ora_record_2 = array(rep(0,sim*p*3),c(p,sim,3))
dist_beta_record_6 = array(rep(0,sim*p*3),c(p,sim,3))
dist_beta_ora_record_6 = array(rep(0,sim*p*3),c(p,sim,3))
dist_beta_record_10 = array(rep(0,sim*p*3),c(p,sim,3))
dist_beta_ora_record_10 = array(rep(0,sim*p*3),c(p,sim,3))
ini_beta = array(rep(0,sim*p*3),c(p,sim,3))


store_matrix_table = matrix(rep(0,15*6),nrow = 15,ncol = 6) #AE SIZE*3
store_sd_table = store_matrix_table

