AE_ini = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
AE_dist_lasso = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
AE_dist = matrix(rep(0,sim*3),nrow =sim,ncol = 3);
AE_dist_ora = matrix(rep(0,sim*3),nrow =sim,ncol = 3);

size_ini = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
size_dist_lasso = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
size_dist = matrix(rep(0,sim*3),nrow =sim,ncol = 3);

lambda_conquer = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_conquer_scad = lambda_conquer
lambda_lasso1 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_lasso5 = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
lambda_dist = matrix(rep(0,sim*3),nrow =sim,ncol = 3)


dist_beta_record_1 = array(rep(0,sim*(p)*3),c((p),sim,3))
dist_beta_record = array(rep(0,sim*(p)*3),c((p),sim,3))
dist_beta_ora_record = array(rep(0,sim*(p)*3),c((p),sim,3))
ini_beta = array(rep(0,sim*(p)*3),c((p),sim,3))


T_IM = matrix(rep(0,sim*3),nrow =sim,ncol = 3); 
T_ora = matrix(rep(0,sim*3),nrow =sim,ncol = 3); 
test_IM  = matrix(rep(0,sim*3),nrow =sim,ncol = 3); 
test_ora  = matrix(rep(0,sim*3),nrow =sim,ncol = 3); 

sigma_IM  = array( 0, c(length( test_v),length( test_v),sim,3) ) ; 
sigma_ora  = array( 0, c(length( test_v),length( test_v),sim,3) ); 

iter_num = matrix(rep(0,sim*3),nrow =sim,ncol = 3)
store_matrix_table = matrix(rep(0,4*3),nrow = 4 ) 

