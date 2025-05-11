AE_ini = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset))
AE_dist = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));

size_ini = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset))
size_dist = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));


lambda_conquer = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset))
lambda_conquer_scad = lambda_conquer
lambda_dist = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset))


dist_beta_record = array(rep(0,sim*(p)*length(Dset)),c((p),sim,length(Dset)))
ini_beta = array(rep(0,sim*(p)*length(Dset)),c((p),sim,length(Dset)))


T_local = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));
T_global = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));
T_IM = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset)); 
test_IM  = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset)); 
sigma_IM  = array( 0, c(length( test_v),length( test_v),sim,length(Dset)) ) ;
test_local = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));
test_global = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));


iter_num = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset))
iter_lasso = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset))
store_matrix_table = matrix(rep(0,5*length(Dset)), ncol = length(Dset) ) 


T_ora = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));
test_ora  = matrix(rep(0,sim*length(Dset)),nrow =sim,ncol = length(Dset));
sigma_ora  = array( 0, c(length( test_v),length( test_v),sim,length(Dset)) )