size_dist_lasso = matrix(rep(0,n_cols), ncol = n_cols)
size_dist_scad1 = matrix(rep(0,n_cols), ncol = n_cols);
size_dist_scad5 = matrix(rep(0,n_cols), ncol = n_cols);
size_ini= matrix(rep(0,n_cols), ncol = n_cols)

lambda_ini = matrix(rep(0,n_cols), ncol = n_cols)
lambda_lasso1 = matrix(rep(0,n_cols), ncol = n_cols)
lambda_lasso5 = matrix(rep(0,n_cols), ncol = n_cols)
lambda_dist1 = matrix(rep(0,n_cols), ncol = n_cols)

check_pred_ini = matrix(rep(0,n_cols), ncol = n_cols) 
check_pred_scad = matrix(rep(0,n_cols), ncol = n_cols) 
check_pred_lasso = matrix(rep(0,n_cols), ncol = n_cols) 
check_pred_lasso1 = matrix(rep(0,n_cols), ncol = n_cols) 

coef_ini_record = matrix(rep(0,(p+1)*n_cols), ncol = n_cols)
coef_record = matrix(rep(0,(p+1)*n_cols), ncol = n_cols)
coeflasso_record = matrix(rep(0,(p+1)*n_cols), ncol = n_cols)

if(ifglobal){
  size_global_scad = matrix(rep(0,n_cols), ncol = n_cols)
  size_global_lasso = matrix(rep(0,n_cols), ncol = n_cols)
  coef_global_record = matrix(rep(0,(p+1)*n_cols), ncol = n_cols)
  coef_global_record_lasso = matrix(rep(0,(p+1)*n_cols), ncol = n_cols)
  check_pred_global = matrix(rep(0,n_cols), ncol = n_cols) 
  check_pred_global_lasso = matrix(rep(0,n_cols), ncol = n_cols) 
  lambda_global_scad = matrix(rep(0,n_cols), ncol = n_cols)
  lambda_global_lasso = matrix(rep(0,n_cols), ncol = n_cols)
}
