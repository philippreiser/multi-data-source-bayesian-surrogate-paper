functions {


  vector gp_conditional_mean(
    array[] vector xw_new,     
    array[] vector xw_obs,     
    vector         f_obs,      
    real           alpha_gp,   
    array[] real   rho,        
    real           jitter
  ) {
    int N_obs = size(xw_obs);
    int N_new = size(xw_new);
    matrix[N_obs, N_obs] K = gp_exp_quad_cov(xw_obs, alpha_gp, rho)
                              + diag_matrix(rep_vector(jitter, N_obs));
    matrix[N_obs, N_obs] L = cholesky_decompose(K);
    vector[N_obs] a         = mdivide_right_tri_low(
                                mdivide_left_tri_low(L, f_obs)', L)';
    return gp_exp_quad_cov(xw_obs, xw_new, alpha_gp, rho)' * a;
  }

  vector gp_conditional_rng(
    array[] vector xw_new,
    array[] vector xw_obs,
    vector         f_obs,
    real           alpha_gp,
    array[] real   rho,
    real           jitter
  ) {
    int N_obs = size(xw_obs);
    int N_new = size(xw_new);
    matrix[N_obs, N_obs] K = gp_exp_quad_cov(xw_obs, alpha_gp, rho)
                              + diag_matrix(rep_vector(jitter, N_obs));
    matrix[N_obs, N_obs] L  = cholesky_decompose(K);
    vector[N_obs] a          = mdivide_right_tri_low(
                                 mdivide_left_tri_low(L, f_obs)', L)';
    matrix[N_obs, N_new] K_cross = gp_exp_quad_cov(xw_obs, xw_new, alpha_gp, rho);
    vector[N_new] mu_cond        = K_cross' * a;
    matrix[N_obs, N_new] V       = mdivide_left_tri_low(L, K_cross);
    matrix[N_new, N_new] K_new   = gp_exp_quad_cov(xw_new, alpha_gp, rho)
                                    + diag_matrix(rep_vector(jitter, N_new));
    return multi_normal_rng(mu_cond, K_new - V' * V);
  }

  vector gp_conditional_noisy_rng(
    array[] vector xw_new,
    array[] vector xw_obs,
    vector         f_obs,
    real           alpha_gp,
    array[] real   rho,
    real           noise_var,
    real           jitter
  ) {
    int N_obs = size(xw_obs);
    int N_new = size(xw_new);
    matrix[N_obs, N_obs] K = gp_exp_quad_cov(xw_obs, alpha_gp, rho)
                              + diag_matrix(rep_vector(noise_var + jitter, N_obs));
    matrix[N_obs, N_obs] L  = cholesky_decompose(K);
    vector[N_obs] a          = mdivide_right_tri_low(
                                 mdivide_left_tri_low(L, f_obs)', L)';
    matrix[N_obs, N_new] K_cross = gp_exp_quad_cov(xw_obs, xw_new, alpha_gp, rho);
    vector[N_new] mu_cond        = K_cross' * a;
    matrix[N_obs, N_new] V       = mdivide_left_tri_low(L, K_cross);
    matrix[N_new, N_new] K_new   = gp_exp_quad_cov(xw_new, alpha_gp, rho)
                                    + diag_matrix(rep_vector(jitter, N_new));
    return multi_normal_rng(mu_cond, K_new - V' * V);
  }

}

data {
  int<lower=1> N_sim;                    
  int<lower=1> N_real;                   
  int<lower=1> N_pred;                   
  int<lower=1> p;                        
  int<lower=1> q;                        

  array[N_sim]  vector[p] x_sim;         
  array[N_sim]  vector[q] w_sim;         
  vector[N_sim] y_sim;                   

  array[N_real] vector[p] x_real;        
  vector[N_real] y_real;                 

  array[N_pred] vector[p] x_pred;        

  real<lower=0, upper=1> w_prior_mean;   
  real<lower=0>          w_prior_sigma;  
}

transformed data {
  real jitter = 1e-8;

  array[N_sim] vector[p+q] xw_sim;
  for (i in 1:N_sim) {
    xw_sim[i, 1:p]         = x_sim[i];
    xw_sim[i, (p+1):(p+q)] = w_sim[i];
  }

  array[N_real] vector[p] x_delta_obs;
  for (j in 1:N_real) x_delta_obs[j] = x_real[j];
}

parameters {
  // Calibration parameter
  real<lower=0, upper=1> w_real;

  // Emulator GP hyperparameters (x+w space)
  array[p+q] real<lower=0> rho_eta;
  real<lower=1e-6> lambda_eta;

  // Discrepancy GP hyperparameters (x space only)
  array[p] real<lower=0> rho_delta;
  real<lower=1e-6> lambda_delta;

  // Noise
  real<lower=1e-6> lambda_e;      // field observation noise
}

transformed parameters {
  vector[N_real] mu_eta_real;
  matrix[N_sim, N_sim] L_eta;
  matrix[N_real, N_real] L_delta;
  array[N_real] vector[p+q] xw_real;
  for (j in 1:N_real) {
    xw_real[j, 1:p] = x_real[j];
    xw_real[j, (p+1):(p+q)] = rep_vector(w_real, q);
  }
  {
    matrix[N_sim, N_sim] K_eta = gp_exp_quad_cov(xw_sim, 1.0/sqrt(lambda_eta), rho_eta)
                                  + diag_matrix(rep_vector(jitter, N_sim));
    L_eta = cholesky_decompose(K_eta);
    vector[N_sim] alpha_vec = mdivide_right_tri_low(
                               mdivide_left_tri_low(L_eta, y_sim)', L_eta)';
    matrix[N_sim, N_real] K_cross = gp_exp_quad_cov(xw_sim, xw_real,
                                     1.0/sqrt(lambda_eta), rho_eta);
    mu_eta_real = K_cross' * alpha_vec;
    matrix[N_real, N_real] K_delta = gp_exp_quad_cov(x_delta_obs, 1.0/sqrt(lambda_delta), rho_delta)
                                      + diag_matrix(rep_vector(1.0/lambda_e + jitter, N_real));
    L_delta = cholesky_decompose(K_delta);
  }
}

model {
  rho_eta   ~ inv_gamma(5, 1);
  lambda_eta ~ gamma(10, 10);
  rho_delta ~ inv_gamma(5, 1);
  lambda_delta ~ gamma(10, 0.3);
  lambda_e ~ gamma(10, 0.03);
  w_real      ~ normal(w_prior_mean, w_prior_sigma);
  target += multi_normal_cholesky_lpdf(y_sim | rep_vector(0.0, N_sim),
                                        L_eta);
  target += multi_normal_cholesky_lpdf(y_real | mu_eta_real, L_delta);
}

generated quantities {

  array[N_pred] vector[p+q] xw_pred;
  for (j in 1:N_pred) {
    xw_pred[j, 1:p]         = x_pred[j];
    xw_pred[j, (p+1):(p+q)] = rep_vector(w_real, q);
  }
  vector[N_pred] f_eta_pred = gp_conditional_mean(
    xw_pred, xw_sim, y_sim, 1.0/sqrt(lambda_eta), rho_eta, jitter);

  array[N_pred] vector[p] x_delta_pred;
  for (j in 1:N_pred) x_delta_pred[j] = x_pred[j];

  vector[N_real] residuals = y_real - mu_eta_real;
  vector[N_pred] f_delta_pred = gp_conditional_noisy_rng(
      x_delta_pred, x_delta_obs, residuals,
      1.0/sqrt(lambda_delta), rho_delta,
      1.0/lambda_e,
      jitter);

  vector[N_pred] mu_pred = f_eta_pred + f_delta_pred;

  vector[N_pred] y_pred;
  for (j in 1:N_pred)
    y_pred[j] = normal_rng(mu_pred[j], 1.0/sqrt(lambda_e));;
}
