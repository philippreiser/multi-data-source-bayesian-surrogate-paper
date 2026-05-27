functions {
  vector gp_conditional_rng(
    array[] vector xw_pred,
    array[] vector xw_sim,
    vector f_sim,
    real alpha_gp,
    array[] real rho,
    real jitter
  ) {
    int N_sim  = size(xw_sim);
    int N_pred = size(xw_pred);

    matrix[N_sim, N_sim] K_sim = gp_exp_quad_cov(xw_sim, alpha_gp, rho);
    K_sim += diag_matrix(rep_vector(jitter, N_sim));
    matrix[N_sim, N_sim] L = cholesky_decompose(K_sim);

    matrix[N_sim, N_pred] K_cross = gp_exp_quad_cov(xw_sim, xw_pred, alpha_gp, rho);

    matrix[N_pred, N_pred] K_pred = gp_exp_quad_cov(xw_pred, alpha_gp, rho);
    K_pred += diag_matrix(rep_vector(jitter, N_pred));

    vector[N_sim] a = mdivide_right_tri_low(
      mdivide_left_tri_low(L, f_sim)', L)';
    vector[N_pred] mu_pred = K_cross' * a;

    matrix[N_sim, N_pred] V = mdivide_left_tri_low(L, K_cross);
    matrix[N_pred, N_pred] cov_pred = K_pred - V' * V;

    return multi_normal_rng(mu_pred, cov_pred);
  }
}

data {
  int<lower=1> N_sim;
  int<lower=1> N_real;
  int<lower=1> N_pred;
  int<lower=1> p;
  int<lower=1> q;

  array[N_real] vector[p] x_real;
  vector[N_real] y_real;

  array[N_pred] vector[p] x_pred;

  array[N_sim] vector[p] x_sim; 
  array[N_sim] vector[q] w_sim;
  vector[N_sim] f_sim;
  array[p + q] real<lower=0> rho;
  real<lower=0> alpha_gp;

  real<lower=0> w_prior_mean;
  real<lower=0, upper=1> w_prior_sigma;
}

transformed data {
  real jitter = 1e-8;

  array[N_sim] vector[p + q] xw_sim;
  for (i in 1:N_sim) {
    xw_sim[i, 1:p]         = x_sim[i];
    xw_sim[i, (p+1):(p+q)] = w_sim[i];
  }

  matrix[N_sim, N_sim] K_base = gp_exp_quad_cov(xw_sim, alpha_gp, rho)
                                 + diag_matrix(rep_vector(jitter, N_sim));
  matrix[N_sim, N_sim] L_base  = cholesky_decompose(K_base);
  vector[N_sim] alpha_vec = mdivide_right_tri_low(
    mdivide_left_tri_low(L_base, f_sim)', L_base)';
}

parameters {
  real<lower=0, upper=1> w_real;

  real<lower=0> sigma;
}

transformed parameters {
  vector[N_real] mu_real;
  {
    array[N_real] vector[p + q] xw_real;
    for (j in 1:N_real) {
      xw_real[j, 1:p]         = x_real[j];
      xw_real[j, (p+1):(p+q)] = rep_vector(w_real, q);
    }
    matrix[N_sim, N_real] K_cross = gp_exp_quad_cov(xw_sim, xw_real, alpha_gp, rho);
    mu_real = K_cross' * alpha_vec;
  }
}

model {
  w_real ~ normal(w_prior_mean, w_prior_sigma);
  sigma  ~ normal(0, 0.5);

  target += normal_lpdf(y_real | mu_real, sigma);
}

generated quantities {
  array[N_pred] vector[p + q] xw_pred;
  for (i in 1:N_pred) {
    xw_pred[i, 1:p]         = x_pred[i];
    xw_pred[i, (p+1):(p+q)] = rep_vector(w_real, q);
  }

  vector[N_pred] f_pred = gp_conditional_rng(
    xw_pred, xw_sim, f_sim, alpha_gp, rho, jitter
  );

  vector[N_pred] y_pred;
  for (j in 1:N_pred)
    y_pred[j] = normal_rng(f_pred[j], sigma);

  vector[N_pred] log_lik;
  for (j in 1:N_pred)
    log_lik[j] = normal_lpdf(y_pred[j] | f_pred[j], sigma);
}
