data {
  // --- Simulation data ---
  int<lower=1> N_sim;
  int<lower=1> p;
  int<lower=1> q;
  array[N_sim] vector[p] x_sim;         
  array[N_sim] vector[q] w_sim;
  vector[N_sim] y_sim;         

  // --- Real data ---
  int<lower=1> N_real;
  array[N_real] vector[p] x_real;
  vector[N_real] y_real;

  real<lower=0> alpha_sim;
  real<lower=0> alpha_real;

  real<lower=0> w_prior_mean;
  real<lower=0, upper=1> w_prior_sigma;
}

transformed data {
  real jitter = 1e-8;
  int N_total = N_sim + N_real;
}

parameters {
  array[p + q] real<lower=0> rho;
  real<lower=0> alpha_gp;

  real<lower=0> sigma;

  vector[N_sim] eta_sim_std;

  real<lower=0, upper=1> w_real;
}

transformed parameters {
  vector[N_sim] f_sim;
  vector[N_real] mu_real;
  {
    array[N_sim] vector[p + q] xw_sim;
    for (i in 1:N_sim) {
      xw_sim[i, 1:p]         = x_sim[i];
      xw_sim[i, (p+1):(p+q)] = w_sim[i];
    }

    array[N_real] vector[p + q] xw_real;
    for (j in 1:N_real) {
      xw_real[j, 1:p]         = x_real[j];
      xw_real[j, (p+1):(p+q)] = rep_vector(w_real, q);
    }

    matrix[N_sim, N_sim] K_sim = gp_exp_quad_cov(xw_sim, alpha_gp, rho);
    K_sim += diag_matrix(rep_vector(jitter, N_sim));
    matrix[N_sim, N_sim] L = cholesky_decompose(K_sim);
    f_sim = L * eta_sim_std;

    matrix[N_sim, N_real] K_cross = gp_exp_quad_cov(xw_sim, xw_real, alpha_gp, rho);

    vector[N_sim] alpha_vec = mdivide_right_tri_low(
      mdivide_left_tri_low(L, f_sim)', L)';
    mu_real = K_cross' * alpha_vec;
  }
}

model {
  // --- Priors ---
  rho      ~ inv_gamma(5, 5);
  alpha_gp ~ std_normal();
  sigma    ~ normal(0, 0.5);
  eta_sim_std ~ std_normal();
  w_real   ~ normal(w_prior_mean, w_prior_sigma);

  // --- Power-scaled likelihoods ---
  target += alpha_sim  * normal_lpdf(y_sim  | f_sim,   sigma);
  target += alpha_real * normal_lpdf(y_real | mu_real, sigma);
}
