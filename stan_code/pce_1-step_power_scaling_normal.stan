#include pce_helpers.stan

data {
  int<lower=1> N_sim;           // Number of input/output simulation pairs
  int<lower=1> N_real;      // Number of measurements per w_real
  int<lower=1> M;               // dimension of input variable
  int<lower=1> L;               // dimension of known input variable
  vector[N_sim] y_sim;          // Output simulation variables
  vector[N_real] y_real;     // Output exp (measured) variables
  int<lower=1> d;               // Degree of polynomials
  matrix[N_sim, L] x_sim;       // (known) input simulation variables
  matrix[N_sim, M-L] w_sim;       // (unknown) input simulation variables
  matrix[N_real, L] x_real;       // Known input experimental variables
  array[L] int x_idxs;           // which dimensions of input are x
  array[M-L] int w_idxs;           // which dimensions of input are w
  real<lower=0> alpha_sim;                 // weighting of simulation data
  real<lower=0> alpha_real;                 // weighting of real data

  int log_link;                   // should a log link be used?
  int normal_likelihood;          // use a normal likelihood?
  int lognormal_likelihood;      // use a log-normal likelihood

  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; // polynomial selection

  vector[M-L] w_prior_mean; // mean of w prior
  vector<lower=0>[M-L] w_prior_sigma; // sigma of w prior
}

transformed data {
  // real sigma = 0.1;
  matrix[N_sim, M] input_sim;
  input_sim[:, x_idxs] = x_sim;
  input_sim[:, w_idxs] = w_sim;
  matrix[N_sim, N_comb] Input_sim = get_PCE(input_sim, d, l_poly_coeffs, comb, N_comb);
}

parameters {
  real c_0;                                      // coefficient for constant "polynomial"
  vector[N_comb] c;                                   // coefficients of non-constant polynomials
  real<lower=0> sigma;         // sigma for both simulation and real data
  vector<lower=-1, upper=1>[M-L] w_real;  // estimated inputs for real (measured) variables from pce
}

model {
  matrix[N_real, M] input_real;
  input_real[:, x_idxs] = x_real;
  input_real[:, w_idxs] = rep_matrix(w_real', N_real);
  vector[N_sim] mu_pred_pce_sim = c_0 + Input_sim*c;
  vector[N_real] mu_pred_pce_real = c_0 + get_PCE(input_real, d,
          l_poly_coeffs, comb, N_comb)*c;
  if (log_link) {
    mu_pred_pce_sim = exp(mu_pred_pce_sim);
    mu_pred_pce_real = exp(mu_pred_pce_real);
  }

  // Prior model
  c_0 ~ normal(0, 5);
  c ~ normal(0, 5);
  // sigma ~ exponential(2);
  sigma ~ normal(0, 0.5);
  w_real ~ normal(w_prior_mean, w_prior_sigma);

  // Observational model
  if (normal_likelihood) {
    target += alpha_sim * normal_lpdf(y_sim | mu_pred_pce_sim, sigma);
    target += alpha_real * normal_lpdf(y_real | mu_pred_pce_real, sigma);
  }
  if (lognormal_likelihood) {
    target += alpha_sim * lognormal_lpdf(y_sim | mu_pred_pce_sim, sigma);
    target += alpha_real * lognormal_lpdf(y_real | mu_pred_pce_real, sigma);
  }
}
