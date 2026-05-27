#include pce_helpers.stan

data {
  int<lower=1> N_real;      // Number of measurements per w_real
  int<lower=1> M;               // dimension of input variable
  int<lower=1> L;               // dimension of known input variable
  vector[N_real] y_real;     // Output exp (measured) variables
  int<lower=1> d;               // Degree of polynomials
  matrix[N_real, L] x_real;       // Known input experimental variables
  array[L] int x_idxs;           // which dimensions of input are x
  array[M-L] int w_idxs;           // which dimensions of input are

  int log_link;                   // should a log link be used?
  int normal_likelihood;          // use a normal likelihood?
  int lognormal_likelihood;      // use a log-normal likelihood

  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; // polynomial selection
  real c_0;                                      // coefficient for constant "polynomial"
  vector[N_comb] c;                                   // coefficients of non-constant polynomials

  vector[M-L] w_prior_mean; // mean of w prior
  vector<lower=0>[M-L] w_prior_sigma; // sigma of w prior
}

parameters {
  real<lower=0> sigma;         // sigma for both simulation and real data
  vector<lower=-1, upper=1>[M-L] w_real;  // estimated inputs for real (measured) variables from pce
}

model {
  matrix[N_real, M] input_real;
  input_real[:, x_idxs] = x_real;
  input_real[:, w_idxs] = rep_matrix(w_real', N_real);
  vector[N_real] mu_pred_pce_real = c_0 + get_PCE(input_real, d,
          l_poly_coeffs, comb, N_comb)*c;
  if (log_link) {
    mu_pred_pce_real = exp(mu_pred_pce_real);
  }

  // sigma ~ exponential(2);
  sigma ~ normal(0, 0.5);
  w_real ~ normal(w_prior_mean, w_prior_sigma);

  if (normal_likelihood) {
    target += normal_lpdf(y_real | mu_pred_pce_real, sigma);
  }
  if (lognormal_likelihood) {
    target += lognormal_lpdf(y_real | mu_pred_pce_real, sigma);
  }
}
