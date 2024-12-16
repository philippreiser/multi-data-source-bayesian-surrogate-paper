functions {
  vector sir(real t, vector y, int N, real beta, real gamma) {
    vector[3] dy_dt;

    real S = y[1];
    real I = y[2];
    real R = y[3];

    dy_dt[1] = -beta * I * S / N;
    dy_dt[2] = beta * I * S / N - gamma * I;
    dy_dt[3] = gamma * I;

    return dy_dt;
  }
}
data {
  int<lower=1> n_days;
  vector[3] y0;
  real t0;
  array[n_days] real ts;
  int N;
  array[n_days] real<lower=0> beta;
  array[n_days] real<lower=0> gamma;
  real<lower=0> phi;
}
generated quantities {
  array[n_days] real pred_cases;
  array[n_days] vector<lower=0>[3] y;
  array[1] real t;
  for (n in 1:n_days) {
    t[1] = ts[n];
    y[n] = to_vector(ode_rk45(sir, y0, t0, t, N,
                 beta[n], gamma[n])[1]);
  }
  pred_cases = neg_binomial_2_rng(to_vector(y[:, 2]) + 1e-9, phi);
}
