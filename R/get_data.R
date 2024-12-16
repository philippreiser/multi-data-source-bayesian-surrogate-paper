library(SobolSequence)
library(vipor)
library(dplyr)
library(here)

source(file.path(here(),"R/utils.R"))

#' Get the ishigami function
#'
#' @param x1 input parameter 1
#' @param x2 input parameter 2
#' @param x3 input parameter 3
#' @param a (optional) control parameter 1
#' @param b (optional) control parameter 2
#'
#' @return The ishigami value
#' @export
#'
#' @examples
#' get_ishigami(0, 0, 0, 7, 0.1)
get_ishigami <- function(x1, x2, x3, a=7, b=0.1) {
  y = a*sin(x2)^2 + (1 + b*(x3^4))*sin(x1)

  return(y)
}

#' Get the sin2d function
#'
#' @param x1 input parameter 1
#' @param x2 input parameter 2
#'
#' @return The sin2d value
#' @export
#'
#' @examples
#' get_sin2d(0, 0)
get_sin2d <- function(x1, x2) {
  y = x2*sin(x1)

  return(y)
}
get_sin2d_real <- function(x1, x2) {
  y <- rep(1, length(x1))
  y[x1 < 0] <- get_sin2d(x1[x1 < 0], x2[x1 < 0])
  y[x1 >= 0] <- get_sin2d(x1[x1 >= 0], 1)
  return(y)
}

get_sin2d_offset <- function(x1, x2, offset=0.1) {
  y = x2*sin(x1)+offset

  return(y)
}

get_linear2d <- function(x1, x2) {
  y = x1*x2

  return(y)
}

get_logistic <- function(x1, x2, k=1, x0=0, A=0) {
  return(x2 / (1 + exp(-k * (x1 - x0))) + A)
}

get_bump <- function(x, scaling=1, start=-1, end=1) {
  scaling = scaling * exp(1)
  return(scaling*exp(-1/(1-(scale_to_1(x, start, end))**2))*((start <= x)&(x <= end)))
}

get_logistic_bump <- function(x1, x2, scaling=0.05, start=5/3, end=10/3) {
  y = get_logistic(x1, x2, k=1, x0=0) + get_bump(x1, scaling=scaling, start=start, end=end)
}

get_log_trend <- function(x1, x2) {
  return(x2 * log(x1) + 0.01 * x1 + 1)
}

get_log_sin <- function(x1, x2) {
  return(get_log_trend(x1, x2) + sin(0.05 * x1))
}

get_log_sin_hf <- function(x1, x2) {
  return(get_log_trend(x1, x2) + sin(0.1 * x1))
}

# SIR model with two input parameters
# x1: time
# x2: beta
# N: total count
# i0: initial infected
get_sir_2d <- function(x1, x2, gamma = 0.55, N=763, i0=1, mean=TRUE) {  
  # initial conditions
  s0 <- N - i0
  r0 <- 0
  y0 <- c(S = s0, I = i0, R = r0)
  
  # fixed parameters
  phi <- 9.6
  
  epred_cases <- rep(NA, length(x1))
  t <- x1
  t0 <- 0
  beta <- x2
  n_days <- length(t)
  gamma <- rep(gamma, n_days)
  data_sir_sim <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N,
                       beta = beta, gamma = gamma, phi = phi)
  sim_model <- cmdstan_model("stan_code/sir_simu_new_interface.stan")
  sim_model_out <- sim_model$sample(data = data_sir_sim, seed = seed,
                                    chains=1, iter_sampling = 1,
                                    fixed_param = TRUE)
  # TODO: replace with (and extract correct column)
  # epred_cases <- matrix(sim_model_out$draws("y"), byrow = FALSE, ncol = 3, nrow = n_days)
  epred_cases <- c(sim_model_out$draws("y")[, , seq(n_days+1, 2*n_days)])
  pred_cases <- c(sim_model_out$draws("pred_cases"))
  if (mean) {
    return(epred_cases)
  } else {
    return(pred_cases)
  }
}

# use the SIR model without noise but return integers.
# the function uses the first decimal place as probability whether to round up
# adapted from matlab code:
# "https://de.mathworks.com/matlabcentral/answers/
# 1946808-how-to-round-numbers-using-the-decimals-as-a-probability"
get_sir_2d_rounded <- function(x1, x2){
  sir_2d_epred <- get_sir_2d(x1, x2, mean = TRUE)
  sir_2d_epred_rounded <- floor(sir_2d_epred) + ((sir_2d_epred -
            floor(sir_2d_epred)) > runif(length(sir_2d_epred), 0, 1))
  return(sir_2d_epred_rounded)
}

get_sir_2d_misspec <- function(x1, x2, gamma=0.7){
  return(get_sir_2d(x1, x2, gamma, mean = TRUE))
}

# SIR model with three input parameters
# x1: time
# x2: beta
# x3: gamma
# N: total count
# i0: initial infected
get_sir <- function(x1, x2, x3, N=763, i0=1, mean=TRUE) {
  # initial conditions
  s0 <- N - i0
  r0 <- 0
  y0 <- c(S = s0, I = i0, R = r0)

  # fixed parameters
  phi <- 9.6

  epred_cases <- rep(NA, length(x1))
  t <- x1
  t0 <- 0
  beta <- x2
  n_days <- length(t)
  gamma <- x3
  data_sir_sim <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N,
                       beta = beta, gamma = gamma, phi = phi)
  sim_model <- cmdstan_model("stan_code/sir_simu_new_interface.stan")
  sim_model_out <- sim_model$sample(data = data_sir_sim, seed = seed,
                                    chains=1, iter_sampling = 1,
                                    fixed_param = TRUE)
  # add one to avoid zero counts
  # TODO: replace with (and extract correct column)
  # epred_cases <- matrix(sim_model_out$draws("y"), byrow = FALSE, ncol = 3, nrow = n_days)
  epred_cases <- c(sim_model_out$draws("y")[, , seq(n_days+1, 2*n_days)])+1
  pred_cases <- c(sim_model_out$draws("pred_cases"))+1
  if (mean) {
    return(epred_cases)
  } else {
    return(pred_cases)
  }
}

# use the SIR model without noise but return integers.
# the function uses the first decimal place as probability whether to round up
# adapted from matlab code:
# "https://de.mathworks.com/matlabcentral/answers/
# 1946808-how-to-round-numbers-using-the-decimals-as-a-probability"
get_sir_epred_rounded <- function(x1, x2, x3){
  sir_epred <- get_sir(x1, x2, x3, mean = TRUE)
  sir_epred_rounded <- floor(sir_epred) + ((sir_epred -
                                                    floor(sir_epred)) > runif(length(sir_epred), 0, 1))
  return(sir_epred_rounded)
}

# scale to [-1, 1] from [lb, ub]
scale_to_1 <- function(x, lb = -pi, ub = pi) {
  x_scaled <- 2/(ub - lb)*(x - lb) - 1

  return(x_scaled)
}

# scale to [lb, ub] from [-1, 1]
scale_from_1 <- function(x, lb = -pi, ub = pi) {
  x_scaled <- (ub - lb)/2*(x + 1) + lb

  return(x_scaled)
}

get_noise <- function(y, noise_model, sigma) {
  N <- length(y)
  if (noise_model == "none"){
    return(y)
  } else if (noise_model == "normal"){
    return(rnorm(N, mean = y, sd = sigma))
  } else if (noise_model == "negbinom"){
    return(rnbinom(N, mu = y, size = sigma))
  }
}

get_ishigami_df <- function(N, sigma=0, a=7, b=0.1, sample_x="sobol"){
  if (sample_x == "uniform"){
    tmp <- matrix(c(runif(N), runif(N), runif(N)), ncol=3)
  } else if (sample_x == "sobol"){
    tmp <- sobolSequence.points(3, count = N)
  } else if (sample_x == "normal"){
    tmp <- matrix(c(rnorm(N), rnorm(N), rnorm(N)), ncol=3)
  } else if (sample_x == "x1_slice"){
    tmp <- matrix(
      c(seq(from=0, to=1,length.out=N), rep(0.5, N), rep(0.5, N)),
      ncol=3)
  } else if (sample_x == "x2_slice"){
    tmp <- matrix(
      c(rep(0.5, N), seq(from=0, to=1,length.out=N), rep(0.5, N)),
      ncol=3)
  } else if (sample_x == "x3_slice"){
    tmp <- matrix(
      c(rep(0.485, N), rep(0.485, N), seq(from=0, to=1,length.out=N)),
      ncol=3)
  }
  tmp <- tmp %>%
    as.data.frame()%>%
    rename(x1 = V1, x2 = V2, x3 = V3) %>%
    mutate_all(~scale_from_1(scale_to_1(., lb = 0, ub = 1))) %>%
    mutate(
      N = N,
      y = get_ishigami(x1, x2, x3, a, b) + rnorm(N, sd=sigma),
      w1 = scale_to_1(x1),
      w2 = scale_to_1(x2),
      w3 = scale_to_1(x3)
    )
  return(tmp)
}

#' Data generation given 2d simulation/real model.
#'
#' @param N number of data points, i.e. model evaluations
#' @param model (string) simulation/real model description that exists as
#'  function when "get_" is added with two input dimensions (x1, x2)
#' @param sigma (optional) standard deviation of additive normal noise model
#' @param sample_x how to sample the input, options:
#'                 "uniform": sample x1 and x2 uniformly
#'                 "sobol": use quasi-random sobol sequence to generate (x1, x2)
#'                 "normal": sample x1 and x2 from a
#'                           Normal(prior_mean, prior_sigma)
#'                 "x1_slice": generate equidistant points along x1 and set
#'                             x_2 = prior_mean
#'                 "x2_slice": generate equidistant points along x2 and set
#'                             x_1 = prior_mean
#'                 "x1_uniform": sample x1 uniformly and set
#'                             x_2 = prior_mean
#'                 "x2_uniform": sample x2 uniformly and set
#'                             x_1 = prior_mean
#'                 "grid": generate (x1, x2) points along a grid
#'                 "x1_vandercorput_x2_normal": generate x1 using the Van-der-
#'                                              Corput sequence and sample x2
#'                                              from a
#'                                              Normal(prior_mean, prior_sigma)
#'                 "x1_uniform_x2_normal": sample x1 uniformly
#'                                         and sample x2
#'                                         from a Normal(prior_mean, prior_sigma)
#'
#' @param prior_mean mean of the Normal-prior for the input
#' @param prior_sigma standard deviation of the Normal-prior for the input
#' @param x_lims list of (x_lb, x_ub), the global lower and upper bounds of x_1
#' @param x_train_percentages list of (x_lower_percentage, x_upper_percentage),
#'                            the lower and upper percentage of the full space
#'                            defined by (x_lb, x_ub) (e.g. to restrict space of
#'                            real data used for training)
#' @param w_lims list of (w_lb, w_ub), the global lower and upper bounds of x_2
#' @param add_one bool that decides whether one should be added to the output
#' (helpful when using a lognormal likelihood to avoid zero counts in the data)
#'
#' @return data frame of simulation/real data points, where
#'         "x1": input 1 in (x_lb, x_ub)
#'         "x2": input 2 in (w_lb, w_ub)
#'         "w1": input 1 scaled to (-1, 1)
#'         "w2": input 2 scaled to (-1, 1)
#'         "N": number of simulation/real data points
#'         "y": (noisy) simulation/real model response
#' @export
#'
#' @examples
get_2d_model_df <- function(N, model, noise_model, sigma=0, sample_x="sobol",
                            prior_mean=0, prior_sigma=1,
                            x_lims=c(-pi, pi), x_train_percentages=c(0, 1),
                            w_lims=c(-1, 1), scale_prior=TRUE,
                            add_one=FALSE
                            ){
  x1_lb <- x_lims[1]
  x1_ub <- x_lims[2]
  x2_lb <- w_lims[1]
  x2_ub <- w_lims[2]
  x1_lower_percentage <- x_train_percentages[1]
  x1_upper_percentage <- x_train_percentages[2]
  if (scale_prior){
    # scale the mean and prior of a normal distribution from [lb, ub] to [0, 1]
    prior_mean <- scale_from_1(scale_to_1(prior_mean, x2_lb, x2_ub), 0, 1)
    prior_sigma <- prior_sigma / (x2_ub - x2_lb)
  }
  get_model <- get(paste0("get_", model))
  # first the input is sampled on the scale [0, 1]
  if (sample_x == "uniform"){
    tmp <- matrix(c(runif(N), runif(N)), ncol=2)
  } else if (sample_x == "sobol"){
    tmp <- sobolSequence.points(2, count = N)
  } else if (sample_x == "normal"){
    tmp <- matrix(c(rnorm(N, prior_mean, prior_sigma),
                    rnorm(N, prior_mean, prior_sigma)), ncol=2)
  } else if (sample_x == "x1_slice"){
    tmp <- matrix(
      c(seq(from=x1_lower_percentage, to=x1_upper_percentage,length.out=N), rep(prior_mean, N)),
      ncol=2)
  } else if (sample_x == "x2_slice"){
    tmp <- matrix(
      c(rep(prior_mean, N), seq(from=0, to=1,length.out=N)),
      ncol=2)
  } else if (sample_x == "x1_uniform"){
    tmp <- matrix(
      c(runif(N, min = x1_lower_percentage, max=x1_upper_percentage), rep(prior_mean, N)),
      ncol=2)
  } else if (sample_x == "x2_uniform"){
    tmp <- matrix(
      c(rep(prior_mean, N), runif(N)),
      ncol=2)
  } else if (sample_x == "grid"){
    tmp <- expand.grid(V1 = seq(0, 1, length.out = N),
                       V2 = seq(0, 1, length.out = N))
  } else if (sample_x == "x1_vandercorput_x2_normal"){
    tmp <- matrix(c(vanDerCorput(N), rnorm(N, prior_mean, prior_sigma)), ncol=2)
  } else if (sample_x == "x1_uniform_x2_normal"){
    tmp <- matrix(c(runif(N, min = x1_lower_percentage, max = x1_upper_percentage), rnorm(N, prior_mean, prior_sigma)), ncol=2)
  }
  tmp <- tmp %>%
    as.data.frame()%>%
    rename(x1 = V1, x2 = V2) %>%
    mutate(
      x1 = scale_from_1(scale_to_1(x1, lb = 0, ub = 1), lb = x1_lb, ub = x1_ub),
      x2 = scale_from_1(scale_to_1(x2, lb = 0, ub = 1), lb = x2_lb, ub = x2_ub)
    )%>%
    mutate(
      N = N,
      y = get_model(x1, x2),
      w1 = scale_to_1(x1, lb = x1_lb, ub = x1_ub), # rename to x1_scaled
      w2 = scale_to_1(x2, lb = x2_lb, ub = x2_ub) # rename to x2_scaled
    ) %>%
    mutate(
      y_noisy = get_noise(y, noise_model, sigma)
    )
  if (add_one) {
    # add one to avoid zero counts
    tmp$y <- tmp$y + 1
    tmp$y_noisy <- tmp$y_noisy + 1
  }
  return(tmp)
}
scale_prior_mean <- function(w_prior_mean, w_lims){
  scale_to_1(w_prior_mean, w_lims[1], w_lims[2])
}
scale_prior_sigma <- function(w_prior_sigma, w_lims){
  2*w_prior_sigma / (w_lims[2] - w_lims[1])
}

get_surrogate_stan_list <- function(df_sim, df_real, x_idxs, w_idxs, beta,
                                   p, M, poly_idx, w_prior_mean, w_prior_sigma,
                                   w_lims,
                                   scale_prior = TRUE,
                                   surrogate_likelihood = "normal",
                                   surrogate_link = "identity",
                                   y_log_transform = FALSE, y_scale = FALSE,
                                   only_x = FALSE){
  # if beta is zero construct a PCE with input dimension equal to the
  # dimension of x
  if (only_x) {
    M <- length(x_idxs)
    w_idxs <- NULL
    w_prior_mean <- numeric()
    w_prior_sigma <- numeric()
  }
  pce_vars <- get_pce_vars(p, M, poly_idx)
  l_poly_coeffs_mat <- pce_vars[[1]]
  comb <- pce_vars[[2]]
  N_real <- nrow(df_real)
  N_sim <- nrow(df_sim)
  y_real <- df_real$y_noisy
  y_sim <- df_sim$y
  alphas <- get_alphas(beta)
  alpha_sim <- alphas[[1]]
  alpha_real <- alphas[[2]]
  if (scale_prior & (M == 2)){
    # scale the mean and prior of a normal distribution from [lb, ub] to [-1, 1]
    w_prior_mean <- scale_prior_mean(w_prior_mean, w_lims)
    w_prior_sigma <- scale_prior_sigma(w_prior_sigma, w_lims)
  } else if (scale_prior & (M == 3)){
    # scale the mean and prior of a normal distribution from [lb, ub] to [-1, 1]
    w_1_prior_mean <- scale_prior_mean(w_prior_mean$w_1, w_lims$w_1)
    w_1_prior_sigma <- scale_prior_sigma(w_prior_sigma$w_1, w_lims$w_1)
    w_2_prior_mean <- scale_prior_mean(w_prior_mean$w_2, w_lims$w_2)
    w_2_prior_sigma <- scale_prior_sigma(w_prior_sigma$w_2, w_lims$w_2)
    w_prior_mean = c(w_1_prior_mean, w_2_prior_mean)
    w_prior_sigma = c(w_1_prior_sigma, w_2_prior_sigma)

  }
  if (y_log_transform) {
    y_real <- df_real$log_y
    y_sim <- df_sim$log_y
  }
  if (y_scale) {
    y_real <- df_real$y_scaled
    y_sim <- df_sim$y_scaled
  }
  log_link <- 0
  if (surrogate_link == "log") {
    log_link <- 1
  }
  normal_likelihood <- 0
  if (surrogate_likelihood == "normal") {
    normal_likelihood <- 1
  }
  lognormal_likelihood <- 0
  if (surrogate_likelihood == "lognormal") {
    lognormal_likelihood <- 1
  }
  data <- list(
    N_sim = nrow(df_sim),
    sigma_sim_lower = 0.0,
    d = p,
    M = M,
    x_sim = as.matrix((df_sim%>%select(starts_with("w")))[, x_idxs]),
    w_sim = as.matrix((df_sim%>%select(starts_with("w")))[, w_idxs]),
    y_sim = y_sim,
    N_real = nrow(df_real),
    y_real = y_real,
    L = length(x_idxs),
    x_real = as.matrix((df_real%>%select(starts_with("w")))[, x_idxs]),
    x_idxs = x_idxs,
    w_idxs = w_idxs,
    beta = beta,
    alpha_sim = alpha_sim,
    alpha_real = alpha_real,
    l_poly_coeffs = t(l_poly_coeffs_mat),
    comb = comb,
    N_comb = nrow(comb),
    w_prior_mean = w_prior_mean,
    w_prior_sigma = w_prior_sigma,
    log_link = log_link,
    normal_likelihood = normal_likelihood,
    lognormal_likelihood = lognormal_likelihood
  )
  return(data)
}

get_data_save_name <- function(results_path, model, data_type,
                               sample_type=NULL) {
  save_name <- file.path(results_path, paste0(model, "_", data_type))
  if (!is.null(sample_type)) {
    save_name <- paste0(save_name, "_", sample_type)
  }
  save_name <- paste0(save_name, ".Rda")
  return(save_name)
}

