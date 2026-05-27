library(orthopolynom)
library(here)
library(rlist)
library(posterior)
Rcpp::sourceCpp(file.path(here(),"R/PCE_helpers.cpp"))

# code adapted from Paul-Christian Bürkner
# https://github.com/paul-buerkner/Bayesian-sparse-PCE
get_pce_vars <- function(p, M, idx=NULL){
  l_polys <- legendre.polynomials(p, normalized=TRUE)
  l_poly_coeffs <- polynomial.coefficients(l_polys)
  l_poly_coeffs_mat <- matrix(0, p+1, p+1)
  for (i in 1:(p+1)){
    for (j in 1:length(l_poly_coeffs[[i]])){
      l_poly_coeffs_mat[i, j] = l_poly_coeffs[[i]][j]
    }
    if (i+1<=p+1){
      l_poly_coeffs_mat[i, seq(i+1, p+1)] <- 0
    }
  }
  # Comb & Poly selection
  comb <- poly_idx_cpp(p, M)
  # first column is the constant polynomial f0
  comb <- comb[-1, , drop = FALSE]
  rownames(comb) <- seq_len(nrow(comb))
  if (!is.null(idx)) {
    # select only desired polynomials
    comb <- comb[idx, , drop = FALSE]
  }
  list(l_poly_coeffs_mat, comb)
}

# code from Paul-Christian Bürkner
# https://github.com/paul-buerkner/Bayesian-sparse-PCE/blob/main/PCE_helpers.R
legendre_polynomials <- function(p) {
  orthopolynom::legendre.polynomials(p, normalized = TRUE)
}

# code from Paul-Christian Bürkner
# https://github.com/paul-buerkner/Bayesian-sparse-PCE/blob/main/PCE_helpers.R
# variance of evaluated polynomials from a number of variables
polynomial_variance_legendre <- function(M) {
  out <- 0.4900143 / (2 ^ (M - 1))
  # if (!is.null(cols)) {
  #   out <- rep(out, length(cols))
  # }
  out
}

# code from Paul-Christian Bürkner
# https://github.com/paul-buerkner/Bayesian-sparse-PCE/blob/main/PCE_helpers.R
PCE <- function(..., p = 10, idx = NULL, scale = TRUE,
                poly = legendre_polynomials) {
  dots <- list(...)
  N <- length(dots[[1]])
  M <- length(dots)
  comb <- poly_idx_cpp(p, M)
  # first column is the constant polynomial f0
  comb <- comb[-1, , drop = FALSE]
  rownames(comb) <- seq_len(nrow(comb))
  if (!is.null(idx)) {
    # select only desired polynomials
    comb <- comb[idx, , drop = FALSE]
  }
  if (is.function(poly)) {
    poly <- poly(p)
    poly <- replicate(M, poly, simplify = FALSE)
  }
  stopifnot(is.list(poly) && length(poly) == M)
  out <- matrix(1, N, nrow(comb))
  for (i in seq_len(nrow(comb))) {
    for (j in seq_len(ncol(comb))) {
      out[, i] <- out[, i] * predict(poly[[j]][[comb[i,j] + 1]], dots[[j]])
    }
  }
  if (scale) {
    # TODO: make more general
    out <- out / sqrt(polynomial_variance_legendre(M))
  }
  colnames(out) <- rownames(comb)
  out
}


get_point_elpd <- function(log_liks){
  point_elpd <- matrixStats::colLogSumExps(log_liks) - log(nrow(log_liks))
  return(point_elpd)
}

#' Convert a list of betas (used for mixture_likelihood) into a list of
#' alpha_sims (alpha_T) and alpha_reals (alpha_R) which are used for
#' power_scaling. beta = 1 means full weight on simulation data and beta = 0
#' means full weight on real data. Scaling factors are calculated as:
#' alpha_T = beta/(1 - beta), if beta < 0.5; alpha_T = 1 else
#' alpha_R = (1 - beta)/beta, if beta > 0.5; alpha_T = 1 else
#'
#' @param betas: vector of betas, where beta is the weighting factor for the
#'               likelihood
#'
#' @return list(alpha_sims, alpha_reals)
#'         - alpha_sim: scaling factor of the likelihood for the simulation data
#'         - alpha_real: scaling factor of the likelihood for the real data
#' @export
#'
#' @examples
get_alphas <- function(betas){
  alpha_sims <- rep(1, length(betas))
  alpha_reals <- rep(1, length(betas))
  alpha_sims[betas < 0.5] <- betas[betas < 0.5] / (1 - betas[betas < 0.5])
  alpha_reals[betas > 0.5] <- (1 - betas[betas > 0.5]) / betas[betas > 0.5]
  return(list(alpha_sims, alpha_reals))
}

#' Get the experiment name which is used as a directory name to store/load
#' results.
#'
#' @param config named list loaded from config-file containing experiment
#'               parameters
#'
#' @return string describing the experiment name
#' @export
#'
#' @examples
get_experiment_name <- function(config){
  experiment_name <- sprintf(paste0("%s_%s_N_sim-%d_likelihood-%s_link-%s_ylogtransform-%s_real_train_upper_%s"),
                             config$data$model_sim,
                             config$data$model_real,
                             config$data$N_sim,
                             config$surrogate_model$likelihood,
                             config$surrogate_model$link,
                             config$data$y_log_transform,
                             config$data$x_train_percentages[2])
  if (config$data$model_real == "covid19_italy_first_wave"){
    experiment_name <- sprintf(paste0("%s_%s_N_sim-%d_N_real-%d_likelihood-%s"),
                               config$data$model_sim,
                               config$data$model_real,
                               config$data$N_sim,
                               config$data$N_real,
                               config$surrogate_model$likelihood
                               )
  }
  return(experiment_name)
}

#' Get the file name of the fitted stan model of the one-step-procedure.
#'
#' @param results_path path to the results
#' @param data_integration_scheme power_scaling or mixture_likelihood
#' @param N_real number of real data points
#' @param beta weighting factor for one-step
#'
#' @return string of file name
#' @export
#'
#' @examples
get_fit_file_name <- function(results_path, data_integration_scheme, N_real,
                              beta){
  file = file.path(results_path, sprintf("surrogate_%s_N_real_%s_beta_%s_model_fit",
                                         data_integration_scheme, N_real, beta))
  return(file)
}

# Calculate the posterior prediction mean + CI given a fitted cmdstan model
#' Sample posterior predictive draws of a fitted cmdstan model and store them
#' in a data frame, either as subsampled draws or as as mean + upper + lower
#' quantile.
#'
#' @param model cmdstan model
#' @param fit fitted cmdstan mdoel
#' @param data stan data on which the model is evaluated on
#' @param pred_var prediction variable
#' @param gt_data_sim data frame for inputs
#' @param var_dim input variable
#' @param prob percentile interval width
#' @param spaghetti bool. Should the sampled draws be returned instead of the
#'                  summary (mean + lower + upper quantile)?
#' @param ndraws number of subsampled draws (only used if spaghetti is TRUE)
#' @param draws_gq (optional). If specified, the model is not evaluated.
#'
#' @return data frame containing the posterior predictions. If spaghetti is
#' FALSE, the summary of the posterior predictive (mean + CI) is returned.
#' If spaghetti is TRUE, subsampled posterior predictive draws along with
#' omega draws are returned.
#' @export
#'
#' @examples
get_posterior_pred <- function(model, fit, data, pred_var, gt_data_sim, var_dim,
                               prob = 0.9, spaghetti = FALSE, ndraws = 100,
                               draws_gq = NULL){
  # TODO: implement here the resampling for posterior predicitve stacking
  # i.e. as arg: stacking, beta
  if (is.null(draws_gq)){
    fit_gq <- model$generate_quantities(fit, data = data, seed = 100,
                                        parallel_chains=fit$num_chains())
    # extract the posterior predictive draws along with the corresponding
    # draws of the prior and posterior of w

    # remove and uncomment below
    # draws_gq <- fit_gq$draws(c(pred_var))
    draws_gq <- bind_draws(fit_gq$draws(c(pred_var, "w_prior_sample")),
                            fit$draws("w_real[1]"))
  }
  if (!spaghetti){
    posterior_pred <- draws_gq %>%
      as_draws_df() %>%
      as_tibble() %>%
      select(starts_with(pred_var)) %>%
      apply(2, quantile, c((1 - prob)/2 , 0.5, 1 - (1 - prob)/2)) %>%
      t() %>%
      data.frame(x_plot = gt_data_sim[[var_dim]], .)
  } else{
    x_plot_df <- data.frame(x_plot = gt_data_sim[[var_dim]])%>%
      mutate(x_id = row_number())

    posterior_pred <- resample_draws(draws_gq, ndraws = ndraws) %>%
      as_draws_df() %>%
      as_tibble() %>%
      mutate(w_id = row_number())%>%
      select(starts_with(paste0(pred_var)), "w_prior_sample", "w_real[1]", "w_id")%>%
      pivot_longer(
        cols = !c("w_prior_sample", "w_real[1]", "w_id"),
        names_to = "x_id",
        names_transform = list(x_id = function(s) as.integer(gsub("\\D", "", s))),
        values_to = "y_pred"
      )%>%
      left_join(y = x_plot_df, by="x_id")
  }

  return(posterior_pred)
}

get_inner_mcmc_draws <- function(fit, variables, chains, iter_sampling, inner_chains){
  if (inner_chains == fit$num_chains()) {
    inner_draws <- as_draws_matrix(fit$draws(variables))
    return(inner_draws)
  }
  stopifnot(
    chains*iter_sampling*inner_chains == fit$num_chains(),
    iter_sampling == nrow(fit$draws())
  )
  inner_draws <- fit$draws(variables)[inner_iter_sampling, seq(inner_chains, chains*iter_sampling*inner_chains, inner_chains), ]
  inner_draws <- as_draws_matrix(inner_draws)
  return(inner_draws)
}

# calculate the posterior pred in R without relying on the generated quantities
# from stan
calc_posterior_pred <- function(surrogate_model_fit, data, gt_data_sim, var_dim,
                                beta,
                                surrogate_inference_model_fit = NULL,
                                ndraws = 100) {
  # extract draws from cmdstanr fit
  c_draws <- as_draws_matrix(surrogate_model_fit$draws("c"))
  c_0_draws <- as_draws_matrix(surrogate_model_fit$draws("c_0"))
  w_real_draws <- as_draws_matrix(surrogate_model_fit$draws("w_real[1]"))
  if (beta > 0.5){
    w_real_draws <- get_inner_mcmc_draws(surrogate_inference_model_fit, "w_real[1]")
  }
  sigma_draws <- as_draws_matrix(surrogate_model_fit$draws("sigma"))
  if (!is.null(surrogate_inference_model_fit)) {
    sigma_two_step_draws <- get_inner_mcmc_draws(surrogate_inference_model_fit, "sigma")
  }
  S <- length(sigma_draws)
  S_subset <- ndraws
  pred_dfs <- vector("list", length(S_subset))
  x1 <- gt_data_sim[[var_dim]]
  for (s in (1:S_subset)){
    error_sigma <- rnorm(1, 0, sigma_draws[s])
    error_sigma_two_step <- rnorm(1, 0, sigma_two_step_draws[s])
    x2 <- rep(w_real_draws[s], nrow(gt_data_sim))
    x_design <- PCE(x1, x2, p=config$surrogate_model$poly_degree, idx=config$surrogate_model$poly_idx, scale=FALSE)
    mu_pred <- x_design%*%t(c_draws[s,])+c_0_draws[s]
    pred <- mu_pred + error_sigma
    pred_dfs[[s]] <- data.frame(x_plot = x1, mu_pred = c(mu_pred),
                                pred = c(pred),
                                w_id = s, w_draw = x2,
                                sigma_draw = sigma_two_step_draws[s],
                                c_0_draw = c_0_draws[s])
    if (!is.null(surrogate_inference_model_fit)) {
      pred_two_step <- mu_pred + error_sigma_two_step
      pred_dfs[[s]]$pred_two_step <- c(pred_two_step)
    }
  }
  pred_df <- bind_rows(pred_dfs)
  return(pred_df)
}

#' Calculate the mean of posterior predictive (mu_pred)
#' and the posterior predictive (pred)
#' given (posterior) draws of c, c_0, omega and sigma.
#'
#' @param c_draws matrix of PCE coefficents
#' @param c_0_draws PCE intercept draws matrix
#' @param w_draws w_draws matrix
#' @param sigma_draws sigma_draws matrix
#' @param x input vector
#' @param number_draws (optional) use only the first number_draws draws instead
#' of all draws
#'
#' @return pred_df
#' @export
#'
#' @examples
calc_post_pred <- function(c_draws, c_0_draws, w_draws, sigma_draws,
    x, number_draws=NULL, link="identity", likelihood="normal",
    y_log_transform=FALSE) {
  if (number_draws == 0) {
    return(data.frame())
  }
  S <- nrow(c_draws)
  pred_dfs <- vector("list", length(S))
  if (!is.null(number_draws)){
    S <- number_draws
  }
  linpred_matrix <- calc_pce_linpred_matrix(c_draws, c_0_draws, w_draws,
                                            sigma_draws, x, S, link,
                                            y_log_transform)
  epred_matrix <- calc_pce_epred_matrix(linpred_matrix, sigma_draws, likelihood)
  sigma_matrix <- matrix(sigma_draws[1:S], nrow = S, ncol = length(x))
  if (likelihood == "normal") {
    errors <- rnorm(S, mean = 0, sd = sigma_draws[1:S,])
    pred_matrix <- epred_matrix + errors
    if (y_log_transform) {
      pred_matrix <- exp(log(epred_matrix) + errors)
    }
  } else if (likelihood == "negbinom") {
    pred_matrix <- matrix( rnbinom(S * ncol(epred_matrix), mu=c(epred_matrix), size = c(sigma_matrix)),
                           nrow=S, ncol=ncol(epred_matrix))
  } else if (likelihood == "lognormal") {
    pred_matrix <- matrix(rlnorm(S * ncol(epred_matrix), meanlog = c(linpred_matrix), sdlog = c(sigma_matrix)),
                          nrow=S, ncol=ncol(epred_matrix))
  }
  x_matrix <- matrix(x, nrow = S, ncol = length(x), byrow = TRUE)
  c_0_matrix <- matrix(c_0_draws[1:S], nrow = S, ncol = length(x))

  if (is.null(w_draws)) {
    w_id_matrix <- matrix(seq(1, S), nrow = S, ncol = length(x))
    w_matrix <- matrix(NA, nrow = S, ncol = length(x))
  } else {
    w_id_matrix <- matrix(seq(1, S), nrow = S, ncol = length(x))
    w_matrix <- matrix(w_draws[1:S], nrow = S, ncol = length(x))
  }
  pred_df <- data.frame(
    x_plot = c(t(x_matrix)),
    mu_pred = c(t(epred_matrix)),
    pred = c(t(pred_matrix)),
    sigma_draw = c(t(sigma_matrix)),
    w_id = c(t(w_id_matrix)),
    w_draw = c(t(w_matrix)),
    c_0_draw = c(t(c_0_matrix))
  )
  return(pred_df)
}

#' Calculate the mean of the posterior predictive of the PCE given
#' linear predictor and the likelihood.
#'
#' @param linpred_matrix linpred_matrix (S, N)
#' @param likelihood likelihood
#'
#' @return epred_matrix (S, N)
#' @export
#'
#' @examples
calc_pce_epred_matrix <- function(linpred_matrix, sigma_draws, likelihood="normal") {
  S <- nrow(linpred_matrix)
  N <- ncol(linpred_matrix)
  if ((likelihood == "lognormal") && (S!=0)){
    sigma_draws <- sigma_draws[1:S]
    sigma_matrix <- matrix(rep(sigma_draws, times=N), ncol=N)
    epred_matrix <- exp(linpred_matrix + sigma_matrix^2 / 2)  
  }
  else {
    epred_matrix <- linpred_matrix
  }
  return(epred_matrix)
}

#' Calculate the linear predictor of the PCE given
#' S (posterior) draws of c, c_0, omega and sigma. Optionally transformed.
#'
#' @param c_draws PCE coefficents draws matrix (S, d)
#' @param c_0_draws PCE intercept draws matrix (S, 1)
#' @param w_draws w_draws matrix (S, 1)
#' @param x input vector (N)
#' @param number_draws (optional) use only the first number_draws draws instead
#' of all draws
#'
#' @return linpred_matrix (S, N)
#' @export
#'
#' @examples
calc_pce_linpred_matrix <- function(c_draws, c_0_draws, w_draws, sigma_draws, x,
                                  number_draws=NULL, link="identity",
                                  y_log_transform=FALSE) {
  S <- nrow(c_draws)
  N <- length(x)
  stopifnot(
    S == nrow(c_0_draws),
    S == nrow(w_draws),
    S >= number_draws
  )
  if (!is.null(w_draws)) {
    stopifnot(S == nrow(w_draws))
  }
  if (!is.null(number_draws)){
    S <- number_draws
  }
  if (S == 0) {
    return(matrix(nrow=0, ncol=N))
  }
  linpred_matrix <- matrix(0, nrow=S, ncol=N)
  for (s in (1:S)){
    if (!is.null(w_draws)){
      if (ncol(w_draws) == 1){
        w <- rep(w_draws[s], length(x))
        x_design <- PCE(x, w, p=config$surrogate_model$poly_degree,
                        idx=config$surrogate_model$poly_idx, scale=FALSE)
      } else if (ncol(w_draws) == 2){
        w1 <- rep(w_draws[s, 1], length(x))
        w2 <- rep(w_draws[s, 2], length(x))
        x_design <- PCE(x, w1, w2, p=config$surrogate_model$poly_degree,
                        idx=config$surrogate_model$poly_idx, scale=FALSE)
      }
      
    } else {
      x_design <- PCE(x, p=config$surrogate_model$poly_degree,
                      idx=config$surrogate_model$poly_idx, scale=FALSE)
    }
    lin_pred <- x_design%*%t(c_draws[s,]) + c_0_draws[s]
    if ((link == "log") || (y_log_transform)){
      lin_pred <- exp(lin_pred)
    }
    linpred_matrix[s, ] <- lin_pred
  }
  return(linpred_matrix)
}

get_lpd_matrix <- function(linpred_matrix, sigma_draws, y, likelihood) {
  S <- nrow(linpred_matrix)
  N <- ncol(linpred_matrix)
  lpd_matrix <- matrix(NA, nrow = S, ncol = N)
  for (s in (1:S)) {
    for (n in (1:N)){
      if (likelihood == "normal"){
        lpd_matrix[s, n] <- dnorm(y[n], linpred_matrix[s, n], sigma_draws[s], log=TRUE)
      } else if (likelihood == "negbinom"){
        lpd_matrix[s, n] <- dnbinom(y[n], mu=linpred_matrix[s, n],
                                    size=sigma_draws[s], log=TRUE)
      } else if (likelihood == "lognormal"){
        lpd_matrix[s, n] <- dlnorm(y[n], linpred_matrix[s, n],
                                    sigma_draws[s], log=TRUE)
      }
    }
  }
  return(lpd_matrix)
}

get_posterior_rmse <- function(epred_matrix, y_matrix) {
  delta_matrix <- epred_matrix - y_matrix
  rmse <- mean(sqrt(colMeans((delta_matrix)^2)))
  return(rmse)
}

#' Calculate the ELPD and RMSE metric of draws and true values
#'
#' @param c_draws c draw matrix
#' @param c_0_draws c_0 draw matrix
#' @param w_draws w draw matrix
#' @param sigma_draws sigma draw matrix
#' @param x input values
#' @param y_true true output values
#' @param y_true_noisy true noisy output values
#'
#' @return list containing the elpd and rmse
#' @export
#'
#' @examples
get_metrics <- function(linpred_matrix, epred_matrix, sigma_draws,
                        y_true, y_true_noisy, likelihood){
  S <- nrow(epred_matrix)
  N <- ncol(epred_matrix)
  lpd_matrix <- get_lpd_matrix(linpred_matrix, sigma_draws, y_true_noisy, likelihood)
  y_matrix <- matrix(y_true, nrow = S, ncol = N, byrow = TRUE)
  rmse <- get_posterior_rmse(epred_matrix, y_matrix)
  ood_metrics <- list(
    "lpd_matrix" = lpd_matrix,
    "elpd" = elpd(lpd_matrix),
    "rmse" = rmse)
  return(ood_metrics)
}

calc_metrics <- function(surrogate_model_fit, sdata_1, beta,
                      surrogate_inference_model_fit=NULL) {
  c_draws <- as_draws_matrix(surrogate_model_fit$draws("c"))
  c_0_draws <- as_draws_matrix(surrogate_model_fit$draws("c_0"))
  if (beta > 0.5) {
    w_real_draws <- get_inner_mcmc_draws(surrogate_inference_model_fit, "w_real[1]")
  } else {
    w_real_draws <- as_draws_matrix(surrogate_model_fit$draws("w_real[1]"))
  }
  sigma_draws <- as_draws_matrix(surrogate_model_fit$draws("sigma"))
  if (!is.null(surrogate_inference_model_fit)) {
    sigma_draws <- get_inner_mcmc_draws(surrogate_inference_model_fit, "sigma")
  }
  S <- length(sigma_draws)
  S_subset <- 1000
  lpd_matrix <- matrix(0, nrow=S_subset, ncol=sdata_1$N_real)
  delta_matrix <- matrix(0, nrow=S_subset, ncol=sdata_1$N_real)
  x1 <- sdata_1$x_real
  y_real <- sdata_1$y_real
  for (s in (1:S_subset)){
    # TODO: replace 4 with number of chains of E-Post
    x2 <- rep(w_real_draws[s], sdata_1$N_real)
    x_design <- PCE(x1, x2, p=config$surrogate_model$poly_degree, idx=config$surrogate_model$poly_idx, scale=FALSE)
    mu_pred <- x_design%*%t(c_draws[s,])+c_0_draws[s]
    sigma_draw <- sigma_draws[s]
    lpd_matrix[s, ] <- dnorm(y_real, mu_pred, sigma_draw, log=TRUE)
    delta_matrix[s, ] <- y_real - mu_pred
  }
  ood_metrics <- list("elpd" = elpd(lpd_matrix),
       "rmse" = mean(sqrt(colMeans((delta_matrix)^2))))
  return(ood_metrics)
}
