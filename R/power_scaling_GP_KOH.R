library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(tidyverse)
library(viridis)
library(latex2exp)
library(matrixStats)
library(loo)
library(patchwork)

set.seed(42)

# ── file paths ─────────────────────────────────────────────────────────

stan_step1_file <- "stan_code/GP_power_scaling_step1.stan"
stan_step2_file <- "stan_code/GP_power_scaling_step2.stan"
stan_koh_file   <- "stan_code/GP_KOH.stan"

data_dir    <- "results/log_trend_log_sin_N_sim-100_likelihood-normal_link-identity_ylogtransform-FALSE_real_train_upper_0.5 _koh"
results_dir <- file.path("results", "gp_power_scaling")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

betas <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# ── MCMC settings ─────────────────────────────────────────────────────────────

# Power-scaling Step 1
chains_1          <- 4
iter_warmup_1     <- 1000
iter_sampling_1   <- 1000
parallel_chains_1 <- 4
adapt_delta       <- 0.95

# Power-scaling Step 2 (E-Post inner MCMC)
S1_sub            <- 200
inner_chains      <- 1
inner_iter_warmup <- 200
inner_iter_sample <- 50 

# KOH
chains_koh          <- 4
iter_warmup_koh     <- 1000
iter_sampling_koh   <- 1000
parallel_chains_koh <- 4

# Plot settings
pp_ndraws <- 200   # spaghetti lines per method / panel

# ── 1. Load data ─────────────────────────────────────────────────────────────
#   w1        – observable input x
#   w2        – calibration input w  (simulator only, known)
#   y_noisy   – noisy model/measurement output

df_sim      <- readRDS(file.path(data_dir, "log_trend_simulation.Rda"))
df_real     <- readRDS(file.path(data_dir, "log_sin_real_x1_slice.Rda"))
df_real_oos <- readRDS(file.path(data_dir, "log_sin_real_oos_ood_x1_uniform.Rda"))

# Optional: subsample for faster iteration during development
set.seed(123)
# df_sim      <- df_sim[sample(nrow(df_sim),  30), ]
# df_real     <- df_real[sample(nrow(df_real), 10), ]
# df_real_oos <- df_real_oos[sample(nrow(df_real_oos), 50), ]

cat(sprintf("Data: N_sim=%d  N_real=%d  N_pred=%d\n",
            nrow(df_sim), nrow(df_real), nrow(df_real_oos)))

# ── Pre-processing ─────────────────────────────────────────────────────────

p <- 1; q <- 1

y_sim  <- df_sim$y_noisy
x_sim  <- matrix(df_sim$w1, ncol = p)
w_sim  <- matrix(df_sim$w2, ncol = q)
y_real <- df_real$y_noisy
x_real <- matrix(df_real$w1, ncol = p)
x_pred <- matrix(df_real_oos$w1, ncol = p)

y_sim_mu <- mean(y_sim, na.rm = TRUE)
y_sim_sd <- sd(y_sim,   na.rm = TRUE)
y_sim_sc  <- (y_sim  - y_sim_mu) / y_sim_sd
y_real_sc <- (y_real - y_sim_mu) / y_sim_sd

x_all <- rbind(x_real, x_sim)
x_min <- min(x_all); x_max <- max(x_all)
scale_x <- function(x) (x - x_min) / (x_max - x_min)
x_sim_sc  <- scale_x(x_sim)
x_real_sc <- scale_x(x_real)
x_pred_sc <- scale_x(x_pred)

w_min <- min(w_sim); w_max <- max(w_sim)
scale_w <- function(w) (w - w_min) / (w_max - w_min)
w_sim_sc <- scale_w(w_sim)

w_prior_mean  <- 0.5
w_prior_sigma <- 0.5

to_vec_array <- function(mat) lapply(seq_len(nrow(mat)), function(i) mat[i, ])

x_sim_arr  <- to_vec_array(x_sim_sc)
w_sim_arr  <- to_vec_array(w_sim_sc)
x_real_arr <- to_vec_array(x_real_sc)
x_pred_arr <- to_vec_array(x_pred_sc)

# ── Helper functions ───────────────────────────────────────────────────────

get_alphas_gp <- function(beta) {
  if (beta == 0)  return(list(alpha_sim = 0.0,               alpha_real = 1.0))
  if (beta == 1)  return(list(alpha_sim = 1.0,               alpha_real = 0.0))
  if (beta < 0.5) return(list(alpha_sim = beta / (1 - beta), alpha_real = 1.0))
  return(list(alpha_sim = 1.0, alpha_real = (1 - beta) / beta))
}

run_epost_step2 <- function(fit_step1, stan_base_s2,
                             N_sim, N_real, N_pred, p, q,
                             w_prior_mean, w_prior_sigma,
                             S1_sub, inner_iter_warmup, inner_iter_sample,
                             inner_chains) {

  draws_s1 <- as_draws_df(fit_step1$draws(
    c("rho", "alpha_gp", "sigma", "w_real",
      paste0("f_sim[", seq_len(N_sim), "]"))
  ))

  idx       <- round(seq(1, nrow(draws_s1), length.out = S1_sub))
  draws_sub <- draws_s1[idx, ]

  rho_cols  <- grep("^rho\\[",   names(draws_sub), value = TRUE)
  fsim_cols <- grep("^f_sim\\[", names(draws_sub), value = TRUE)

  w_real_out <- numeric(S1_sub)
  sigma_out  <- numeric(S1_sub)
  y_pred_mat <- matrix(NA_real_, S1_sub, N_pred)
  f_pred_mat <- matrix(NA_real_, S1_sub, N_pred)

  cat(sprintf("  E-Post Step 2: %d inner fits ...\n", S1_sub))
  pb <- txtProgressBar(min = 0, max = S1_sub, style = 3)

  for (s in seq_len(S1_sub)) {
    row <- draws_sub[s, ]
    stan_data_s2 <- c(
      stan_base_s2,
      list(
        f_sim         = as.numeric(row[fsim_cols]),
        rho           = as.numeric(row[rho_cols]),
        alpha_gp      = as.numeric(row[["alpha_gp"]]),
        w_prior_mean  = w_prior_mean,
        w_prior_sigma = w_prior_sigma
      )
    )

    init_s2 <- list(list(
      w_real = max(0.01, min(0.99, as.numeric(row[["w_real"]]))),
      sigma  = max(0.01,           as.numeric(row[["sigma"]]))
    ))

    suppressMessages({
      fit_s2 <- mod_step2$sample(
        data          = stan_data_s2,
        init          = init_s2,
        chains        = inner_chains,
        iter_warmup   = inner_iter_warmup,
        iter_sampling = inner_iter_sample,
        refresh       = 0,
        show_messages = FALSE
      )
    })

    d2   <- as_draws_df(fit_s2$draws())
    last <- nrow(d2)
    w_real_out[s] <- as.numeric(d2[last, "w_real"])
    sigma_out[s]  <- as.numeric(d2[last, "sigma"])

    yp_cols <- grep("^y_pred\\[", names(d2), value = TRUE)
    fp_cols <- grep("^f_pred\\[", names(d2), value = TRUE)
    y_pred_mat[s, ] <- as.numeric(d2[last, yp_cols])
    f_pred_mat[s, ] <- as.numeric(d2[last, fp_cols])

    setTxtProgressBar(pb, s)
  }
  close(pb)

  list(w_real_draws  = w_real_out,
       sigma_draws   = sigma_out,
       y_pred_matrix = y_pred_mat,
       f_pred_matrix = f_pred_mat)
}

#' ELPD/obs and RMSE on scaled test data.
compute_metrics <- function(f_pred_matrix, sigma_draws, y_test_sc) {
  S <- nrow(f_pred_matrix); N <- ncol(f_pred_matrix)
  ll <- matrix(NA_real_, S, N)
  for (s in seq_len(S))
    ll[s, ] <- dnorm(y_test_sc, f_pred_matrix[s, ], sigma_draws[s], log = TRUE)
  elpd_per_obs <- mean(matrixStats::colLogSumExps(ll) - log(S))
  y_mat <- matrix(y_test_sc, S, N, byrow = TRUE)
  rmse  <- sqrt(mean((f_pred_matrix - y_mat)^2))
  list(elpd_per_obs = elpd_per_obs, rmse = rmse)
}

# ── Compile Stan models ────────────────────────────────────────────────────

cat("Compiling Stan models ...\n")
mod_step1 <- cmdstan_model(stan_step1_file)
mod_step2 <- cmdstan_model(stan_step2_file)
mod_koh   <- cmdstan_model(stan_koh_file)

# ── Fixed base data for Step 2 ────────────────────────────────────────────

stan_base_step2 <- list(
  N_sim  = nrow(df_sim),
  N_real = nrow(df_real),
  N_pred = nrow(df_real_oos),
  p = p, q = q,
  x_sim  = x_sim_arr, w_sim = w_sim_arr,
  x_real = x_real_arr, y_real = y_real_sc,
  x_pred = x_pred_arr
)

# ── Power-scaling loop over betas ─────────────────────────────────────────

y_test_sc <- (df_real_oos$y_noisy - y_sim_mu) / y_sim_sd

ps_results <- list()

for (beta in betas) {
  cat(sprintf("\n══ Power-Scaling  beta = %.2f ══\n", beta))

  alph <- get_alphas_gp(beta)
  cat(sprintf("  alpha_sim=%.4f  alpha_real=%.4f\n",
              alph$alpha_sim, alph$alpha_real))

  # Step 1: joint power-scaled training
  stan_data_s1 <- list(
    N_sim = nrow(df_sim), N_real = nrow(df_real),
    p = p, q = q,
    x_sim = x_sim_arr, w_sim = w_sim_arr, y_sim = y_sim_sc,
    x_real = x_real_arr, y_real = y_real_sc,
    alpha_sim  = alph$alpha_sim,
    alpha_real = alph$alpha_real,
    w_prior_mean  = w_prior_mean,
    w_prior_sigma = w_prior_sigma
  )

  cat("  Step 1 ...\n")
  fit_s1_path <- file.path(
    results_dir,
    sprintf("gp_ps_step1_beta_%.2f", beta)
  )
  if (file.exists(fit_s1_path)) {
    fit_s1 <- readRDS(fit_s1_path)
  } else {
    fit_s1 <- mod_step1$sample(
      data            = stan_data_s1,
      seed            = 42,
      chains          = chains_1,
      parallel_chains = parallel_chains_1,
      iter_warmup     = iter_warmup_1,
      iter_sampling   = iter_sampling_1,
      adapt_delta     = adapt_delta,
      refresh         = 200,
      init = 0
    )
    fit_s1$save_object(fit_s1_path)
    # Diagnostics
    np_s1 <- nuts_params(fit_s1)
    ggsave(
      file.path(results_dir, sprintf("trace_step1_beta_%.2f.png", beta)),
      mcmc_trace(fit_s1$draws(),
                 pars = c("alpha_gp", "sigma", "w_real"), regex_pars = "^rho",
                 np = np_s1),
      width = 10, height = 8)
    rh <- rhat(fit_s1)
    ggsave(
      file.path(results_dir, sprintf("rhat_step1_beta_%.2f.png", beta)),
      mcmc_rhat(rh[grepl("alpha_gp|sigma|w_real|rho", names(rh))]) +
        yaxis_text(hjust = 1),
      width = 8, height = 6)
  }

  # Step 2: E-Post refinement of (w_real, sigma) on real data only
  fit_s2_path <- file.path(results_dir,
                           sprintf("gp_ps_epost_beta_%.2f.rds", beta))
  if (file.exists(fit_s2_path)) {
    epost <- readRDS(fit_s2_path)
  } else {
    epost <- run_epost_step2(
      fit_step1 = fit_s1, stan_base_s2 = stan_base_step2,
      N_sim = nrow(df_sim), N_real = nrow(df_real), N_pred = nrow(df_real_oos),
      p = p, q = q,
      w_prior_mean = w_prior_mean, w_prior_sigma = w_prior_sigma,
      S1_sub = S1_sub,
      inner_iter_warmup = inner_iter_warmup,
      inner_iter_sample = inner_iter_sample,
      inner_chains      = inner_chains
    )
    saveRDS(epost, fit_s2_path)
  }

  metrics <- compute_metrics(epost$f_pred_matrix, epost$sigma_draws, y_test_sc)
  cat(sprintf("  ELPD/obs = %.4f   RMSE = %.4f\n",
              metrics$elpd_per_obs, metrics$rmse))

  ps_results[[as.character(beta)]] <- list(
    beta = beta,
    alpha_sim = alph$alpha_sim, alpha_real = alph$alpha_real,
    w_real_mean  = mean(epost$w_real_draws),
    w_real_sd    = sd(epost$w_real_draws),
    sigma_mean   = mean(epost$sigma_draws),
    sigma_sd     = sd(epost$sigma_draws),
    elpd_per_obs = metrics$elpd_per_obs,
    rmse         = metrics$rmse,
    epost        = epost
  )
}

# ── KOH model ─────────────────────────────────────────────────────────────

cat("\n══ KOH ══\n")

fit_koh_path <- file.path(results_dir, "gp_koh_fit")
if (file.exists(fit_koh_path)) {
  fit_koh <- readRDS(fit_koh_path)
} else {
  stan_data_koh <- list(
    N_sim = nrow(df_sim), N_real = nrow(df_real), N_pred = nrow(df_real_oos),
    p = p, q = q,
    x_sim  = x_sim_arr, w_sim = w_sim_arr, y_sim = y_sim_sc,
    x_real = x_real_arr, y_real = y_real_sc,
    x_pred = x_pred_arr,
    w_prior_mean  = w_prior_mean,
    w_prior_sigma = w_prior_sigma
  )
  
  fit_koh <- mod_koh$sample(
    data            = stan_data_koh,
    seed            = 42,
    chains          = chains_koh,
    parallel_chains = parallel_chains_koh,
    iter_warmup     = iter_warmup_koh,
    iter_sampling   = iter_sampling_koh,
    adapt_delta     = adapt_delta,
    refresh         = 200,
    init            = 0
  )
  fit_koh$save_object(fit_koh_path) 
}

variables_koh <- c(
  "w_real", "lambda_eta", "lambda_delta", "lambda_e"
)
np_koh <- nuts_params(fit_koh)
ggsave(file.path(results_dir, "trace_koh.png"),
  mcmc_trace(fit_koh$draws(),
    pars        = variables_koh,
    regex_pars  = c("rho_eta", "rho_delta"),
    np          = np_koh),
  width = 12, height = 10)
rh_koh <- rhat(fit_koh)
ggsave(file.path(results_dir, "rhat_koh.png"),
  mcmc_rhat(rh_koh[grepl(
    "lambda_eta|lambda_delta|lambda_e|w_real|rho_eta|rho_delta",
    names(rh_koh))]) + yaxis_text(hjust = 1),
  width = 8, height = 8)

cat("KOH summary:\n")
print(fit_koh$summary(variables_koh))

# Extract KOH posterior draws needed for plots and metrics
koh_draws      <- as_draws_df(fit_koh$draws())
koh_fpred_cols <- grep("^mu_pred\\[", names(koh_draws), value = TRUE)
koh_ypred_cols <- grep("^y_pred\\[", names(koh_draws), value = TRUE)
koh_sigma_vec  <- 1 / sqrt(koh_draws[["lambda_e"]])
koh_w_vec      <- koh_draws[["w_real"]]

koh_f_mat <- as.matrix(koh_draws[, koh_fpred_cols])
koh_y_mat <- as.matrix(koh_draws[, koh_ypred_cols])
koh_metrics <- compute_metrics(koh_f_mat, koh_sigma_vec, y_test_sc)
cat(sprintf("KOH  ELPD/obs = %.4f   RMSE = %.4f\n",
            koh_metrics$elpd_per_obs, koh_metrics$rmse))

# ── Posterior Predictive Weighting ───────────────────────────────────────

cat("\n══ Posterior Predictive Weighting ══\n")

# --- Extract posterior draws from the two boundary fits ----------------

ep_real <- ps_results[["0"]]$epost   # beta = 0  (real-data-only)
ep_sim  <- ps_results[["1"]]$epost   # beta = 1  (sim-data-only)

S_real <- nrow(ep_real$f_pred_matrix)
S_sim  <- nrow(ep_sim$f_pred_matrix)

# --- Find optimal stacking weights on the test set --------------------

get_point_lpd <- function(f_mat, sigma_vec, y_test) {
  S <- nrow(f_mat); N <- ncol(f_mat)
  ll <- matrix(NA_real_, S, N)
  for (s in seq_len(S))
    ll[s, ] <- dnorm(y_test, f_mat[s, ], sigma_vec[s], log = TRUE)
  matrixStats::colLogSumExps(ll) - log(S)
}

lpd_real_pts <- get_point_lpd(ep_real$f_pred_matrix, ep_real$sigma_draws, y_test_sc)
lpd_sim_pts  <- get_point_lpd(ep_sim$f_pred_matrix,  ep_sim$sigma_draws,  y_test_sc)

# stacking_weights expects an N x K matrix of point log-predictive densities
lpd_matrix_stack <- cbind(lpd_real_pts, lpd_sim_pts)   # N x 2
stacking_wts     <- loo::stacking_weights(lpd_matrix_stack)
pw_optimal_beta  <- stacking_wts[2]
cat(sprintf("  Optimal stacking weight: %.4f\n", pw_optimal_beta))

pw_results <- list()

for (beta in betas) {
  n_real_draws <- round((1 - beta) * S_real)
  n_sim_draws  <- round(beta       * S_sim)

  # subsample draws from each surrogate
  idx_real <- if (n_real_draws > 0) seq_len(n_real_draws) else integer(0)
  idx_sim  <- if (n_sim_draws  > 0) seq_len(n_sim_draws)  else integer(0)

  f_real <- if (length(idx_real) > 0) ep_real$f_pred_matrix[idx_real, , drop = FALSE] else matrix(nrow = 0, ncol = ncol(ep_real$f_pred_matrix))
  f_sim  <- if (length(idx_sim)  > 0) ep_sim$f_pred_matrix[idx_sim,  , drop = FALSE] else matrix(nrow = 0, ncol = ncol(ep_sim$f_pred_matrix))
  s_real <- if (length(idx_real) > 0) ep_real$sigma_draws[idx_real]  else numeric(0)
  s_sim  <- if (length(idx_sim)  > 0) ep_sim$sigma_draws[idx_sim]    else numeric(0)

  y_real_draws <- if (length(idx_real) > 0) ep_real$y_pred_matrix[idx_real, , drop = FALSE] else matrix(nrow = 0, ncol = ncol(ep_real$y_pred_matrix))
  y_sim_draws  <- if (length(idx_sim)  > 0) ep_sim$y_pred_matrix[idx_sim,  , drop = FALSE] else matrix(nrow = 0, ncol = ncol(ep_sim$y_pred_matrix))
  w_real_draws_real <- if (length(idx_real) > 0) ep_real$w_real_draws[idx_real] else numeric(0)
  w_real_draws_sim  <- if (length(idx_sim)  > 0) ep_sim$w_real_draws[idx_sim]   else numeric(0)

  f_pw <- rbind(f_real, f_sim)
  y_pw <- rbind(y_real_draws, y_sim_draws)
  s_pw <- c(s_real, s_sim)

  metrics_pw <- compute_metrics(f_pw, s_pw, y_test_sc)
  cat(sprintf("  beta=%.2f  ELPD/obs=%.4f  RMSE=%.4f\n",
              beta, metrics_pw$elpd_per_obs, metrics_pw$rmse))

  pw_results[[as.character(beta)]] <- list(
    beta          = beta,
    f_pred_matrix = f_pw,
    y_pred_matrix = y_pw,
    sigma_draws   = s_pw,
    w_real_draws  = c(w_real_draws_real, w_real_draws_sim),
    elpd_per_obs  = metrics_pw$elpd_per_obs,
    rmse          = metrics_pw$rmse
  )
}

# ── Summary table ─────────────────────────────────────────────────────────

ps_metrics_df <- bind_rows(lapply(ps_results, function(r) data.frame(
  method       = "Power-Scaling",
  beta         = r$beta,
  w_real_mean  = r$w_real_mean,
  w_real_sd    = r$w_real_sd,
  sigma_mean   = r$sigma_mean,
  sigma_sd     = r$sigma_sd,
  elpd_per_obs = r$elpd_per_obs,
  rmse         = r$rmse
)))

koh_row <- data.frame(
  method       = "KOH",
  beta         = NA_real_,
  w_real_mean  = mean(koh_w_vec),
  w_real_sd    = sd(koh_w_vec),
  sigma_mean   = mean(koh_sigma_vec),
  sigma_sd     = sd(koh_sigma_vec),
  elpd_per_obs = koh_metrics$elpd_per_obs,
  rmse         = koh_metrics$rmse
)

pw_metrics_df <- bind_rows(lapply(pw_results, function(r) data.frame(
  method       = "Post-Pred Weighting",
  beta         = r$beta,
  w_real_mean  = NA_real_,
  w_real_sd    = NA_real_,
  sigma_mean   = NA_real_,
  sigma_sd     = NA_real_,
  elpd_per_obs = r$elpd_per_obs,
  rmse         = r$rmse
)))

all_metrics_df <- bind_rows(ps_metrics_df, koh_row, pw_metrics_df)
cat("\n──── Results summary ────\n")
print(all_metrics_df)
saveRDS(all_metrics_df, file.path(results_dir, "all_metrics.rds"))
write.csv(all_metrics_df, file.path(results_dir, "all_metrics.csv"),
          row.names = FALSE)

# ── Plots ──────────────────────────────────────────────────────────────────

method_colors <- c("Power-Scaling"       = "#66CCEE",
                   "KOH"                 = "#228833",
                   "Post-Pred Weighting" = "#EE6677")
theme_set(theme_bw(base_size = 13))

# ── ELPD and RMSE ────────────────────────────────────────────────────────

koh_hline_df <- data.frame(method = "KOH",
                           elpd   = koh_metrics$elpd_per_obs,
                           rmse   = koh_metrics$rmse)

linetype_scale <- c("Power-Scaling"       = "solid",
                    "Post-Pred Weighting" = "solid",
                    "KOH"                 = "dashed")

p_elpd <- ggplot(all_metrics_df,
                 aes(x = beta, y = elpd_per_obs,
                     color = method, group = method)) +
  geom_hline(
    data = koh_hline_df,
    aes(yintercept = elpd, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_line(aes(linetype = method), linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, name = NULL) +
  scale_linetype_manual(values = linetype_scale, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(x = TeX("$\\beta$"), y = "ELPD")

p_rmse <- ggplot(all_metrics_df,
                 aes(x = beta, y = rmse,
                     color = method, group = method)) +
  geom_hline(
    data = koh_hline_df,
    aes(yintercept = rmse, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_line(aes(linetype = method), linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, name = NULL) +
  scale_linetype_manual(values = linetype_scale, name = NULL) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(x = TeX("$\\beta$"), y = "RMSE")

p_metrics <- (p_elpd | p_rmse) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "gp_elpd_rmse.pdf"),
       plot = p_metrics, height = 4, width = 9)
ggsave(file.path(results_dir, "gp_elpd_rmse.png"),
       plot = p_metrics, height = 4, width = 9)

# ── Posterior predictive spaghetti ───────────────────────────────────────
pp_base_size  <- 18
pp_strip_size <- 20

theme_pp <- theme_bw(base_size = pp_base_size) %+replace%
  theme(strip.text = element_text(size = pp_strip_size))

betas_sel <- c(0, 0.1, 0.5, 1)

x_pred_orig  <- df_real_oos$w1
df_real_plot <- df_real     %>% transmute(x = w1, y = y_noisy)
df_test_plot <- df_real_oos %>% transmute(x = w1, y = y_noisy)

make_spaghetti_df <- function(f_mat, w_draws, method_label, beta_label) {
  idx <- sample(nrow(f_mat), min(pp_ndraws, nrow(f_mat)))
  bind_rows(lapply(seq_along(idx), function(k) data.frame(
    x       = x_pred_orig,
    mu_pred = f_mat[idx[k], ] * y_sim_sd + y_sim_mu,
    w_draw  = w_draws[idx[k]] * (w_max - w_min) + w_min,
    draw_id = k, beta = beta_label, method = method_label
  )))
}

beta_levels <- sprintf("beta==%.2f", betas_sel)

pspw_df <- bind_rows(
  bind_rows(lapply(betas_sel, function(b)
    make_spaghetti_df(ps_results[[as.character(b)]]$epost$f_pred_matrix,
                      ps_results[[as.character(b)]]$epost$w_real_draws,
                      "Power Scaling", sprintf("beta==%.2f", b)))),
  bind_rows(lapply(betas_sel, function(b)
    make_spaghetti_df(pw_results[[as.character(b)]]$f_pred_matrix,
                      pw_results[[as.character(b)]]$w_real_draws,
                      "Post Pred Weighting", sprintf("beta==%.2f", b))))
) %>% mutate(
  beta   = factor(beta,   levels = beta_levels),
  method = factor(method, levels = c("Power Scaling", "Post Pred Weighting"))
)

koh_pp_df <- make_spaghetti_df(koh_f_mat, koh_w_vec, "KOH", "KOH")

w_lims_plot <- range(c(pspw_df$w_draw, koh_pp_df$w_draw), na.rm = TRUE)

data_points <- list(
  geom_point(data = df_real_plot, aes(x = x, y = y),
             color = "blue", size = 1.5, inherit.aes = FALSE),
  geom_point(data = df_test_plot, aes(x = x, y = y),
             color = "lightblue", size = 1, inherit.aes = FALSE)
)

spaghetti_aes <- list(
  geom_line(alpha = 0.15),
  data_points,
  scale_color_viridis_c(TeX("$\\omega_R$"), limits = w_lims_plot,
                        option = "F", na.value = "grey60"),
  coord_cartesian(ylim = c(-0.5, 10)),
  labs(x = TeX("$x$"), y = TeX("$\\mu_R$")),
  theme(strip.text = element_text(size = 18))
)

p_pspw <- ggplot(pspw_df, aes(x, mu_pred, color = w_draw,
                              group = interaction(draw_id, beta, method))) +
  spaghetti_aes +
  facet_grid(rows = vars(method), cols = vars(beta),
             labeller = labeller(beta = label_parsed, method = label_value)) +
  theme(legend.position = "none")

p_koh <- ggplot(koh_pp_df, aes(x, mu_pred, color = w_draw, group = draw_id)) +
  spaghetti_aes +
  facet_wrap(vars(beta)) +
  guides(color = "none") +
  labs(y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_pp <- (p_pspw | wrap_plots(plot_spacer(), p_koh, plot_spacer(),
                             ncol = 1, heights = c(0.25, 1, 0.25))) +
  plot_layout(widths = c(length(betas_sel), 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "gp_posterior_pred.pdf"),
       plot = p_pp, height = 8, width = 4 * (length(betas_sel) + 1))
ggsave(file.path(results_dir, "gp_posterior_pred.png"),
       plot = p_pp, height = 8, width = 4 * (length(betas_sel) + 1))

# ── Posterior predictive CI plot ─────────────────────────────────────────

x_pred_orig  <- df_real_oos$w1
df_real_plot <- df_real     %>% transmute(x = w1, y = y_noisy)
df_test_plot <- df_real_oos %>% transmute(x = w1, y = y_noisy)

make_ci_df <- function(f_mat, method_label, beta_label,
                       probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  f_orig <- f_mat * y_sim_sd + y_sim_mu
  tibble(
    x      = x_pred_orig,
    lo95   = apply(f_orig, 2, quantile, probs[1]),
    lo50   = apply(f_orig, 2, quantile, probs[2]),
    median = apply(f_orig, 2, quantile, probs[3]),
    hi50   = apply(f_orig, 2, quantile, probs[4]),
    hi95   = apply(f_orig, 2, quantile, probs[5]),
    beta   = beta_label,
    method = method_label
  )
}

pspw_ci_df <- bind_rows(
  bind_rows(lapply(betas_sel, function(b)
    make_ci_df(ps_results[[as.character(b)]]$epost$y_pred_matrix,
               "Power Scaling", sprintf("beta==%.2f", b)))),
  bind_rows(lapply(betas_sel, function(b)
    make_ci_df(pw_results[[as.character(b)]]$y_pred_matrix,
               "Post Pred Weighting", sprintf("beta==%.2f", b))))
) %>% mutate(
  beta   = factor(beta,   levels = beta_levels),
  method = factor(method, levels = c("Power Scaling", "Post Pred Weighting"))
)

koh_ci_df <- make_ci_df(koh_y_mat, "KOH", "KOH")

ci_aes <- list(
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, fill = "red"),
  geom_ribbon(aes(ymin = lo50, ymax = hi50), alpha = 0.25, fill = "red"),
  geom_line(aes(y = median), color = "red", linewidth = 0.8),
  geom_point(data = df_real_plot, aes(x = x, y = y),
             color = "blue", size = 1.5, inherit.aes = FALSE),
  geom_point(data = df_test_plot, aes(x = x, y = y),
             color = "lightblue", size = 1, inherit.aes = FALSE),
  coord_cartesian(ylim = c(-0.5, 10)),
  labs(x = TeX("$x$"), y = TeX("$y_R$")),
  theme(legend.position = "none", strip.text = element_text(size = 18))
)

p_pspw_ci <- ggplot(pspw_ci_df, aes(x = x)) +
  ci_aes +
  facet_grid(rows = vars(method), cols = vars(beta),
             labeller = labeller(beta = label_parsed, method = label_value))

p_koh_ci <- ggplot(koh_ci_df, aes(x = x)) +
  ci_aes +
  facet_wrap(vars(beta)) +
  labs(y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_pp_ci <- (p_pspw_ci | wrap_plots(plot_spacer(), p_koh_ci, plot_spacer(),
                                   ncol = 1, heights = c(0.25, 1, 0.25))) +
  plot_layout(widths = c(length(betas_sel), 1))

ggsave(file.path(results_dir, "gp_posterior_pred_ci.pdf"),
       plot = p_pp_ci, height = 8, width = 4 * (length(betas_sel) + 1))
ggsave(file.path(results_dir, "gp_posterior_pred_ci.png"),
       plot = p_pp_ci, height = 8, width = 4 * (length(betas_sel) + 1))

# ── Posterior densities of w_real and sigma ───────────────────────────────

ps_ws_df <- bind_rows(lapply(betas, function(beta) {
  ep <- ps_results[[as.character(beta)]]$epost
  data.frame(
    beta   = factor(beta),
    w_real = ep$w_real_draws * (w_max - w_min) + w_min,
    sigma  = ep$sigma_draws  * y_sim_sd
  )
}))
koh_ws_df <- data.frame(
  beta   = factor(NA),
  w_real = koh_w_vec   * (w_max - w_min) + w_min,
  sigma  = koh_sigma_vec * y_sim_sd
)

p_w <- ggplot(ps_ws_df, aes(x = w_real, color = beta)) +
  geom_density(linewidth = 0.7) +
  geom_density(data = koh_ws_df, aes(x = w_real),
               color = method_colors["KOH"], linetype = "dashed",
               linewidth = 0.9, inherit.aes = FALSE) +
  annotate("text", x = Inf, y = Inf,
           label = "dashed = KOH", hjust = 1.05, vjust = 1.5,
           color = method_colors["KOH"], size = 3) +
  scale_color_viridis_d(TeX("$\\beta$")) +
  facet_wrap(vars(beta), labeller = label_bquote(beta == .(beta))) +
  labs(x = TeX("$\\omega_R$ (orig. scale)"), y = "density",
       title = "Posterior of calibration parameter") +
  theme(legend.position = "none")

p_s <- ggplot(ps_ws_df, aes(x = sigma, color = beta)) +
  geom_density(linewidth = 0.7) +
  geom_density(data = koh_ws_df, aes(x = sigma),
               color = method_colors["KOH"], linetype = "dashed",
               linewidth = 0.9, inherit.aes = FALSE) +
  scale_color_viridis_d(TeX("$\\beta$")) +
  facet_wrap(vars(beta), labeller = label_bquote(beta == .(beta))) +
  labs(x = TeX("$\\sigma_R$ (orig. scale)"), y = "density",
       title = "Posterior of noise parameter") +
  theme(legend.position = "none")

ggsave(file.path(results_dir, "gp_w_sigma_posterior.pdf"),
       plot = p_w / p_s, height = 10, width = 10)
ggsave(file.path(results_dir, "gp_w_sigma_posterior.png"),
       plot = p_w / p_s, height = 10, width = 10)

# ── RMSE ratio vs KOH ────────────────────────────────────────────────────

ratio_df <- bind_rows(ps_metrics_df, pw_metrics_df) %>%
  filter(beta > 0) %>%
  filter(method == "Power-Scaling") %>%
  mutate(rmse_ratio = rmse / koh_metrics$rmse)

  
best_rmse_ratio <- min(ratio_df$rmse_ratio)

p_rmse_ratio <- ggplot(ratio_df, aes(x = beta, y = rmse_ratio,
                                      color = method, group = method)) +
  geom_hline(yintercept = 1, color = method_colors["KOH"],
             linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(breaks = sort(c(scales::pretty_breaks()(ratio_df$rmse_ratio),
                                     best_rmse_ratio)),
                     labels = function(x) ifelse(x == best_rmse_ratio,
                                                 round(best_rmse_ratio, 2), x)) +
  labs(x = TeX("$\\beta$"), y = "RMSE / RMSE(KOH)") +
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "gp_rmse_ratio.pdf"), plot = p_rmse_ratio,
       height = 4, width = 5)
ggsave(file.path(results_dir, "gp_rmse_ratio.png"), plot = p_rmse_ratio,
       height = 4, width = 5)

cat("\nAll done. Results in:", results_dir, "\n")
