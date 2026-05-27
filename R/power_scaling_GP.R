# power_scaling_GP.R
#
# Implements the two-step power-scaling procedure for GP surrogates,
# mirroring the PCE power-scaling logic in main.R / eval.R, but using
# GP_power_scaling_step1.stan and GP_power_scaling_step2.stan.
#
# For each beta in {0, 0.1, 0.5, 1}:
#   Step 1 – joint GP training with power-scaled sim + real likelihoods
#   Step 2 – E-Post: for each Step-1 draw of (rho, alpha_gp, f_sim),
#             run a short inner MCMC to refine (w_real, sigma) on real
#             data only (full likelihood).
#
# Outputs (saved to results/gp_power_scaling/):
#   • Stan fit objects for both steps (one per beta)
#   • MCMC diagnostics (trace, Rhat plots)
#   • Posterior predictive spaghetti plot coloured by w_real draw
#   • ELPD and RMSE table across betas
# ---------------------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(tidyverse)
library(viridis)
library(latex2exp)
library(matrixStats)
library(loo)

set.seed(42)

# ── 0. User settings ─────────────────────────────────────────────────────────

# Paths to the two Stan models (adjust if needed)
stan_step1_file <- "stan_code/GP_power_scaling_step1.stan"
stan_step2_file <- "stan_code/GP_power_scaling_step2.stan"

# Data files produced by the existing PCE pipeline
data_dir     <- "results/log_trend_log_sin_N_sim-100_likelihood-normal_link-identity_ylogtransform-FALSE_real_train_upper_0.5 _koh"
results_dir  <- file.path("results", "gp_power_scaling")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Betas to evaluate
betas <- c(0.5, 1)

# MCMC settings
chains_1          <- 4          # Step-1 chains
iter_warmup_1     <- 1000
iter_sampling_1   <- 1000
parallel_chains_1 <- 4

# Step-2 inner MCMC (E-Post): run once per outer draw but keep it cheap
# We subsample S1_sub draws from Step 1 to keep Step 2 tractable.
S1_sub            <- 200        # number of Step-1 draws used as Step-2 seeds
inner_chains      <- 1          # single chain per outer draw
inner_iter_warmup <- 200
inner_iter_sample <- 50         # we only keep the last draw

adapt_delta       <- 0.95

# Prediction / evaluation settings
pp_ndraws <- 200   # spaghetti lines in the posterior predictive plot

# ── 1. Load data ─────────────────────────────────────────────────────────────
# Adjust these readRDS calls / file names to match your actual data objects.
# The script assumes the same data frame structure as the PCE pipeline:
#   df_sim$w1   – observable input x  (simulator)
#   df_sim$w2   – calibration input w  (simulator, known)
#   df_sim$y_noisy – noisy output       (simulator)
#   df_real$w1  – observable input x   (real)
#   df_real$y_noisy – noisy output      (real)

df_sim      <- readRDS(file.path(data_dir, "log_trend_simulation.Rda"))
df_real     <- readRDS(file.path(data_dir, "log_sin_real_x1_slice.Rda"))
df_real_oos <- readRDS(file.path(data_dir, "log_sin_real_oos_ood_x1_uniform.Rda"))

# Optional: subsample for faster iteration during development
set.seed(123)
df_sim      <- df_sim[sample(nrow(df_sim),  30), ]
df_real     <- df_real[sample(nrow(df_real), 10), ]
df_real_oos <- df_real_oos[sample(nrow(df_real_oos), 50), ]

cat(sprintf("Data loaded: N_sim=%d  N_real=%d  N_pred=%d\n",
            nrow(df_sim), nrow(df_real), nrow(df_real_oos)))

# ── 2. Pre-process ────────────────────────────────────────────────────────────

p <- 1   # dimension of x (observable input)
q <- 1   # dimension of w (calibration parameter)

y_sim  <- df_sim$y_noisy
x_sim  <- matrix(df_sim$w1,  ncol = p)
w_sim  <- matrix(df_sim$w2,  ncol = q)
y_real <- df_real$y_noisy
x_real <- matrix(df_real$w1, ncol = p)
x_pred <- matrix(df_real_oos$w1, ncol = p)

# --- Standardise outputs (z-score using sim mean/sd) ---
y_sim_mu <- mean(y_sim, na.rm = TRUE)
y_sim_sd <- sd(y_sim,   na.rm = TRUE)
y_sim_scaled  <- (y_sim  - y_sim_mu) / y_sim_sd
y_real_scaled <- (y_real - y_sim_mu) / y_sim_sd

# --- Scale observable inputs to [0, 1] ---
x_all <- rbind(x_real, x_sim)
x_min <- min(x_all); x_max <- max(x_all)
scale_x <- function(x) (x - x_min) / (x_max - x_min)
x_sim_sc  <- scale_x(x_sim)
x_real_sc <- scale_x(x_real)
x_pred_sc <- scale_x(x_pred)

# --- Scale calibration parameter to [0, 1] ---
w_min <- min(w_sim); w_max <- max(w_sim)
scale_w <- function(w) (w - w_min) / (w_max - w_min)
w_sim_sc <- scale_w(w_sim)

# Prior for w_real: uniform over [0,1] approximated by N(0.5, 0.5) truncated.
# Adjust to match domain knowledge.
w_prior_mean  <- 0.5
w_prior_sigma <- 0.5

# Helper: convert beta -> (alpha_sim, alpha_real) per Eq. 21
get_alphas_gp <- function(beta) {
  if (beta < 0.5) {
    alpha_sim  <- beta / (1 - beta)
    alpha_real <- 1.0
  } else {
    alpha_sim  <- 1.0
    alpha_real <- (1 - beta) / beta
  }
  # Edge cases
  if (beta == 0) { alpha_sim  <- 0.0; alpha_real <- 1.0 }
  if (beta == 1) { alpha_sim  <- 1.0; alpha_real <- 0.0 }
  list(alpha_sim = alpha_sim, alpha_real = alpha_real)
}

# Build arrays of vectors for Stan (Stan expects array[] vector[d])
to_vec_array <- function(mat) {
  lapply(seq_len(nrow(mat)), function(i) mat[i, ])
}

x_sim_arr  <- to_vec_array(x_sim_sc)
w_sim_arr  <- to_vec_array(w_sim_sc)
x_real_arr <- to_vec_array(x_real_sc)
x_pred_arr <- to_vec_array(x_pred_sc)

# ── 3. Compile Stan models ────────────────────────────────────────────────────

cat("Compiling Step-1 model ...\n")
mod_step1 <- cmdstan_model(stan_step1_file)

cat("Compiling Step-2 model ...\n")
mod_step2 <- cmdstan_model(stan_step2_file)

# ── 4. Helper: collect Step-2 draws (E-Post) ─────────────────────────────────

#' Run Step-2 (E-Post inner MCMC) for each of S1_sub draws from Step 1.
#'
#' Returns a list with:
#'   w_real_draws  – numeric vector of length S1_sub
#'   sigma_draws   – numeric vector of length S1_sub
#'   y_pred_matrix – (S1_sub x N_pred) matrix of posterior predictive draws
#'   f_pred_matrix – (S1_sub x N_pred) matrix of latent GP draws
run_epost_step2 <- function(fit_step1, stan_base, x_sim_arr, w_sim_arr,
                             x_real_arr, x_pred_arr,
                             y_real_scaled, N_sim, N_real, N_pred, p, q,
                             w_prior_mean, w_prior_sigma,
                             S1_sub, inner_iter_warmup, inner_iter_sample,
                             inner_chains) {

  draws_s1 <- as_draws_df(fit_step1$draws(
    c("rho", "alpha_gp", "sigma", "w_real", paste0("f_sim[", 1:N_sim, "]"))
  ))

  # Subsample S1_sub rows uniformly
  idx <- round(seq(1, nrow(draws_s1), length.out = S1_sub))
  draws_sub <- draws_s1[idx, ]

  w_real_draws  <- numeric(S1_sub)
  sigma_draws   <- numeric(S1_sub)
  y_pred_matrix <- matrix(NA_real_, nrow = S1_sub, ncol = N_pred)
  f_pred_matrix <- matrix(NA_real_, nrow = S1_sub, ncol = N_pred)

  rho_cols   <- grep("^rho\\[", names(draws_sub), value = TRUE)
  fsim_cols  <- grep("^f_sim\\[", names(draws_sub), value = TRUE)

  cat(sprintf("  Running E-Post Step 2 (%d inner fits) ...\n", S1_sub))
  pb <- txtProgressBar(min = 0, max = S1_sub, style = 3)

  for (s in seq_len(S1_sub)) {
    row   <- draws_sub[s, ]
    rho_s <- as.numeric(row[rho_cols])
    alpha_gp_s <- as.numeric(row[["alpha_gp"]])
    f_sim_s    <- as.numeric(row[fsim_cols])

    # Use Step-1 estimates of w_real and sigma as init for inner chain
    w_init     <- as.numeric(row[["w_real"]])
    sigma_init <- as.numeric(row[["sigma"]])

    stan_data_s2 <- c(
      stan_base,
      list(
        f_sim      = f_sim_s,
        rho        = rho_s,
        alpha_gp   = alpha_gp_s,
        w_prior_mean  = w_prior_mean,
        w_prior_sigma = w_prior_sigma
      )
    )

    init_s2 <- list(list(
      w_real = max(0.01, min(0.99, w_init)),
      sigma  = max(0.01, sigma_init)
    ))

    suppressMessages({
      fit_s2 <- mod_step2$sample(
        data            = stan_data_s2,
        init            = init_s2,
        chains          = inner_chains,
        iter_warmup     = inner_iter_warmup,
        iter_sampling   = inner_iter_sample,
        refresh         = 0,
        show_messages   = FALSE
      )
    })

    draws_s2 <- as_draws_df(fit_s2$draws())

    # Keep only the final draw (last iteration)
    last <- nrow(draws_s2)
    w_real_draws[s]  <- as.numeric(draws_s2[last, "w_real"])
    sigma_draws[s]   <- as.numeric(draws_s2[last, "sigma"])

    y_pred_cols <- grep("^y_pred\\[", names(draws_s2), value = TRUE)
    f_pred_cols <- grep("^f_pred\\[", names(draws_s2), value = TRUE)
    y_pred_matrix[s, ] <- as.numeric(draws_s2[last, y_pred_cols])
    f_pred_matrix[s, ] <- as.numeric(draws_s2[last, f_pred_cols])

    setTxtProgressBar(pb, s)
  }
  close(pb)

  list(
    w_real_draws  = w_real_draws,
    sigma_draws   = sigma_draws,
    y_pred_matrix = y_pred_matrix,
    f_pred_matrix = f_pred_matrix
  )
}

# ── 5. Helper: compute ELPD and RMSE ─────────────────────────────────────────

#' Compute ELPD (log score) and RMSE against test data.
#'
#' @param f_pred_matrix  (S x N_test) matrix of latent GP predictions
#' @param sigma_draws    length-S vector of sigma draws
#' @param y_test         length-N_test vector of (scaled) observed values
#' @return list(elpd_per_obs, rmse)
compute_metrics <- function(f_pred_matrix, sigma_draws, y_test) {
  S      <- nrow(f_pred_matrix)
  N_test <- ncol(f_pred_matrix)
  stopifnot(length(sigma_draws) == S, length(y_test) == N_test)

  # Log-likelihood matrix (S x N_test)
  ll_matrix <- matrix(NA_real_, S, N_test)
  for (s in seq_len(S)) {
    ll_matrix[s, ] <- dnorm(y_test,
                            mean = f_pred_matrix[s, ],
                            sd   = sigma_draws[s],
                            log  = TRUE)
  }

  # ELPD: log of Monte Carlo average over S draws, summed over observations
  log_mean_lik <- matrixStats::colLogSumExps(ll_matrix) - log(S)
  elpd_per_obs <- mean(log_mean_lik)

  # RMSE based on posterior mean prediction
  mu_pred <- colMeans(f_pred_matrix)  # E[f | data], shape N_test
  y_mat   <- matrix(y_test, nrow = S, ncol = N_test, byrow = TRUE)
  rmse    <- sqrt(mean((f_pred_matrix - y_mat)^2))

  list(elpd_per_obs = elpd_per_obs, rmse = rmse)
}

# ── 6. Main loop over betas ───────────────────────────────────────────────────

# Fixed part of the Step-2 Stan data (same for all betas)
stan_base_step2 <- list(
  N_sim    = nrow(df_sim),
  N_real   = nrow(df_real),
  N_pred   = nrow(df_real_oos),
  p        = p,
  q        = q,
  x_sim    = x_sim_arr,
  w_sim    = w_sim_arr,
  x_real   = x_real_arr,
  y_real   = y_real_scaled,
  x_pred   = x_pred_arr
)

results_list <- list()  # will hold one entry per beta

for (beta in betas) {
  cat(sprintf("\n══════════════════════════════════════\n"))
  cat(sprintf("  beta = %.2f\n", beta))
  cat(sprintf("══════════════════════════════════════\n"))

  alphas     <- get_alphas_gp(beta)
  alpha_sim  <- alphas$alpha_sim
  alpha_real <- alphas$alpha_real
  cat(sprintf("  alpha_sim=%.4f  alpha_real=%.4f\n", alpha_sim, alpha_real))

  # --- Step 1 Stan data ---
  stan_data_s1 <- list(
    N_sim         = nrow(df_sim),
    N_real        = nrow(df_real),
    p             = p,
    q             = q,
    x_sim         = x_sim_arr,
    w_sim         = w_sim_arr,
    y_sim         = y_sim_scaled,
    x_real        = x_real_arr,
    y_real        = y_real_scaled,
    alpha_sim     = alpha_sim,
    alpha_real    = alpha_real,
    w_prior_mean  = w_prior_mean,
    w_prior_sigma = w_prior_sigma
  )

  # --- Step 1: joint training ---
  cat("  Step 1: joint GP training ...\n")
  fit_s1 <- mod_step1$sample(
    data            = stan_data_s1,
    seed            = 42,
    chains          = chains_1,
    parallel_chains = parallel_chains_1,
    iter_warmup     = iter_warmup_1,
    iter_sampling   = iter_sampling_1,
    adapt_delta     = adapt_delta,
    refresh         = 200
  )

  # Save Step-1 fit
  fit_s1_path <- file.path(results_dir,
    sprintf("gp_ps_step1_beta_%.2f", beta))
  fit_s1$save_object(fit_s1_path)
  cat(sprintf("  Step-1 fit saved -> %s\n", fit_s1_path))

  # MCMC diagnostics for Step 1
  np_s1 <- nuts_params(fit_s1)
  p_trace <- mcmc_trace(fit_s1$draws(),
                        pars  = c("alpha_gp", "sigma", "w_real"),
                        regex_pars = "^rho",
                        np    = np_s1)
  ggsave(file.path(results_dir,
    sprintf("trace_step1_beta_%.2f.png", beta)),
    plot = p_trace, width = 10, height = 8)

  rhats_s1 <- rhat(fit_s1)
  p_rhat <- mcmc_rhat(rhats_s1[grepl("alpha_gp|sigma|w_real|rho",
                                      names(rhats_s1))]) +
    yaxis_text(hjust = 1)
  ggsave(file.path(results_dir,
    sprintf("rhat_step1_beta_%.2f.png", beta)),
    plot = p_rhat, width = 8, height = 6)

  cat("  Step-1 summary (key params):\n")
  print(fit_s1$summary(c("alpha_gp", "sigma", "w_real")))

  # --- Step 2: E-Post refinement of (w_real, sigma) ---
  epost <- run_epost_step2(
    fit_step1        = fit_s1,
    stan_base        = stan_base_step2,
    x_sim_arr        = x_sim_arr,
    w_sim_arr        = w_sim_arr,
    x_real_arr       = x_real_arr,
    x_pred_arr       = x_pred_arr,
    y_real_scaled    = y_real_scaled,
    N_sim            = nrow(df_sim),
    N_real           = nrow(df_real),
    N_pred           = nrow(df_real_oos),
    p                = p,
    q                = q,
    w_prior_mean     = w_prior_mean,
    w_prior_sigma    = w_prior_sigma,
    S1_sub           = S1_sub,
    inner_iter_warmup = inner_iter_warmup,
    inner_iter_sample = inner_iter_sample,
    inner_chains      = inner_chains
  )

  # Save E-Post results
  epost_path <- file.path(results_dir,
    sprintf("gp_ps_epost_beta_%.2f.rds", beta))
  saveRDS(epost, epost_path)
  cat(sprintf("  E-Post results saved -> %s\n", epost_path))

  # Diagnostics: posterior of w_real and sigma after Step 2
  cat(sprintf("  w_real (Step2): mean=%.3f  sd=%.3f\n",
              mean(epost$w_real_draws), sd(epost$w_real_draws)))
  cat(sprintf("  sigma  (Step2): mean=%.3f  sd=%.3f\n",
              mean(epost$sigma_draws),  sd(epost$sigma_draws)))

  # --- Evaluation on test set (OOS) ---
  y_test_scaled <- (df_real_oos$y_noisy - y_sim_mu) / y_sim_sd
  metrics <- compute_metrics(
    f_pred_matrix = epost$f_pred_matrix,
    sigma_draws   = epost$sigma_draws,
    y_test        = y_test_scaled
  )

  cat(sprintf("  ELPD/obs = %.4f   RMSE = %.4f\n",
              metrics$elpd_per_obs, metrics$rmse))

  results_list[[as.character(beta)]] <- list(
    beta          = beta,
    alpha_sim     = alpha_sim,
    alpha_real    = alpha_real,
    w_real_mean   = mean(epost$w_real_draws),
    w_real_sd     = sd(epost$w_real_draws),
    sigma_mean    = mean(epost$sigma_draws),
    sigma_sd      = sd(epost$sigma_draws),
    elpd_per_obs  = metrics$elpd_per_obs,
    rmse          = metrics$rmse,
    epost         = epost
  )
}

# ── 7. Summary table ──────────────────────────────────────────────────────────

metrics_df <- bind_rows(lapply(results_list, function(r) {
  data.frame(
    beta         = r$beta,
    alpha_sim    = r$alpha_sim,
    alpha_real   = r$alpha_real,
    w_real_mean  = r$w_real_mean,
    w_real_sd    = r$w_real_sd,
    sigma_mean   = r$sigma_mean,
    sigma_sd     = r$sigma_sd,
    elpd_per_obs = r$elpd_per_obs,
    rmse         = r$rmse
  )
}))

cat("\n══════════════════════════════════════\n")
cat("  Results summary\n")
cat("══════════════════════════════════════\n")
print(metrics_df)

saveRDS(metrics_df, file.path(results_dir, "gp_ps_metrics.rds"))
write.csv(metrics_df, file.path(results_dir, "gp_ps_metrics.csv"),
          row.names = FALSE)

# ── 8. Plots ──────────────────────────────────────────────────────────────────

# 8a. ELPD and RMSE over beta
p_elpd <- ggplot(metrics_df, aes(x = beta, y = elpd_per_obs)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = TeX("$\\beta$"), y = "ELPD / obs",
       title = "GP Power-Scaling: Predictive Performance") +
  theme_bw(base_size = 14)

p_rmse <- ggplot(metrics_df, aes(x = beta, y = rmse)) +
  geom_point(size = 3) +
  geom_line() +
  scale_y_log10() +
  labs(x = TeX("$\\beta$"), y = "RMSE (log scale)") +
  theme_bw(base_size = 14)

library(patchwork)
p_perf <- p_elpd / p_rmse
ggsave(file.path(results_dir, "gp_ps_elpd_rmse.pdf"),
       plot = p_perf, height = 8, width = 6)
ggsave(file.path(results_dir, "gp_ps_elpd_rmse.png"),
       plot = p_perf, height = 8, width = 6)

# 8b. Posterior predictive spaghetti plots (one panel per beta)
pp_list <- list()
for (beta in betas) {
  r <- results_list[[as.character(beta)]]
  ep <- r$epost

  # Subsample pp_ndraws from available draws for the spaghetti
  n_avail  <- nrow(ep$y_pred_matrix)
  idx_pp   <- sample(n_avail, min(pp_ndraws, n_avail))

  # x_pred on original scale for plotting
  x_pred_orig <- df_real_oos$w1

  pp_df <- do.call(rbind, lapply(seq_along(idx_pp), function(k) {
    s <- idx_pp[k]
    data.frame(
      x_plot   = x_pred_orig,
      y_pred   = ep$y_pred_matrix[s, ] * y_sim_sd + y_sim_mu,
      mu_pred  = ep$f_pred_matrix[s, ] * y_sim_sd + y_sim_mu,
      w_draw   = ep$w_real_draws[s],
      draw_id  = k,
      beta     = beta
    )
  }))
  pp_list[[as.character(beta)]] <- pp_df
}

pp_all <- bind_rows(pp_list)

# w_real draws back to original scale for colour axis
pp_all <- pp_all %>%
  mutate(w_real_orig = w_draw * (w_max - w_min) + w_min)

p_pp <- ggplot() +
  geom_line(
    data  = pp_all,
    aes(x = x_plot, y = mu_pred,
        color = w_real_orig, group = interaction(draw_id, beta)),
    alpha = 0.15
  ) +
  geom_point(
    data  = df_real,
    aes(x = w1, y = y_noisy),
    color = "blue", size = 1.5
  ) +
  geom_point(
    data  = df_real_oos,
    aes(x = w1, y = y_noisy),
    color = "lightblue", size = 1
  ) +
  facet_wrap(vars(beta),
             labeller = label_bquote(beta == .(beta)),
             nrow = 1) +
  scale_color_viridis_c(TeX("$\\omega_R$ (orig. scale)"), option = "F") +
  labs(x = TeX("$x$"), y = TeX("$\\hat{y}_R$ (posterior mean)"),
       title = "GP Power-Scaling: Posterior Predictive") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(results_dir, "gp_ps_posterior_pred.pdf"),
       plot = p_pp, height = 5, width = 14)
ggsave(file.path(results_dir, "gp_ps_posterior_pred.png"),
       plot = p_pp, height = 5, width = 14)

# 8c. Posterior of w_real and sigma over beta
w_sigma_df <- bind_rows(lapply(betas, function(beta) {
  r <- results_list[[as.character(beta)]]
  ep <- r$epost
  data.frame(
    beta   = beta,
    w_real = ep$w_real_draws * (w_max - w_min) + w_min,  # back to original scale
    sigma  = ep$sigma_draws   * y_sim_sd                 # back to original scale
  )
}))

p_w <- ggplot(w_sigma_df, aes(x = w_real, color = factor(beta))) +
  geom_density() +
  scale_color_viridis_d(TeX("$\\beta$")) +
  labs(x = TeX("$\\omega_R$ (orig. scale)"), y = "density",
       title = "Posterior of calibration parameter") +
  facet_wrap(vars(beta), labeller = label_bquote(beta == .(beta))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

p_s <- ggplot(w_sigma_df, aes(x = sigma, color = factor(beta))) +
  geom_density() +
  scale_color_viridis_d(TeX("$\\beta$")) +
  labs(x = TeX("$\\sigma_R$ (orig. scale)"), y = "density",
       title = "Posterior of noise parameter") +
  facet_wrap(vars(beta), labeller = label_bquote(beta == .(beta))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

p_ws <- p_w / p_s
ggsave(file.path(results_dir, "gp_ps_w_sigma_posterior.pdf"),
       plot = p_ws, height = 10, width = 10)
ggsave(file.path(results_dir, "gp_ps_w_sigma_posterior.png"),
       plot = p_ws, height = 10, width = 10)

cat("\nAll done. Results in:", results_dir, "\n")
