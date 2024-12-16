library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(stringr)
library(loo)
library(forcats)
library(truncnorm)
library(purrr)

source(file.path(here(),"R/get_data.R"))
source(file.path(here(),"R/plot.R"))
source(file.path(here(),"R/utils.R"))
source(file.path(here(), "R/load_config_args.R"))


results_path <- file.path("results", experiment_name)

set.seed(seed)
# read data
# variable changing between 1d/2d
df_sim <- readRDS(file=get_data_save_name(results_path, model_sim, "simulation"))
df_real <- readRDS(file=get_data_save_name(results_path, model_real, "real", sample_x_real))
df_real_ood <- readRDS(file=get_data_save_name(results_path, model_real, "real_ood", sample_x_real_ood))
df_real_ood <- df_real_ood[df_real_ood$y_noisy>0, ]
df_real_oos <- readRDS(file=get_data_save_name(results_path, model_real, "real_oos", sample_x_real_ood))
df_real_oos <- df_real_oos[df_real_oos$y_noisy>0, ]
df_real_oos_ood <- readRDS(file=get_data_save_name(results_path, model_real, "real_oos_ood", sample_x_real_ood))
df_real_oos_ood <- df_real_oos_ood[df_real_oos_ood$y_noisy>0, ]
df_real_oos_ood_2 <- readRDS(file=get_data_save_name(results_path, model_real, "real_oos_ood_2", sample_x_real_ood))
df_real_oos_ood_2 <- df_real_oos_ood_2[df_real_oos_ood_2$y_noisy>0, ]
real_test_list <- list("OOD" = df_real_ood,
                       "OOS" = df_real_oos,
                       "OOS+OOD" = df_real_oos_ood,
                       "OOS+OOD2" = df_real_oos_ood_2)

x_plot <- seq(from=-1, to=1,length.out=100)

pp_list <- list()
surrogate_coeffs_list <- list()
inference_params_list <- list()
test_list <- list()

for (data_integration_scheme in data_integration_schemes){
  for (i in seq_along(betas)){
    # read fit
    beta <- betas[[i]]
    surrogate_model_fit <- readRDS(
      file = get_fit_file_name(results_path, data_integration_scheme, N_real,
                               beta))
    ## store surrogate coefficient in data frame
    surrogate_coeffs_df <- as_draws_df(surrogate_model_fit$draws())
    surrogate_coeffs_df$data_integration_scheme <- data_integration_scheme
    surrogate_coeffs_df$beta <- beta
    surrogate_coeffs_list[[length(surrogate_coeffs_list)+1]] <- surrogate_coeffs_df

    # use additional omega/sigma inference
    surrogate_inference_model_fit <- readRDS(
      file = get_fit_file_name(results_path, paste0(data_integration_scheme, "_inference"), N_real,
                               beta))

    ## store inference parameters in data frame
    w_draws <- get_inner_mcmc_draws(surrogate_inference_model_fit, "w_real",
                                    chains, iter_sampling, inner_chains)
    sigma_draws <- get_inner_mcmc_draws(surrogate_inference_model_fit, "sigma",
                                        chains, iter_sampling, inner_chains)
    inference_params_df <- as_draws_df(w_draws)
    inference_params_df$sigma <- c(sigma_draws)
    if (inner_mcmc_type == "e_post") {
      inference_params_df$rhat <- sapply(
        seq(1, chains*iter_sampling*inner_chains, inner_chains),
        function(j) posterior::rhat(
          surrogate_inference_model_fit$draws("w_real[1]")[, j+(0:(inner_chains-1)), ]))
    }
    inference_params_df$data_integration_scheme <- data_integration_scheme
    inference_params_df$beta <- beta
    inference_params_list[[length(inference_params_list)+1]] <- inference_params_df

    posterior_pred <- calc_post_pred(
      c_draws = as_draws_matrix(surrogate_model_fit$draws("c")),
      c_0_draws = as_draws_matrix(surrogate_model_fit$draws("c_0")),
      w_draws = w_draws,
      sigma_draws = sigma_draws,
      x = x_plot,
      number_draws = pp_ndraws,
      link = surrogate_link,
      likelihood = surrogate_likelihood,
      y_log_transform = y_log_transform)

    posterior_pred$data_integration_scheme <- data_integration_scheme
    posterior_pred$beta <- beta
    pp_list[[length(pp_list)+1]] <- posterior_pred

    # model evaluation on test sets
    for (test_scenario in names(real_test_list)){
      df_real_test <- real_test_list[[test_scenario]]
      linpred_matrix <- calc_pce_linpred_matrix(
         c_draws = as_draws_matrix(surrogate_model_fit$draws("c")),
         c_0_draws = as_draws_matrix(surrogate_model_fit$draws("c_0")),
         w_draws = w_draws,
         sigma_draws = sigma_draws,
         x = df_real_test$w1,
         link = surrogate_link,
         y_log_transform = y_log_transform
      )
      epred_matrix <- calc_pce_epred_matrix(linpred_matrix, sigma_draws,
                                            surrogate_likelihood)
      test_metrics <- get_metrics(linpred_matrix, epred_matrix, sigma_draws,
                                  df_real_test$y, df_real_test$y_noisy,
                                  surrogate_likelihood)
      test_elpd <- test_metrics$elpd$estimates
      test_rmse <- test_metrics$rmse
      test_list[[length(test_list)+1]] <- data.frame(
        "real_test_elpd" = test_elpd[[1]]/nrow(df_real_test),
        "real_test_elpd_se" = test_elpd[[3]]/nrow(df_real_test),
        "real_test_rmse" = test_rmse,
        "beta" = beta,
        "test_scenario" = test_scenario,
        "data_integration_scheme" = data_integration_scheme)
    }
    # diagnostics plots
    surrogate_fit_np <- nuts_params(surrogate_model_fit)
    mcmc_trace(surrogate_model_fit$draws(), regex_pars=c("c", "sigma", "w_real"),
                     np = surrogate_fit_np)
    ggsave(file.path(results_path, sprintf("trace-surrogate_%s-N_real_%s-beta_%.2f.png", data_integration_scheme, N_real, beta)), width=10, height=8)
    surrogate_fit_rhats <- rhat(surrogate_model_fit)
    mcmc_rhat(surrogate_fit_rhats[grepl( "c|sigma|w_real" , names( surrogate_fit_rhats ) ) ]) + yaxis_text(hjust = 1)
    ggsave(file.path(results_path, sprintf("rhat-surrogate_%s-N_real_%s-beta_%.2f.png", data_integration_scheme, N_real, beta)), width=10, height=8)
  }
}

# surrogate posterior predictive stacking
data_integration_scheme <- "pp_weighting"
surrogate_model_fit_real <- readRDS(
  file = get_fit_file_name(results_path, "pp_weighting", N_real,
                           0))
surrogate_model_fit_sim <- readRDS(
  file = get_fit_file_name(results_path, "power_scaling", N_real,
                           1))
surrogate_inference_model_fit_sim <- readRDS(
  file = get_fit_file_name(results_path, "power_scaling_inference", N_real,
                           1))
draws_real <- bind_draws(
  merge_chains(as_draws_array(surrogate_model_fit_real$draws(
    c("c", "c_0", "sigma")))))

draws_sim <- bind_draws(
  merge_chains(as_draws_array(surrogate_model_fit_sim$draws(
    c("c", "c_0")))),
  merge_chains(get_inner_mcmc_draws(surrogate_inference_model_fit_sim, "w_real",
                                    chains, iter_sampling, inner_chains)),
  merge_chains(get_inner_mcmc_draws(surrogate_inference_model_fit_sim, "sigma",
                                    chains, iter_sampling, inner_chains)))

## calculate optimal weights via stacking
optimal_beta_stacking <- data.frame( nrow = 1)
for (test_scenario in names(real_test_list)){
  df_real_test <- real_test_list[[test_scenario]]
  linpred_matrix_real <- calc_pce_linpred_matrix(
    c_draws = as_draws_matrix(subset_draws(draws_real, "c")),
    c_0_draws = as_draws_matrix(subset_draws(draws_real, "c_0")),
    w_draws = NULL,
    sigma_draws = as_draws_matrix(subset_draws(draws_real, "sigma")),
    x = df_real_test$w1,
    link = surrogate_link,
    y_log_transform = y_log_transform)
  epred_matrix_real <- calc_pce_epred_matrix(linpred_matrix_real, as_draws_matrix(subset_draws(draws_real, "sigma")),
                                        surrogate_likelihood)
  lpd_matrix_real <- get_metrics(
    linpred_matrix_real, epred_matrix_real,
    sigma_draws = as_draws_matrix(subset_draws(draws_real, "sigma")),
    y_true = df_real_test$y, y_true_noisy = df_real_test$y_noisy,
    likelihood = surrogate_likelihood)$lpd_matrix
  linpred_matrix_sim <- calc_pce_linpred_matrix(
    c_draws = as_draws_matrix(subset_draws(draws_sim, "c")),
    c_0_draws = as_draws_matrix(subset_draws(draws_sim, "c_0")),
    w_draws = as_draws_matrix(subset_draws(draws_sim, "w_real")),
    sigma_draws = as_draws_matrix(subset_draws(draws_sim, "sigma")),
    x = df_real_test$w1,
    link = surrogate_link,
    y_log_transform = y_log_transform)
  epred_matrix_sim <- calc_pce_epred_matrix(linpred_matrix_sim, as_draws_matrix(subset_draws(draws_sim, "sigma")),
                                             surrogate_likelihood)
  lpd_matrix_sim <- get_metrics(
    linpred_matrix_sim, epred_matrix_sim,
    sigma_draws = as_draws_matrix(subset_draws(draws_sim, "sigma")),
    y_true = df_real_test$y, y_true_noisy = df_real_test$y_noisy,
    likelihood = surrogate_likelihood)$lpd_matrix

  lpd_point_test <- cbind(
    get_point_elpd(lpd_matrix_real),
    get_point_elpd(lpd_matrix_sim)
  )
  stacking_wts <- stacking_weights(lpd_point_test)
  optimal_beta_stacking[, test_scenario] <- stacking_wts[2]
  elpd_df <- data.frame(
    x = df_real_test$w1,
    elpd_real = get_point_elpd(lpd_matrix_real),
    elpd_sim = get_point_elpd(lpd_matrix_sim),
    y = df_real_test$y_noisy
  )
  p1 <- ggplot(elpd_df)+
    geom_point(aes(x = x, y = y, color=elpd_real))+
    xlim(-1, 1)+
    scale_color_viridis_c(limits=NULL, option="F")
  p2 <- ggplot(elpd_df)+
    geom_point(aes(x = x, y = y, color=elpd_sim))+
    xlim(-1, 1)+
    scale_color_viridis_c(limits=NULL, option="F")
  print(p1 / p2)
}

for (i in seq_along(betas)){
  beta <- betas[[i]]
  posterior_pred_real <- calc_post_pred(
    c_draws = as_draws_matrix(subset_draws(draws_real, "c")),
    c_0_draws = as_draws_matrix(subset_draws(draws_real, "c_0")),
    w_draws = NULL,
    sigma_draws = as_draws_matrix(subset_draws(draws_real, "sigma")),
    x = x_plot,
    number_draws = round((1 - beta) * pp_ndraws),
    link = surrogate_link,
    likelihood = surrogate_likelihood,
    y_log_transform = y_log_transform
  )
  posterior_pred_sim <- calc_post_pred(
    c_draws = as_draws_matrix(subset_draws(draws_sim, "c")),
    c_0_draws = as_draws_matrix(subset_draws(draws_sim, "c_0")),
    w_draws = as_draws_matrix(subset_draws(draws_sim, "w_real")),
    sigma_draws =as_draws_matrix(subset_draws(draws_sim, "sigma")),
    x = x_plot,
    number_draws = round(beta * pp_ndraws),
    link = surrogate_link,
    likelihood = surrogate_likelihood,
    y_log_transform = y_log_transform
  )
  posterior_pred_sim$w_id <- posterior_pred_sim$w_id + round((1 - beta) * pp_ndraws)
  posterior_pred_weighted <- rbind(posterior_pred_real, posterior_pred_sim)
  posterior_pred_weighted$data_integration_scheme <- data_integration_scheme
  posterior_pred_weighted$beta <- beta
  pp_list[[length(pp_list)+1]] <- posterior_pred_weighted
  # calculate metrics on different test scenarios

  for (test_scenario in names(real_test_list)){
    df_real_test <- real_test_list[[test_scenario]]
    linpred_matrix_real <- calc_pce_linpred_matrix(
      c_draws = as_draws_matrix(subset_draws(draws_real, "c")),
      c_0_draws = as_draws_matrix(subset_draws(draws_real, "c_0")),
      w_draws = NULL,
      sigma_draws = as_draws_matrix(subset_draws(draws_real, "sigma")),
      x = df_real_test$w1,
      number_draws = round((1 - beta) * nrow(draws_real)),
      link = surrogate_link,
      y_log_transform = y_log_transform)
    epred_matrix_real <- calc_pce_epred_matrix(linpred_matrix_real, as_draws_matrix(subset_draws(draws_real, "sigma")),
                                               surrogate_likelihood)
    linpred_matrix_sim <- calc_pce_linpred_matrix(
      c_draws = as_draws_matrix(subset_draws(draws_sim, "c")),
      c_0_draws = as_draws_matrix(subset_draws(draws_sim, "c_0")),
      w_draws = as_draws_matrix(subset_draws(draws_sim, "w_real")),
      sigma_draws = as_draws_matrix(subset_draws(draws_sim, "sigma")),
      x = df_real_test$w1,
      number_draws = round((beta * nrow(draws_real))),
      link = surrogate_link,
      y_log_transform = y_log_transform)
    epred_matrix_sim <- calc_pce_epred_matrix(linpred_matrix_sim, as_draws_matrix(subset_draws(draws_sim, "sigma")),
                                              surrogate_likelihood)
    linpred_matrix_weighted <- rbind(linpred_matrix_real, linpred_matrix_sim)
    epred_matrix_weighted <- rbind(epred_matrix_real, epred_matrix_sim)
    if (beta < 1){
      sigma_draws_real <- as_draws_matrix(subset_draws(
        draws_real, "sigma", iteration = 1:((1 - beta)*nrow(draws_real))))
    } else{
      sigma_draws_real <- matrix(nrow = 0, ncol= 1)
    }
    if (beta > 0){
      sigma_draws_sim <- as_draws_matrix(subset_draws(
        draws_sim, "sigma", iteration = 1:(beta*nrow(draws_real))))
    } else{
      sigma_draws_sim <- matrix(nrow = 0, ncol= 1)
    }
    sigma_draws_weighted <- rbind(sigma_draws_real, sigma_draws_sim)
    test_metrics <- get_metrics(linpred_matrix_weighted, epred_matrix_weighted,
                                sigma_draws_weighted, y_true = df_real_test$y,
                                y_true_noisy = df_real_test$y_noisy,
                                likelihood = surrogate_likelihood)
    test_elpd <- test_metrics$elpd$estimates
    test_rmse <- test_metrics$rmse
    test_list[[length(test_list)+1]] <- data.frame(
      "real_test_elpd" = test_elpd[[1]]/nrow(df_real_test),
      "real_test_elpd_se" = test_elpd[[3]]/nrow(df_real_test),
      "real_test_rmse" = test_rmse,
      "beta" = beta,
      "test_scenario" = test_scenario,
      "data_integration_scheme" = data_integration_scheme)
  }
}
# save results
pp_df <- bind_rows(pp_list)
saveRDS(pp_df, file.path(results_path, paste0("pp_df_", seed, ".rds")))
# saveRDS(pp_df, file.path(results_path, paste0("pp_df_", seed, "_ndraws1000.rds")))
surrogate_coeffs_df <- bind_rows(surrogate_coeffs_list)
saveRDS(surrogate_coeffs_df, file.path(results_path, paste0("surrogate_coeffs_df_", seed, ".rds")))
inference_params_df <- bind_rows(inference_params_list)
saveRDS(inference_params_df, file.path(results_path, paste0("inference_params_df_", seed, ".rds")))
test_df <- bind_rows(test_list)
saveRDS(test_df, file.path(results_path, paste0("test_df_", seed, ".rds")))

# read results
pp_df <- readRDS(file.path(results_path, paste0("pp_df_", seed, ".rds")))
surrogate_coeffs_df <- readRDS(file.path(results_path, paste0("surrogate_coeffs_df_", seed, ".rds")))
inference_params_df <- readRDS(file.path(results_path, paste0("inference_params_df_", seed, ".rds")))
test_df <- readRDS(file.path(results_path, paste0("test_df_", seed, ".rds")))

levels_data_integration <- c(
  "Power Scaling" = "power_scaling",
  "Post Pred Weighting" = "pp_weighting"
)


# create posterior epred plot over beta and data_integration_scheme
source(file.path(here(),"R/plot.R"))
pred_var <- "mu_pred"
plot_pp_pce_facet(pp_df, pred_var,
  betas=betas_plot,
  ylims=ylims_plot,
  w_breaks=w_breaks)+
  xlab(xlab)+
  ylab(ylab_epred)

ggsave(file.path(results_path, paste0(
  "posterior_pred_spaghetti_TRUE_predvar_",
  pred_var, "_", case_study, ".png")),
  height=8, width=15)
ggsave(file.path(results_path, paste0(
  "posterior_pred_spaghetti_TRUE_predvar_",
  pred_var, "_", case_study, ".pdf")),
  height=8, width=15)

# to receive the CIs for more than 200 posterior draws, the above code has to be rerun
# with pp_ndraws <- 1000
pp_1000_df <- readRDS(file.path(results_path, paste0("pp_df_", seed, "_ndraws1000.rds")))
p <- c(0.005, 0.05, 0.25, 0.5, 0.75, 0.95, 0.995)
p_names <- map_chr(p, ~paste0(.x*100, "%"))

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

pp_df_ci_quantiles <- pp_1000_df %>%
  filter(beta %in% betas_plot)%>%
  mutate(
    data_integration_scheme = fct_recode(data_integration_scheme,
                                          !!!levels_data_integration
    ))%>%
  group_by(x_plot, beta, data_integration_scheme) %>%
  summarize(across(c(pred, mu_pred), p_funs))

pp_df_ci_quantiles <- pp_df_ci_quantiles %>%
  mutate(
    x_plot = scale_from_1(x_plot, lb = 1, ub = config$data$x_lims[2])
  )

pp_df_ci_quantiles %>%
  ggplot(aes(x = x_plot)) +
  geom_ribbon(aes(ymin = pmax(`pred_0.5%`, ylims_plot[1]), ymax = pmin(`pred_99.5%`, ylims_plot[2]), fill = "99% CI"), alpha = 0.15) +
  geom_ribbon(aes(ymin = pmax(`pred_5%`, ylims_plot[1]), ymax = pmin(`pred_95%`, ylims_plot[2]), fill = "90% CI"), alpha = 0.25) +
  geom_ribbon(aes(ymin = pmax(`pred_25%`, ylims_plot[1]), ymax = pmin(`pred_75%`, ylims_plot[2]), fill = "50% CI"), alpha = 0.35) +
  geom_line(aes(y = `pred_50%`, color = "Median Prediction")) +
  facet_grid(cols = vars(beta), rows = vars(data_integration_scheme),
             labeller = label_bquote(cols = beta == .(beta)))+
  geom_point(aes(x = x1, y = y_noisy), df_real, color = "blue", alpha=1, size=1)+
  geom_point(aes(x = x1, y = y_noisy), df_real_ood, color = "lightblue", alpha=1, size=1)+
  scale_fill_manual(values = c("99% CI" = "red", "90% CI" = "red", "50% CI" = "red"), name=NULL) +
  scale_color_manual(values = c("Median Prediction" = "red"), name=NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
  # scale_x_continuous(breaks = c(1, 7, 14))+
  ylim(ylims_plot)+
  labs(y = ylab_pred, x= xlab)+
  theme(legend.position = "bottom")
ggsave(file.path(results_path, paste0(
  "posterior_pred_spaghetti_FALSE_", case_study, ".pdf")),
  height=8, width=15)

levels_test_scenarios <- c("OOS", "OOD", 
            "OOS/OOD (1/1)", "OOS/OOD (5/1)"
            )
test_df$test_scenario[test_df$test_scenario == "OOS+OOD"] <- "OOS/OOD (1/1)"
test_df$test_scenario[test_df$test_scenario == "OOS+OOD2"] <- "OOS/OOD (5/1)"
test_df_plot <- test_df %>%
  mutate(
    data_integration_scheme = fct_recode(data_integration_scheme,
                                         !!!levels_data_integration
    ),
    test_scenario = factor(test_scenario, levels=levels_test_scenarios))
if (model_real == "sir_2d_misspec") {
  test_df_plot <- test_df_plot %>%
    filter(test_scenario %in% c("OOS", "OOD"))
}

# Predictive performance on test sets via ELPD
elpd_plot <- ggplot(data = test_df_plot, aes(x=beta, y=real_test_elpd, color=data_integration_scheme))+
  facet_grid(vars(test_scenario), scales = "free_y")+
  geom_point()+
  geom_line()+
  # geom_errorbar(aes(ymin = real_test_elpd-real_test_elpd_se, ymax = real_test_elpd+real_test_elpd_se))+
  # geom_vline(data=filter(test_df_plot, test_scenario=="OOD"), aes(xintercept=optimal_beta_stacking$OOD)) +
  # geom_vline(data=filter(test_df_plot, test_scenario=="OOS"), aes(xintercept=optimal_beta_stacking$OOS)) +
  # geom_vline(data=filter(test_df_plot, test_scenario=="OOS+OOD"), aes(xintercept=optimal_beta_stacking$`OOS+OOD`)) +
  # ggtitle(paste0("N_ood = ", nrow(df_real_ood), " N_oos = ", nrow(df_real_oos),
  #                " N_ood_oos = ", nrow(df_real_oos_ood)))+
  labs(x = TeX("$beta$"), y = TeX("ELPD"))+
  theme(legend.title = element_blank())

# Predictive performance on test sets via RMSE
rmse_plot <- ggplot(data = test_df_plot, aes(x=beta, y=real_test_rmse, color=data_integration_scheme))+
  facet_grid(vars(test_scenario), scales = "free_y")+
  geom_point()+
  geom_line()+
  scale_y_continuous(trans='log10')+
  labs(x = TeX("$beta$"), y = TeX("RMSE"))+
  theme(legend.title = element_blank())

(elpd_plot + rmse_plot) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(file.path(results_path, paste0("elpd_rmse_", case_study, ".png")), height=9, width=15)
ggsave(file.path(results_path, paste0("elpd_rmse_", case_study, ".pdf")), height=9, width=15)

test_ps_df <- test_df %>%
  filter(
    test_scenario %in% c("OOD"),
    beta > 0,
    data_integration_scheme == "power_scaling"
  )
rmse_plot_zoom <- ggplot(data = test_ps_df, aes(x=beta, y=real_test_rmse, color=data_integration_scheme))+
  geom_point()+
  geom_line()+
  labs(x = TeX("$beta$"), y = TeX("$RMSE$ on real data"))+
  theme(legend.title = element_blank(), legend.position = 'bottom')

data <- get_surrogate_stan_list(df_sim, df_real, x_idxs, w_idxs, 1,
                               poly_degree, M, poly_idx,
                               w_prior_mean_sim, w_prior_sigma_sim,
                               w_lims, surrogate_link = surrogate_link,
                               y_log_transform = y_log_transform,
                               y_scale = y_scale)
# posterior distribution of omega_R
omega_posterior_plot <- ggplot(
  inference_params_df %>%
    filter(beta %in% c(0, 0.1, 0.4, 0.5, 0.75, 1))
    )+
  geom_density(aes(`w_real[1]`, color = beta))+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(-1, 1, length.out=100),
                y = dtruncnorm(seq(-1, 1, length.out=100), -1, 1, data$w_prior_mean, data$w_prior_sigma)))+
  geom_vline(xintercept = df_real$w2)+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'density', x=TeX('$omega_R$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))

# posterior distribution of sigma_R
sigma_posterior_plot <- ggplot(
  inference_params_df %>%
    filter(beta %in% c(0, 0.1, 0.4, 0.5, 0.75, 1)))+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(0, max(inference_params_df$sigma), length.out=100),
                y = dtruncnorm(seq(0, max(inference_params_df$sigma), length.out=100), 0, Inf, 0, 0.5)))+
  geom_density(aes(`sigma`, color = beta))+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'density', x=TeX('$sigma_R$'))+
  facet_grid(vars(beta),
            labeller = label_bquote(rows = beta == .(beta)))

(omega_posterior_plot + sigma_posterior_plot) +
  plot_layout(guides = "collect") & theme(legend.position = 'none')
ggsave(file.path(results_path, "omega_R_sigma_R_posterior_over_beta.png"), height=8, width=15)
ggsave(file.path(results_path, "omega_R_sigma_R_posterior_over_beta.pdf"), height=8, width=15)


## additional plots ##

omega_tilde_posterior_plot <- ggplot(
  surrogate_coeffs_df
)+
  geom_density(aes(`w_real[1]`, color = beta))+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(-1, 1, length.out=100),
                y=dnorm(seq(-1, 1, length.out=100), data$w_prior_mean, data$w_prior_sigma)))+
  geom_vline(xintercept = df_real$w2)+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$\\tilde{omega}_R$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))

sigma_tilde_posterior_plot <- ggplot(
  surrogate_coeffs_df)+
  geom_density(aes(sigma, color = beta))+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$\\tilde{sigma}_R$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))

((omega_tilde_posterior_plot+omega_posterior_plot) / (sigma_tilde_posterior_plot+sigma_posterior_plot)) +
  plot_layout(guides = "collect") & theme(legend.position = 'none')
ggsave(file.path(results_path, "tilde_omega_R_sigma_R_posterior_over_beta.png"), height=16, width=15)
ggsave(file.path(results_path, "tilde_omega_R_sigma_R_posterior_over_beta.pdf"), height=16, width=15)


# posterior distribution of sigma_R
sigma_tilde_posterior_plot <- ggplot(
    surrogate_coeffs_df %>%
    # filter(beta == 1)
    filter(!beta %in% c(0))
)+
  geom_density(aes(`sigma`, color = beta))+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$sigma_R$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))

pp_df[seq(1, nrow(pp_df), 100), ] %>%
  ggplot()+
  geom_density(aes(x=w_draw, color=beta))+
  scale_colour_viridis_c(TeX("$beta$"))+
  facet_grid(vars(beta))

ggplot(surrogate_coeffs_df)+
  geom_density(aes(`c[2]`, color = factor(beta)))+
  facet_grid(vars(beta))+
  scale_colour_viridis_d(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('Surrogate paramter $c_1$'))
#scale_colour_gradient(low = "yellow", high = "red", na.value = NA)
ggsave(file.path(results_path, "surrogate_paramter_c1_posterior_over_beta.png"), height=5, width=8)

# uncomment to draw a circle in the scatter plot with a radius 2 * sigma
# of the prior over the surrogate paramters c
# prior_c <- data.frame(
#   x0 = 0,
#   y0 = 0,
#   r = 2*5
# )

ggplot(surrogate_coeffs_df)+
  geom_point(aes(x=`c[1]`, y=`c[2]`, color = factor(beta)))+
  facet_grid(vars(beta))+
  scale_colour_viridis_d(TeX("$beta$"))+
  labs(x = TeX('Surrogate paramter $c_1$'), y=TeX('Surrogate paramter $c_2$'))+
  coord_fixed()
ggsave(file.path(results_path, "scatter_c_1_c_2_posterior_over_beta.png"), height=10, width=5)

scatter_tilde_omega_sigma <- ggplot(surrogate_coeffs_df)+
  geom_point(aes(x=`w_real[1]`, y=sigma), alpha=0.2)+
  facet_grid2(vars(beta), scales = "free", independent = "all",
              labeller = label_bquote(rows = beta == .(beta)))+
  labs(x = TeX('$\\tilde{omega}_R$'), y=TeX('$\\tilde{sigma}$'))+
  ggtitle("First Step")

scatter_omega_sigma <- ggplot(inference_params_df)+
  geom_point(aes(x=`w_real[1]`, y=sigma, color = rhat>1.05), alpha=0.2)+
  scale_color_manual(labels = c("Rhat<=1.05", "Rhat>1.05"), values = c("blue", "red")) +
  facet_grid2(vars(beta), scales = "free", independent = "all",
             labeller = label_bquote(rows = beta == .(beta)))+
  labs(x = TeX('$omega_R$'), y=TeX('$sigma_R$'))+
  ggtitle("Second Step")

(scatter_tilde_omega_sigma + scatter_omega_sigma)+
  plot_layout(guides = "collect")

ggsave(file.path(results_path, "scatter_omega_sigma_draws_over_beta.png"), height=15, width=10)
ggsave(file.path(results_path, "scatter_omega_sigma_draws_over_beta.pdf"), height=15, width=10)

# rhat over beta
inference_params_df %>%
  group_by(beta)%>%
  mutate(
    percenatage_bad_rhats = sum(rhat > 1.05)/1000
  ) %>%
  ggplot()+
  geom_point(aes(x=beta, y=percenatage_bad_rhats))

# Posterior plots of sigma over beta using only one-step
# (without second inference step)
ggplot(surrogate_coeffs_df)+
  geom_density(aes(`sigma`, color = beta))+
  facet_grid(vars(data_integration_scheme))+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$sigma$'))+
  facet_grid(vars(beta))+
  xlim(0, 1.2)+
  ylim(0, 18)+
  ggtitle("One-Step (Power Scaling)")
#xlim(0, 0.2)
ggsave(file.path(results_path, "sigma_posterior_over_beta_facett.png"), height=10, width=10)

# Posterior plots of omega_R over beta using only one-step
# (without second inference step)
surrogate_coeffs_df%>%
  ggplot()+
  geom_density(aes(`w_real[1]`, color = beta, group = beta))+
  facet_grid(vars(beta))+
  scale_colour_viridis_c(TeX("$beta$"), limits=c(0, 1))+
  labs(y = 'posterior distribution', x=TeX('$omega_R$'))+
  geom_vline(xintercept = df_real$w2)+
  xlim(-1, 1)
