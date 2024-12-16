library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(future.apply)

source(file.path(here(),"R/get_data.R"))
source(file.path(here(),"R/plot.R"))
source(file.path(here(), "R/utils.R"))
source(file.path(here(), "R/inference.R"))
source(file.path(here(), "R/load_config_args.R"))

results_path <- file.path("results", experiment_name)
if (!dir.exists(results_path)){
  dir.create(results_path, recursive=TRUE)
}
file.copy(file.path(config_dir, config_file_name), results_path)

df_sim <- get_2d_model_df(N_sim,
                          model = model_sim,
                          noise_model = noise_model_sim,
                          sigma = sigma_sim,
                          sample_x = sample_x_sim,
                          prior_mean = w_prior_mean_sim,
                          prior_sigma = w_prior_sigma_sim,
                          x_lims = x_lims, w_lims = w_lims,
                          add_one = add_one)
saveRDS(df_sim, file=get_data_save_name(results_path, model_sim, "simulation"))

df_real <- get_2d_model_df(N_real,
                           model = model_real,
                           noise_model = noise_model_real,
                           sigma = sigma_real,
                           sample_x = sample_x_real,
                           prior_mean = w_prior_mean_real,
                           prior_sigma = w_prior_sigma_real,
                           x_lims = x_lims, w_lims = w_lims,
                           x_train_percentages = x_train_percentages,
                           add_one = add_one)

saveRDS(df_real, file=get_data_save_name(results_path, model_real, "real", sample_x_real))
df_real_ood <- get_2d_model_df(100,
                               model = model_real,
                               noise_model = noise_model_real,
                               sigma = sigma_real,
                               sample_x = sample_x_real_ood,
                               prior_mean = w_prior_mean_real,
                               prior_sigma = w_prior_sigma_real,
                               x_lims = x_lims, w_lims = w_lims,
                               x_train_percentages = c(x_train_percentages[2], 1),
                               add_one = add_one)
saveRDS(df_real_ood, file=get_data_save_name(results_path, model_real, "real_ood", sample_x_real_ood))
df_real_oos <- get_2d_model_df(100,
                               model = model_real,
                               noise_model = noise_model_real,
                               sigma = sigma_real,
                               sample_x = sample_x_real_ood,
                               prior_mean = w_prior_mean_real,
                               prior_sigma = w_prior_sigma_real,
                               x_lims = x_lims, w_lims = w_lims,
                               x_train_percentages = x_train_percentages,
                               add_one = add_one
)
saveRDS(df_real_oos, file=get_data_save_name(results_path, model_real, "real_oos", sample_x_real_ood))

df_real_oos_ood <- rbind(df_real_oos, df_real_ood)
saveRDS(df_real_oos_ood, file=get_data_save_name(results_path, model_real, "real_oos_ood", sample_x_real_ood))


# sample directly from df_real_oos and df_real_ood
df_real_ood_2 <- sample_n(df_real_ood, 17)
df_real_oos_2 <- sample_n(df_real_oos, 83)
# Combine the samples
df_real_oos_ood_2 <- bind_rows(df_real_ood_2, df_real_oos_2)

saveRDS(df_real_oos_ood_2, file=get_data_save_name(results_path, model_real, "real_oos_ood_2", sample_x_real_ood))

plot_sim_real_model(config, quantiles=c(0.05, 0.5, 0.95))
  
ggsave(file.path(results_path, "sim_real_model_real_train_test.pdf"), height=6, width=12)

for (data_integration_scheme in data_integration_schemes){
  for (i in seq_along(betas)){
    beta <- betas[[i]]
    w_idxs <- setdiff(c(1:M), x_idxs)
    data <- get_surrogate_stan_list(df_sim, df_real, x_idxs, w_idxs, beta,
                                   poly_degree, M, poly_idx,
                                   w_prior_mean_sim, w_prior_sigma_sim,
                                   w_lims,
                                   surrogate_likelihood = surrogate_likelihood,
                                   surrogate_link = surrogate_link,
                                   y_log_transform = y_log_transform,
                                   y_scale = y_scale)
    if (surrogate_likelihood == "lognormal"){
      surrogate_file <- file.path(
        here(), paste0("stan_code/pce_1-step_", data_integration_scheme, "_", "normal", ".stan"))
    } else{
      surrogate_file <- file.path(
        here(), paste0("stan_code/pce_1-step_", data_integration_scheme, "_", surrogate_likelihood, ".stan"))
    }
    surrogate_model <- cmdstan_model(surrogate_file)
    surrogate_model_fit <- surrogate_model$sample(
      data = data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      adapt_delta = adapt_delta,
      init = 0
    )
    # save fit
    surrogate_model_fit$save_object(
      file = get_fit_file_name(results_path, data_integration_scheme, N_real,
                              beta))

    # mcmc diagnostics
    surrogate_fit_np <- nuts_params(surrogate_model_fit)
    mcmc_trace(surrogate_model_fit$draws(), regex_pars=c("c", "sigma", "w_real"),
                     np = surrogate_fit_np)
    ggsave(file.path(results_path, sprintf("trace-surrogate_%s-N_real_%s-beta_%.2f.png", data_integration_scheme, N_real, beta)), width=10, height=8)
    surrogate_fit_rhats <- rhat(surrogate_model_fit)
    mcmc_rhat(surrogate_fit_rhats[grepl( "c|sigma|w_real" , names( surrogate_fit_rhats ) ) ]) + yaxis_text(hjust = 1)
    ggsave(file.path(results_path, sprintf("rhat-surrogate_%s-N_real_%s-beta_%.2f.png", data_integration_scheme, N_real, beta)), width=10, height=8)

    ### inference step for sigma/omega using real data ###
    # inference for sigma and omega using real data
    if (surrogate_likelihood == "lognormal"){
      surrogate_inference_file <- file.path(here(), paste0("stan_code/pce_1-step_power_scaling_omega_sigma_inference_",
                                                          "normal", ".stan"))
    } else{
      surrogate_inference_file <- file.path(here(), paste0("stan_code/pce_1-step_power_scaling_omega_sigma_inference_",
                                                          surrogate_likelihood, ".stan"))
    }
    surrogate_inference_model <- cmdstan_model(surrogate_inference_file)
    surrogate_inference_model_fit <- inference_step(
      surrogate_model_fit,
      surrogate_inference_model,
      NULL,
      data,
      number_samples_1 = iter_sampling*chains,
      number_samples_2 = inner_iter_sampling,
      number_chains_2 = inner_chains,
      type = inner_mcmc_type
    )

    surrogate_inference_model_fit$save_object(
      file = get_fit_file_name(results_path, paste0(data_integration_scheme, "_inference"), N_real,
                               beta))
    print(mcmc_dens(surrogate_model_fit$draws(), "sigma") /
            mcmc_dens(surrogate_inference_model_fit$draws(), "sigma"))
    print(mcmc_dens(surrogate_model_fit$draws(), "w_real[1]") /
            mcmc_dens(surrogate_inference_model_fit$draws(), "w_real[1]"))
  }
}
# posterior predictive weighting for beta = 0 without omega input
data <- get_surrogate_stan_list(df_sim, df_real, x_idxs, w_idxs, 0,
  poly_degree, M, poly_idx,
  w_prior_mean_sim, w_prior_sigma_sim,
  w_lims,
  surrogate_likelihood = surrogate_likelihood,
  surrogate_link = surrogate_link,
  y_log_transform = y_log_transform,
  y_scale = y_scale,
  only_x = TRUE)
surrogate_model_fit <- surrogate_model$sample(
  data = data,
  seed = seed,
  chains = chains,
  parallel_chains = parallel_chains,
  iter_sampling = iter_sampling,
  iter_warmup = iter_warmup,
  adapt_delta = adapt_delta
)
# save fit
surrogate_model_fit$save_object(
  file = get_fit_file_name(results_path, "pp_weighting", N_real,
  0))
