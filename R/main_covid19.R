library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(future.apply)
library(forcats)

source(file.path(here(),"R/get_data.R"))
source(file.path(here(),"R/get_data_covid19.R"))
source(file.path(here(),"R/plot.R"))
source(file.path(here(), "R/utils.R"))
source(file.path(here(), "R/inference.R"))
source(file.path(here(), "R/load_config_args.R"))

results_path <- file.path("results", experiment_name)
if (!dir.exists(results_path)){
  dir.create(results_path, recursive=TRUE)
}
file.copy(file.path(config_dir, config_file_name), results_path)

df_real <- get_covid19_df(x_lims, 1, N_real)

saveRDS(df_real, file=get_data_save_name(results_path, model_real, "real", sample_x_real))

df_real_ood <- get_covid19_df(x_lims, N_real+1, x_lims[2])
saveRDS(df_real_ood, file=get_data_save_name(results_path, model_real, "real_ood", sample_x_real_ood))

colors <- c("real_train" = "blue", "sim" = "black", "real_test" = "lightblue")

df_sim <- get_sir_model_df(N_sim,
                           model = model_sim,
                           #noise_model = noise_model_sim,
                           sigma = sigma_sim,
                           sample_x = sample_x_sim,
                           prior_mean = w_prior_mean_sim,
                           prior_sigma = w_prior_sigma_sim,
                           x_lims = x_lims, w_lims = w_lims,
                           population = df_real$population[1], i0 = df_real$y_noisy[1])
saveRDS(df_sim, file=get_data_save_name(results_path, model_sim, "simulation"))
ggplot()+
  geom_point(aes(x = w1, y = y_noisy, color = "real_train"), data=df_real) +
  geom_point(aes(x = w1, y = y, color = "sim"), data=df_sim) +
  geom_point(aes(x = w1, y = y_noisy, color = "real_test"), data=df_real_ood)+
  labs(x = "x",
       y = "y",
       color = "Data") +
  scale_color_manual(values = colors)
ggsave(file.path(results_path, "sim_real_real_ood_data.png"), height=5, width=10)
ggplot()+
  geom_point(aes(x = x1, y = y_noisy, color = "real_train"), data=df_real) +
  geom_point(aes(x = x1, y = y_noisy, color = "real_test"), data=df_real_ood)+
  labs(x = "Day",
       y = "# Infections * 10^5",
       color = "Data") +
  scale_color_manual(values = colors)
ggsave(file.path(results_path, "real_real_ood_data.png"), height=5, width=10)

spaghetti <- TRUE
pp_ndraws <- 200
x_plot <- seq(from=-1, to=1,length.out=100)
levels_data_integration <- c(
  "Power Scaling" = "power_scaling",
  "Posterior Pred Weighting" = "pp_weighting"
)
var_dim <- "w1"
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
    if (beta > 0){
      print(mcmc_dens(surrogate_model_fit$draws(c("w_real"))))
    }
    print(mcmc_trace(surrogate_model_fit$draws(), regex_pars=c("c", "sigma", "w_real")))
    
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
    surrogate_inference_model_fit <- inference_e_post(surrogate_model_fit,
                                                     surrogate_inference_model,
                                                     NULL,
                                                     data,
                                                     number_samples_1 = iter_sampling*chains,
                                                     number_samples_2 = inner_iter_sampling,
                                                     number_chains_2 = inner_chains)

    surrogate_inference_model_fit$save_object(
      file = get_fit_file_name(results_path, paste0(data_integration_scheme, "_inference"), N_real,
                               beta))
    print(mcmc_dens(surrogate_model_fit$draws(), "sigma") /
            mcmc_dens(surrogate_inference_model_fit$draws(), "sigma"))
    print(mcmc_dens(surrogate_model_fit$draws(c("w_real"))) /
            mcmc_dens(surrogate_inference_model_fit$draws(c("w_real"))))
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
  adapt_delta = adapt_delta
)
# save fit
surrogate_model_fit$save_object(
  file = get_fit_file_name(results_path, "pp_weighting", N_real, 0))
