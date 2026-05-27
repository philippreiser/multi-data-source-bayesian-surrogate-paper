args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  # Case Study 1.1
  # config_file_name <- "log_sin.yml"
  
  # Case Study 1.2
  # config_file_name <- "log_sin_train_0.7.yml"

  # Case Study 2.1
  # config_file_name <- "sir_2d.yml"

  # Case Study 2.2
  # config_file_name <- "sir_covid19.yml"

  # debug
  config_file_name <- "log_sin_debug.yml"
} else{
  config_file_name <- args[1]
}
config_dir <- file.path(here(), "config")
config <- yaml::read_yaml(file.path(config_dir, config_file_name), eval.expr=TRUE)

seed <- config$seed
set.seed(seed)

# data parameters
N_sim <- config$data$N_sim
N_real <- config$data$N_real
sigma_sim <- config$data$sigma_sim
sigma_real <- config$data$sigma_real
model_sim <- config$data$model_sim
model_real <- config$data$model_real
noise_model_sim <- config$data$noise_model_sim
noise_model_real <- config$data$noise_model_real
x_idxs <- config$data$x_idxs
x_lims <- config$data$x_lims
x_train_percentages <- config$data$x_train_percentages
w_lims <- config$data$w_lims
sample_x_sim <- config$data$sample_x_sim
sample_x_real <- config$data$sample_x_real # sprintf("x%i_slice", x_idxs)
sample_x_real_ood <- config$data$sample_x_real_ood # sprintf("x%i_uniform", x_idxs)
y_log_transform <- config$data$y_log_transform
y_scale <- config$data$y_scale
w_prior_mean_sim <- config$data$w_prior_mean_sim
w_prior_sigma_sim <- config$data$w_prior_sigma_sim
w_prior_mean_real <- config$data$w_prior_mean_real
w_prior_sigma_real <- config$data$w_prior_sigma_real
offset_real <- config$data$offset_real
add_one <- config$data$add_one

# surrogate parameters
data_integration_schemes <- config$surrogate_model$data_integration_schemes
betas <- config$surrogate_model$betas
M <- config$surrogate_model$M
poly_degree <- config$surrogate_model$poly_degree
poly_idx <- config$surrogate_model$poly_idx
surrogate_likelihood <- config$surrogate_model$likelihood
surrogate_link <- config$surrogate_model$link

# mcmc paramters
chains <- config$mcmc$chains
parallel_chains <- config$mcmc$parallel_chains
iter_sampling <- config$mcmc$iter_sampling
iter_warmup <- config$mcmc$iter_warmup
adapt_delta <- config$mcmc$adapt_delta

# inner mcmc parameters
inner_chains <- config$inner_mcmc$chains
inner_parallel_chains <- config$inner_mcmc$parallel_chains
inner_iter_sampling <- config$inner_mcmc$iter_sampling
inner_iter_warmup <- config$inner_mcmc$iter_warmup
inner_adapt_delta <- config$inner_mcmc$adapt_delta
inner_mcmc_type <- config$inner_mcmc$type

experiment_name <- get_experiment_name(config)

# plotting paramters
pp_ndraws <- config$plots$pp_ndraws
betas_plot <- config$plots$betas_plot
ylims_plot <- config$plots$ylims_plot
xlab <- config$plots$xlab
ylab_epred <- config$plots$ylab_epred
ylab_pred <- config$plots$ylab_pred
case_study <- config$plots$case_study
w_breaks <- config$plots$w_breaks
