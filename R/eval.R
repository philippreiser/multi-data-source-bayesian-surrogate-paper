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
source(file.path(here(),"R/eval_helpers.R"))
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

eval_df_list <- eval_surrogates()

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
# pp_1000_df <- readRDS(file.path(results_path, paste0("pp_df_", seed, "_ndraws1000.rds")))
plot_pp_pce_ci(pp_df)
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
