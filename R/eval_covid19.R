library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(stringr)
library(loo)
library(forcats)
library(truncnorm)

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
df_real_ood$y <- df_real_ood$y_noisy
real_test_list <- list("OOD" = df_real_ood)
# run eval.R

# plot elpd
elpd_plot +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(file.path(results_path, paste0("elpd_", case_study, ".pdf")), height=5, width=7)


# plot inference posteriors
# posterior distribution of omega_R
omega_1_posterior_plot <- ggplot(
  inference_params_df%>%
  filter(beta %in% betas_plot)
)+
  geom_density(aes(`w_real[1]`, color = beta))+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(-1, 1, length.out=100),
                y = dtruncnorm(seq(-1, 1, length.out=100), -1, 1, data$w_prior_mean[1], data$w_prior_sigma[1])))+
  geom_vline(xintercept = df_real$w2)+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$omega_{R1}$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))

omega_2_posterior_plot <- ggplot(
  inference_params_df%>%
    filter(beta %in% betas_plot)
)+
  geom_density(aes(`w_real[2]`, color = beta))+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(-1, 1, length.out=100),
                y = dtruncnorm(seq(-1, 1, length.out=100), -1, 1, data$w_prior_mean[1], data$w_prior_sigma[1])))+
  geom_vline(xintercept = df_real$w2)+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$omega_{R2}$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))


# posterior distribution of sigma_R
sigma_posterior_plot <- ggplot(
  inference_params_df%>%
    filter(beta %in% betas_plot))+
  # geom_density(aes(x = value), linetype = "dashed",
  #           data =
  #             data.frame(
  #               value=rexp(4000, rate = 10)))+
  geom_density(aes(`sigma`, color = beta))+
  scale_colour_viridis_c(TeX("$beta$"))+
  labs(y = 'posterior distribution', x=TeX('$sigma_R$'))+
  facet_grid(vars(beta),
             labeller = label_bquote(rows = beta == .(beta)))
#  ylim(0, 13)

(omega_1_posterior_plot + omega_2_posterior_plot + sigma_posterior_plot) +
  plot_layout(guides = "collect") & theme(legend.position = 'none')
ggsave(file.path(results_path, "omega_R_sigma_R_posterior_over_beta.png"), height=8, width=15)
ggsave(file.path(results_path, "omega_R_sigma_R_posterior_over_beta.pdf"), height=8, width=15)

# plot omega priors
omega_1_prior_plot <- ggplot()+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(-1, 1, length.out=100),
                y = dtruncnorm(seq(-1, 1, length.out=100), -1, 1, data$w_prior_mean[1], data$w_prior_sigma[1])))+
  labs(y = 'prior distribution', x=TeX('$omega_{R1}$'))
omega_2_prior_plot <- ggplot()+
  geom_line(aes(x = x, y = y), linetype = "dashed",
            data =
              data.frame(
                x=seq(-1, 1, length.out=100),
                y = dtruncnorm(seq(-1, 1, length.out=100), -1, 1, data$w_prior_mean[2], data$w_prior_sigma[2])))+
  labs(y = 'prior distribution', x=TeX('$omega_{R2}$'))
omega_1_prior_plot + omega_2_prior_plot
