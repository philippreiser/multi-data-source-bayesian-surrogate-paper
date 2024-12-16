library(ggplot2)
library(tidyr)
library(patchwork)
library(latex2exp)
library(ggh4x)
source(file.path(here(), "R/utils.R"))
theme_set(theme_bw())
theme_update(text = element_text(size = 23))

plot_pp_pce <- function(posterior_pred,
                        w_var = "w_prior_sample",
                        w_lims=NULL, y_pred = "y_pred", alpha=0.1,
                        w_breaks = c(0.8, 1)){
  pfit_1 <- ggplot() +
    geom_line(aes(x = x_plot, y = !!sym(y_pred), color=!!sym(w_var), group = w_id),
              data = posterior_pred, alpha=alpha)+
    scale_color_viridis_c(TeX("$\\omega_R$"),
                          limits=w_lims, option="F",
                          breaks=w_breaks)
  pfit_1
}

plot_pp_pce_facet <- function(pp_df, pred_var, alpha=0.1,
                              betas_plot=c(0, 0.1, 0.4, 0.75, 1.0),
                              ylims=c(-0.5, 9),
                              scale_inputs=FALSE,
                              w_breaks=c(0.8, 1)) {
  pp_df <- pp_df%>%
    filter(beta %in% betas_plot)%>%
    mutate(
      data_integration_scheme = fct_recode(data_integration_scheme,
                                            !!!levels_data_integration
      ))
  if (!scale_inputs) {
    pp_df <- pp_df %>%
      mutate(
        x_plot = scale_from_1(x_plot, lb = config$data$x_lims[1], ub = config$data$x_lims[2]),
        w_draw = scale_from_1(w_draw, lb = config$data$w_lims[1], ub = config$data$w_lims[2])
      )
    pp_plot <- plot_pp_pce(pp_df, w_var="w_draw",
                           y_pred = pred_var, alpha = alpha,
                           w_breaks = w_breaks) +
      facet_grid(cols = vars(beta), rows = vars(data_integration_scheme),
                 labeller = label_bquote(cols = beta == .(beta))) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 4))+
      # scale_x_continuous(breaks = c(1, 7, 14))+
      labs(y = 'y', x= TeX("$x$"))+
      theme(legend.position = "bottom")+
      ylim(ylims)+
      geom_point(aes(x = x1, y = y_noisy), df_real, color = "blue", alpha=1, size=1)+
      geom_point(aes(x = x1, y = y_noisy), df_real_ood, color = "lightblue", alpha=1, size=1)
  } else {
    pp_plot <- plot_pp_pce(pp_df, w_var="w_draw",
                           y_pred = pred_var, alpha = alpha,
                           w_breaks = w_breaks) +
      facet_grid(cols = vars(beta), rows = vars(data_integration_scheme),
                 labeller = label_bquote(cols = beta == .(beta))) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3))+
      labs(y = 'y', x= TeX("$x$ (scaled)"))+
      theme(legend.position = "bottom")+
      ylim(ylims)+
      geom_point(aes(x = w1, y = y_noisy), df_real, color = "blue", alpha=1, size=1)+
      geom_point(aes(x = w1, y = y_noisy), df_real_ood, color = "lightblue", alpha=1, size=1)
  }
  return(pp_plot)
}

plot_sim_real_model <- function(config, quantiles=c(0.05, 0.5, 0.95)){
  df_list <- list()
  for (i in seq_along(quantiles)){
    quantile <- quantiles[i]
    omega_i_sim <- qnorm(quantile, config$data$w_prior_mean_sim,
                     config$data$w_prior_sigma_sim)
    # df_sim <- get_2d_model_df(100, model = config$data$model_sim, sigma = 0,
    #                           sample_x = "x1_slice",
    #                           prior_mean = omega_i_sim, prior_sigma = 0,
    #                           x_lims = config$data$x_lims, w_lims = config$data$w_lims)
    df_sim <- get_2d_model_df(100,
                              model = model_sim,
                              noise_model = noise_model_sim,
                              sigma = 0,
                              sample_x = "x1_slice",
                              prior_mean = omega_i_sim,
                              prior_sigma = 0,
                              x_lims = x_lims, w_lims = w_lims,
                              add_one = add_one)
    df_sim <- data.frame(
      y = df_sim$y,
      x = df_sim$x1,
      #omega = df_sim$w2,
      quantile = quantile,
      data_source = "Simulation Model"
    )
    df_list[[length(df_list)+1]] <- df_sim
    omega_i_real <- qnorm(quantile, config$data$w_prior_mean_real,
                         config$data$w_prior_sigma_real)
    df_synthetic_truth <- get_2d_model_df(100,
                               model = model_real,
                               noise_model = "none",
                               sigma = 0,
                               sample_x = "x1_slice",
                               prior_mean = omega_i_real,
                               prior_sigma = 0,
                               x_lims = x_lims, w_lims = w_lims,
                               add_one = add_one)
    df_synthetic_truth <- data.frame(
      y = df_synthetic_truth$y,
      x = df_synthetic_truth$x1,
      quantile = quantile,
      data_source = "Synthetic Truth"
    )
    df_list[[length(df_list)+1]] <- df_synthetic_truth
  }
  df_plot <- bind_rows(df_list)
  df_plot <- df_plot %>%
    pivot_wider(names_from = quantile, values_from = y)
  colors <- c(
    "Simulation Model" = "black", "Synthetic Truth" = "blue",
    "Real Train Data" = "blue", "Real Test Data (OOS)" = "violet",
    "Real Test Data (OOD)" = "lightblue")
  p_sim_real <- ggplot(data = df_plot)+
    geom_line(aes(x = x, y=`0.5`, color=data_source))+
    geom_ribbon(aes(x = x, ymin = `0.05`, ymax = `0.95`, group=data_source), fill="grey", alpha=0.2)+
    geom_point(aes(x = x1, y = y_noisy, color = "Real Train Data"), data=df_real) +
    geom_point(aes(x = x1, y = y_noisy, color = "Real Test Data (OOS)"), data=df_real_oos)+
    geom_point(aes(x = x1, y = y_noisy, color = "Real Test Data (OOD)"), data=df_real_ood)+
    theme(legend.position = "bottom",
          legend.text = element_text(size = 16))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
    guides(fill=guide_legend(nrow = 2, byrow = TRUE))+
    scale_color_manual(values = colors,
                       name = NULL,
                       breaks = c(
                         "Simulation Model", 
                         "Synthetic Truth", 
                         "Real Train Data", 
                         "Real Test Data (OOS)", 
                         "Real Test Data (OOD)"))+
    labs(y="y", x="x")
  return(p_sim_real)
}

costum_colors_1 <- c('#66CCEE', '#EE6677', '#228833', '#CCBB44', '#4477AA', '#AA3377')
