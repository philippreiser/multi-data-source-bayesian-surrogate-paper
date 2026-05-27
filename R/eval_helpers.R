eval_surrogates <- function() {
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
  return(
    list(
      pp_df,
      surrogate_coeffs_df,
      surrogate_coeffs_df,
      inference_params_df,
      test_df
    )
  )
}