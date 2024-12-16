evaluate_test_scenario <- function(df_real_test){
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

evaluate_test_scenarios <- function(real_test_list, linpred_matrix, epred_matrix) {
  for (test_scenario in names(real_test_list)){
    df_real_test <- real_test_list[[test_scenario]]
    evaluate_test_scenario
  }
}