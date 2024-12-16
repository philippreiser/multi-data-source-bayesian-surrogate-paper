library(here)
library(future)

inference_step <- function(surrogate_model_fit, inference_model,
                           inference_model_fit=NULL,
                           data=NULL,
                           number_samples_1 = 100,
                           number_samples_2 = 250,
                           number_chains_2 = 1,
                           type = "e_post") {
  if (type == "e_post") {
    inference_fit <- inference_e_post(surrogate_model_fit, inference_model,
      inference_model_fit,
      data,
      number_samples_1,
      number_samples_2,
      number_chains_2)
  } else {
    inference_fit <- inference_point(surrogate_model_fit, inference_model,
      inference_model_fit,
      data,
      number_samples_2,
      number_chains_2)
  }
}

inference_e_post <- function(surrogate_model_fit, inference_model,
                             inference_model_fit=NULL,
                             data=NULL,
                             number_samples_1 = 100,
                             number_samples_2 = 250,
                             number_chains_2 = 1) {
  fits_csv_files <- c()
  csv_dir <- file.path(here(), "fitted_models/_imodel_fits_cmdstan", strsplit(tempdir(), "/")[[1]][3])
  dir.create(csv_dir, showWarnings = FALSE, recursive=TRUE)
  rhats <- vector("list", length(number_samples_1))
  fit_e_post <- function(i) {
    print(paste0("fitting model: ", i))
    data$c_0 <- surrogate_model_fit$draws("c_0")[i]
    data$c <- as.vector(surrogate_model_fit$draws("c", format="matrix")[i, ])
    w_real <- as.vector(surrogate_model_fit$draws("w_real", format="matrix")[i, ])
    init_value <- if (beta == 1) {
      0
    } else {
      function() {list(w_real = w_real)}
    }
    inference_model_fit <- inference_model$sample(
      data = data,
      seed = i,
      chains = number_chains_2,
      parallel_chains = 1,
      iter_sampling = number_samples_2,
      iter_warmup = iter_warmup,
      adapt_delta = adapt_delta,
      output_dir=csv_dir,
      init = init_value,
      refresh=0
    )
    return(inference_model_fit$output_files())
  }
  plan(multisession, workers=4)
  fits_epost <- future.apply::future_lapply(
    seq_along(1:number_samples_1),
    fit_e_post, future.seed = TRUE)
  fits_csv_files <- do.call(c, fits_epost)

  inference_fit <- as_cmdstan_fit(fits_csv_files)
  unlink(csv_dir, recursive = TRUE)
  return(inference_fit)
}

inference_point <- function(surrogate_model_fit, inference_model,
                             inference_model_fit=NULL,
                             data=NULL,
                             number_samples_2 = 250,
                             number_chains_2 = 4) {
  data$c_0 <- mean(surrogate_model_fit$draws("c_0"))
  data$c <- as.vector(colMeans(surrogate_model_fit$draws("c", format="matrix")))
  w_real <- as.vector(colMeans(surrogate_model_fit$draws("w_real", format="matrix")[i, ]))
  init_value <- if (beta == 1) {
    0
  } else {
    function() {list(w_real = w_real)}
  }
  inference_model_fit <- inference_model$sample(
    data = data,
    seed = seed,
    chains = number_chains_2,
    parallel_chains = number_chains_2,
    iter_sampling = number_samples_2,
    iter_warmup = iter_warmup,
    adapt_delta = adapt_delta,
    init = init_value,
    refresh=0
  )
  return(inference_model_fit)
}

inference_elik <- function(surrogate_model_fit, inference_model,
                            inference_model_fit=NULL,
                            data=NULL,
                            number_samples_2 = 250,
                            number_chains_2 = 4) {
  data$c_0 <- as.vector(surrogate_model_fit$draws("c_0", format="matrix"))
  data$c <- surrogate_model_fit$draws("c", format="matrix")
  data$sample_weights <- rep(1, length(data$c_0))/length(data$c_0)
  data$N_samples <- length(data$sample_weights)
  w_real <- as.array(mean(surrogate_model_fit$draws("w_real", format="matrix")[i, ]))
  init_value <- if (beta == 1) {
    0
  } else {
    function() {list(w_real = w_real)}
  }
  inference_model_fit <- inference_model$sample(
    data = data,
    seed = seed,
    chains = number_chains_2,
    parallel_chains = number_chains_2,
    iter_sampling = number_samples_2,
    iter_warmup = iter_warmup,
    adapt_delta = adapt_delta,
    init = init_value
  )
  return(inference_model_fit)
}
