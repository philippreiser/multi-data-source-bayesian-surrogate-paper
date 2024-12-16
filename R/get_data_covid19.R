library("COVID19")

get_covid19_df <- function(x_lims, x_start, x_end){
  x1_lb <- x_lims[1]
  x1_ub <- x_lims[2]
  t_start <- x_start
  t_end <- x_end
  date_start <- as.Date("2020-02-24") + x_start - 1
  date_end <- as.Date("2020-02-24") + x_end - 1
  covid_data <- covid19("Italy", level = 1, start = date_start, end = date_end, verbose = FALSE, dir="data_covid19")
  infected <- covid_data$confirmed - covid_data$recovered
  N_t <- t_end - t_start + 1
  covid19_df <- data.frame(
    x1 = seq(from = x_start, to = x_end))
  covid19_df <- covid19_df %>%
    mutate(
      N = N_t,
      y_noisy = infected,
      w1 = scale_to_1(x1, lb = x1_lb, ub = x1_ub), # rename to x1_scaled
      population = covid_data$population[1]
    ) %>%
    mutate(
      log_y = log(y_noisy),
      y_scaled = scale_to_1(y_noisy, lb=0, ub=600)
    )
  covid19_df$y_noisy <- covid19_df$y_noisy/covid19_df$population * 100000
  covid19_df$population <- 100000
  return(covid19_df)
}

get_sir_model_df <- function(N, model, sigma=0, sample_x="sobol",
                            prior_mean=0, prior_sigma=1,
                            x_lims=c(-pi, pi), x_train_percentages=c(0, 1),
                            w_lims=c(-1, 1), scale_prior=TRUE,
                            population=1000, i0=1
){
  x1_lb <- x_lims[1]
  x1_ub <- x_lims[2]
  x2_lb <- w_lims$w_1[1]
  x2_ub <- w_lims$w_1[2]
  x3_lb <- w_lims$w_2[1]
  x3_ub <- w_lims$w_2[2]
  x1_lower_percentage <- x_train_percentages[1]
  x1_upper_percentage <- x_train_percentages[2]
  if (scale_prior){
    # scale the mean and prior of a normal distribution from [lb, ub] to [0, 1]
    w_1_prior_mean <- scale_from_1(scale_to_1(prior_mean$w_1, x2_lb, x2_ub), 0, 1)
    w_1_prior_sigma <- prior_sigma$w_1 / (x2_ub - x2_lb)
    w_2_prior_mean <- scale_from_1(scale_to_1(prior_mean$w_2, x3_lb, x3_ub), 0, 1)
    w_2_prior_sigma <- prior_sigma$w_2 / (x3_ub - x3_lb)
  }
  get_model <- get(paste0("get_", model))
  # first the input is sampled on the scale [0, 1]
  if (sample_x == "uniform"){
    tmp <- matrix(c(runif(N), runif(N), runif(N)), ncol=3)
  } else if (sample_x == "sobol"){
    tmp <- sobolSequence.points(3, count = N)
  } else if (sample_x == "x1_slice"){
    tmp <- matrix(
      c(seq(from=x1_lower_percentage, to=x1_upper_percentage,length.out=N),
        rep(w_1_prior_mean, N),
        rep(w_2_prior_mean, N)),
      ncol=3)
  } else if (sample_x == "x1_uniform"){
    tmp <- matrix(
      c(runif(N, min = x1_lower_percentage, max=x1_upper_percentage),
        rep(w_1_prior_mean, N),
        rep(w_2_prior_mean, N)),
      ncol=3)
  }
  tmp <- tmp %>%
    as.data.frame()%>%
    rename(x1 = V1, x2 = V2, x3 = V3) %>%
    mutate(
      x1 = scale_from_1(scale_to_1(x1, lb = 0, ub = 1), lb = x1_lb, ub = x1_ub),
      x2 = scale_from_1(scale_to_1(x2, lb = 0, ub = 1), lb = x2_lb, ub = x2_ub),
      x3 = scale_from_1(scale_to_1(x3, lb = 0, ub = 1), lb = x3_lb, ub = x3_ub)
    )%>%
    mutate(
      N = N,
      y = get_sir(x1, x2, x3, N=population, i0=i0) + rnorm(N, sd=sigma),
      w1 = scale_to_1(x1, lb = x1_lb, ub = x1_ub), # rename to x1_scaled
      w2 = scale_to_1(x2, lb = x2_lb, ub = x2_ub), # rename to x2_scaled
      w3 = scale_to_1(x3, lb = x3_lb, ub = x3_ub) # rename to x3_scaled
    ) %>%
    mutate(
      log_y = log(y),
      y_scaled = scale_to_1(y, lb=0, ub=600) + rnorm(N, sd=sigma)
    )
  return(tmp)
}
