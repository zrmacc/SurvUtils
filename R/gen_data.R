# Purpose: Simulate Data.
# Updated: 2024-08-04

#' Generate Data
#' 
#' @param base_event_rate Baseline arrival rate for events.
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param covariates Numeric design matrix.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n Number of subjects. Overwritten by `nrow(covariates)` if covariates are provided.
#' @param simple Return only the index, time, and status? If FALSE, returns additional data.
#' @param tau Truncation time.
#' @return Data.frame.
#' @export
GenData <- function(
  base_event_rate = 1.0,
  beta_event = NULL,
  censoring_rate = 0.25,
  covariates = NULL,
  frailty_variance = 0.0,
  min_event_rate = 0.05,
  n = 100,
  simple = TRUE,
  tau = 4.0
) {
  
  # Generate subject-specific data frame.
  if (is.null(covariates)) {
    covariates <- data.matrix(rep(1, n))
  } else {
    covariates <- data.matrix(covariates)
    n <- nrow(covariates)
  }
  df <- data.frame(idx = seq_len(n), covariates)
  
  # Calculate subject-specific event rate.
  if (is.null(beta_event)) {beta_event <- rep(0, ncol(covariates))}
  df$true_event_rate <- base_event_rate * exp(covariates %*% beta_event)
  df$true_event_rate <- pmax(df$true_event_rate, min_event_rate)
  
  # Apply frailty.
  if (frailty_variance > 0) {
    theta <- 1 / frailty_variance
    df$frailty <- stats::rgamma(n = n, shape = theta, rate = theta)
  } else {
    df$frailty <- 1
  }
  df$true_event_rate <- df$frailty * df$true_event_rate 
  
  # Generate event times.
  df$event_time <- stats::rexp(n = n, rate = df$true_event_rate)
  df$censor_time <- stats::rexp(n = n, rate = censoring_rate)
  df$time <- pmin(df$event_time, df$censor_time, tau)
  df$status <- 1 * (df$time == df$event_time)
  
  # Output.
  if (simple) {
    idx <- time <- status <- NULL
    out <- df %>% dplyr::select(idx, time, status)
  } else {
    out <- df
  }
  return(out)
}


#' Generate Competing Risks Data
#' 
#' @param base_death_rate Baseline arrival for for death (the competing risk).
#' @param beta_death Numeric vector of log rate ratios for the death rate.
#' @param base_event_rate Baseline arrival rate for events.
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param covariates Numeric design matrix.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_death_rate Minimum subject-specific event rate. Must be non-negative.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n Number of subjects. Overwritten by `nrow(covariates)` if covariates are provided.
#' @param simple Return only the index, time, and status? If FALSE, returns additional data.
#' @param tau Truncation time.
#' @return Data.frame.
#' @export
GenCRData <- function(
  base_death_rate = 0.25,
  beta_death = NULL,
  base_event_rate = 1.0,
  beta_event = NULL,
  censoring_rate = 0.25,
  covariates = NULL,
  frailty_variance = 0.0,
  min_death_rate = 0.00,
  min_event_rate = 0.05,
  n = 100,
  simple = TRUE,
  tau = 4.0
) {
  
  # Generate subject-specific data frame.
  if (is.null(covariates)) {
    covariates <- data.matrix(rep(1, n))
  } else {
    covariates <- data.matrix(covariates)
    n <- nrow(covariates)
  }
  df <- data.frame(idx = seq_len(n), covariates)

  # Calculate subject-specific death rate.
  if (is.null(beta_death)) {beta_death <- rep(0, ncol(covariates))}
  df$true_death_rate <- base_death_rate * exp(covariates %*% beta_death)
  df$true_death_rate <- pmax(df$true_death_rate, min_death_rate)
  
  # Calculate subject-specific event rate.
  if (is.null(beta_event)) {beta_event <- rep(0, ncol(covariates))}
  df$true_event_rate <- base_event_rate * exp(covariates %*% beta_event)
  df$true_event_rate <- pmax(df$true_event_rate, min_event_rate)
  
  # Apply frailty.
  if (frailty_variance > 0) {
    theta <- 1 / frailty_variance
    df$frailty <- stats::rgamma(n = n, shape = theta, rate = theta)
  } else {
    df$frailty <- 1
  }
  df$true_death_rate <- df$frailty * df$true_death_rate
  df$true_event_rate <- df$frailty * df$true_event_rate 
  
  # Overall hazard.
  overall_hazard <- df$true_death_rate + df$true_event_rate
  
  # Arrivals.
  event_times <- stats::rexp(n = n, rate = overall_hazard)
  
  # Arrival type.
  pi <- cbind(df$true_event_rate, df$true_death_rate)
  pi <- t(apply(pi, 1, function(x){x / sum(x)}))
  event_status <- sapply(seq_len(n), function(i) {
    apply(stats::rmultinom(n = 1, size = 1, prob = pi[i, ]), 2, which.max)
  })
  
  df$event_time <- event_times
  df$event_type <- event_status

  # Censoring.
  if (censoring_rate > 0) {
    censor_time <- stats::rexp(n = n, rate = censoring_rate)
    censor_time <- pmin(censor_time, tau)
  }  else {
    censor_time <- tau
  }
  
  df$censor_time <- censor_time
  df$time <- pmin(event_times, censor_time)
  df$status <- (event_times <= censor_time) * event_status
  
  # Output.
  if (simple) {
    idx <- time <- status <- NULL
    out <- df %>% dplyr::select(idx, time, status)
  } else {
    out <- df
  }
  return(out)
}

