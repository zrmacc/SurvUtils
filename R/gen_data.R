#' Purpose: Simulate Data.
#' Updated: 2022-03-19

#' Generate Data
#' 
#' @param base_event_rate Baseline arrival rate for events.
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param covariates Numeric design matrix.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n Number of subjects. Overwritten by `nrow(covariates)` if covariates are provided.
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
  return(df)
}
