#' Purpose: Kaplan-Meier estimator.
#' Updated: 2023-09-26


#' Tabulate Kaplan-Meier
#'
#' Builds a table of unique event/censoring times with counts, number at risk,
#' hazard, cumulative hazard, survival, and pointwise confidence intervals.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @return Data.frame with columns \code{time}, \code{censor}, \code{events},
#'   \code{nar}, \code{haz}, \code{cum_haz}, \code{cum_haz_var}, \code{cum_haz_lower},
#'   \code{cum_haz_upper}, \code{surv}, \code{surv_var}, \code{surv_lower}, \code{surv_upper}.
#' @export
TabulateKM <- function(
    data,
    status_name = "status",
    time_name = "time",
    alpha = 0.05
) {
  
  # Format data.
  time <- status <- NULL
  df <- data %>%
    dplyr::rename(
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
    )
  
  # Events table.
  out <- df %>%
    dplyr::arrange(time) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      censor = sum(status == 0),
      events = sum(status == 1)
    ) %>% dplyr::ungroup()
  
  # Add initial row.
  out <- rbind(c(0.00, 0, 0), out)
  
  # Number at risk.
  n <- nrow(data)
  n_removed <- cumsum(out$censor + out$events)
  out$nar <- n - c(0, n_removed[1:(nrow(out) - 1)])
  
  # Critical value.
  z <- stats::qnorm(p = 1 - alpha / 2)
  
  # Hazard.
  out$haz <- out$events / out$nar
  out$cum_haz <- cumsum(out$haz)
  out$cum_haz_var <- cumsum(out$events / out$nar^2)
  
  out$cum_haz_lower <- out$cum_haz * 
    exp(-z * sqrt(out$cum_haz_var) / out$cum_haz)
  out$cum_haz_upper <- out$cum_haz *
    exp(+z * sqrt(out$cum_haz_var) / out$cum_haz)
  
  out$cum_haz_lower <- ifelse(out$cum_haz == 0.0, 0, out$cum_haz_lower)
  out$cum_haz_upper <- ifelse(out$cum_haz == 0.0, 0, out$cum_haz_upper)
  
  # Survival.
  out$surv <- cumprod(1 - out$haz)
  out$surv_var <- (out$surv^2) * out$cum_haz_var
  
  out$surv_lower <- out$surv ^ exp(
    -z * sqrt(out$surv_var) / (out$surv * log(out$surv))
  )
  out$surv_upper <- out$surv ^ exp(
    +z * sqrt(out$surv_var) / (out$surv * log(out$surv))
  )
  
  out$surv_lower <- ifelse(out$surv == 1.0, 1, out$surv_lower)
  out$surv_upper <- ifelse(out$surv == 1.0, 1, out$surv_upper)
  
  out$surv_lower <- ifelse(out$surv == 0.0, 0, out$surv_lower)
  out$surv_upper <- ifelse(out$surv == 0.0, 0, out$surv_upper)
  
  return(out)
}


# -----------------------------------------------------------------------------


#' Generate Survival and Hazard Curves
#'
#' Fits the Kaplan-Meier estimator and returns an object containing step
#' functions for the cumulative hazard, survival, their variances, and
#' confidence bounds, plus number at risk. Intended for a single sample.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @importFrom stats stepfun
#' @return Object of class \code{OneSampleSurv} with slots \code{CumHaz}, \code{Surv},
#'   \code{NAR}, their variance and CI functions, and \code{tmax}.
#' @export
SurvCurves <- function(
    data,
    alpha = 0.05,
    status_name = "status",
    time_name = "time"
) {
  
  # Format data.
  df <- data %>%
    dplyr::rename(
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
    )
  
  # Tabulate Kaplan-Meier.
  km <- TabulateKM(df, alpha = alpha)
  
  # Define step functions.
  fn <- list()
  fn$cum_haz <- stepfun(x = km$time, y = c(0, km$cum_haz))
  fn$cum_haz_var <- stepfun(x = km$time, y = c(0, km$cum_haz_var))
  fn$cum_haz_lower <- stepfun(x = km$time, y = c(0, km$cum_haz_lower))
  fn$cum_haz_upper <- stepfun(x = km$time, y = c(0, km$cum_haz_upper))
  fn$nar <- stepfun(x = km$time, y = c(km$nar, 0), right = TRUE)
  fn$surv <- stepfun(x = km$time, y = c(1, km$surv))
  fn$surv_var <- stepfun(x = km$time, y = c(0, km$surv_var))
  fn$surv_lower <- stepfun(x = km$time, y = c(1, km$surv_lower))
  fn$surv_upper <- stepfun(x = km$time, y = c(1, km$surv_upper))
  for (i in 1:length(fn)) {
    f <- fn[[i]]
    class(f) <- "function"
    fn[[i]] <- f
  }
  
  # Output.
  out <- methods::new(
    "OneSampleSurv",
    CumHaz = fn[["cum_haz"]],
    CumHazVar = fn[["cum_haz_var"]],
    CumHazLower = fn[["cum_haz_lower"]],
    CumHazUpper = fn[["cum_haz_upper"]],
    NAR = fn[["nar"]],
    Surv = fn[["surv"]],
    SurvVar = fn[["surv_var"]],
    SurvLower = fn[["surv_lower"]],
    SurvUpper = fn[["surv_upper"]],
    tmax = max(km$time)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Kaplan-Meier Influence
#'
#' Calculates the influence function contribution for the Kaplan-Meier
#' survival estimate at \code{tau} for each observation. Used for standard
#' error estimation and the perturbation bootstrap.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param tau Time at which the survival probability is evaluated. Defaults to
#'   the maximum observation time if NULL or if tau exceeds it.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with the same rows as \code{data}, plus column \code{influence}.
#' @export
KMInfluence <- function(
    data,
    tau = NULL,
    status_name = "status",
    time_name = "time"
) {
  
  # Format data.
  time <- status <- NULL
  df <- data %>%
    dplyr::rename(
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
    )
  
  # Evaluation time.
  tmax <- max(df$time)
  if (is.null(tau) || tau > tmax) {
    if (!is.null(tau) && tau > tmax) {
      warning("tau cannot exceed the maximum observation time.")
    }
    tau <- max(df$time)
  }
  
  # Calculate influence function.
  influence <- InfluenceKM(
    status = df$status,
    time = df$time,
    trunc_time = tau
  )
  
  data$influence <- as.numeric(influence)
  return(data)
}

