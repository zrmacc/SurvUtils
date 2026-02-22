# Purpose: One-sample analysis.
# Updated: 2024-08-04

# -----------------------------------------------------------------------------

#' One Sample Event Rate at a Time Point
#'
#' Estimates the event (failure) probability at time \code{tau}, i.e. 1 - S(tau),
#' with standard error and confidence interval from the Kaplan-Meier estimator.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param tau Time at which to evaluate the rate. Defaults to the maximum
#'   observation time if NULL.
#' @param alpha Type I error level for the confidence interval (default 0.05).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with columns \code{tau}, \code{rate}, \code{se}, \code{lower}, \code{upper}.
#' @export
OneSampleRates <- function(
    data,
    tau = NULL,
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
  if (is.null(tau)) {
    tau <- max(df$time)
  }

  km <- SurvCurves(df, alpha = alpha)
  # Event (failure) rate = 1 - S(tau); CI bounds reversed for 1 - S.
  out <- data.frame(
    tau = tau,
    rate = 1 - km@Surv(tau),
    se = sqrt(km@SurvVar(tau)),
    lower = 1 - km@SurvUpper(tau),
    upper = 1 - km@SurvLower(tau)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' One Sample Cumulative Incidence at a Time Point
#'
#' Estimates the cumulative incidence of the event of interest at time \code{tau}
#' in the presence of a competing risk. Status must be 0 = censored, 1 = event of
#' interest, 2 = competing risk (e.g. death).
#'
#' @param data Data.frame with time and status.
#' @param tau Time at which to evaluate the cumulative incidence. Defaults to
#'   the maximum observation time if NULL.
#' @param alpha Type I error level for the confidence interval (default 0.05).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with columns \code{tau}, \code{rate}, \code{se}, \code{lower}, \code{upper}.
#' @export
OneSampleCIC <- function(
    data,
    tau = NULL,
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
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  
  cic <- CICurves(df, alpha = alpha)
  out <- data.frame(
    tau = tau,
    rate = cic@CIC(tau),
    se = sqrt(cic@CICVar(tau)),
    lower = cic@CICLower(tau),
    upper = cic@CICUpper(tau)
  )
  return(out)
}


# -----------------------------------------------------------------------------


#' Extract Percentile
#'
#' For a survival function, extract the first time the probability falls
#' below the percentile q.
#'
#' @param time Time.
#' @param prob Probability. 
#' @param q Percentile.
#' @return Numeric time.
#' @noRd
GetPercentile <- function(time, prob, q = 0.5) {
  
  # Check tail.
  if (min(prob) > q) {return(Inf)}
  
  # First time at which prob <= percentile.
  out <- time[min(which(prob <= q))]
  return(out)
}


#' One Sample Survival Percentiles
#'
#' Estimates quantiles of the survival distribution (e.g. median survival time)
#' with confidence intervals from the Kaplan-Meier estimator.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param p Percentile(s), e.g. 0.5 for median. Can be a vector.
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with columns \code{prob}, \code{time}, \code{lower}, \code{upper}
#'   (one row per value of \code{p}).
#' @importFrom dplyr "%>%"
#' @export
OneSamplePercentiles <- function(
    data,
    p = 0.5,
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
  
  tab <- TabulateKM(df, alpha = alpha)
  
  results <- lapply(p, function(x) {
    sub <- data.frame(
      prob = x,
      time = GetPercentile(tab$time, tab$surv, q = x),
      lower = GetPercentile(tab$time, tab$surv_lower, q = x),
      upper = GetPercentile(tab$time, tab$surv_upper, q = x)
    )
    return(sub)
  })
  out <- do.call(rbind, results)
  return(out)
}


# -----------------------------------------------------------------------------


#' One Sample Restricted Mean Survival Time
#'
#' Estimates the restricted mean survival time (RMST) up to \code{tau}, i.e. the
#' area under the survival curve from 0 to tau, with standard error and
#' confidence interval.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param tau Restriction time. Defaults to the maximum observation time if NULL.
#' @param alpha Type I error level for the confidence interval (default 0.05).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with columns \code{tau}, \code{auc} (RMST), \code{se},
#'   \code{lower}, \code{upper}.
#' @export
OneSampleRMST <- function(
    data,
    tau = NULL,
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
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  
  # Table.
  time <- NULL
  tab <- TabulateKM(df, alpha = alpha) %>%
    dplyr::filter(time <= tau)
  auc <- RMST(status = df$status, time = df$time, extend = FALSE, tau = tau)
  
  # Variance calculation.
  # Mu_{\tau}(t) = \int_{t}^{\tau}S(u)du.
  # The variance is the sum over *event* times <= tau of
  # mu_t^2 * dN_t / Y_t where mu_t is defined as above,
  # dN_t is the number of events and Y_t is the number at risk.
  
  event_times <- tab$time
  surv <- tab$surv
  
  delta_t <- diff(c(event_times, tau))
  mu_t <- rev(cumsum(rev(delta_t * surv)))
  
  var <- sum(mu_t^2 * tab$haz / tab$nar)
  
  # Output.
  z <- stats::qnorm(p = 1 - alpha / 2)
  out <- data.frame(
    tau = tau,
    auc = auc,
    se = sqrt(var)
  )
  out$lower <- out$auc - z * out$se
  out$upper <- out$auc + z * out$se
  return(out)
}
