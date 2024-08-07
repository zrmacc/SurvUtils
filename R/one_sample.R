# Purpose: One-sample analysis.
# Updated: 2024-08-04

# -----------------------------------------------------------------------------

#' One Sample Rates
#' 
#' @param data Data.frame.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
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
      status = {{status_name}},
      time = {{time_name}}
    )
  if (is.null(tau)) {
    tau <- max(df$time)
  }

  km <- SurvCurves(df, alpha = alpha)
  out <- data.frame(
    tau = tau,
    rate = km@Surv(tau),
    se = sqrt(km@SurvVar(tau)),
    lower = km@SurvLower(tau),
    upper = km@SurvUpper(tau)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' One Sample Cumulative Incidence Curve
#' 
#' For calculating cumulative incidence in the presence in the presence of a 
#' competing risk, status is assumed to have values of 0 for censoring, 1 for
#' the event of interest, and 2 for the competing risk.
#' 
#' @param data Data.frame. 
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
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
      status = {{status_name}},
      time = {{time_name}}
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


#' One Sample Percentiles
#' 
#' @param data Data.frame.
#' @param p Percentile. 
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
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
      status = {{status_name}},
      time = {{time_name}}
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


#' One Sample RMST
#' 
#' @param data Data.frame.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
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
      status = {{status_name}},
      time = {{time_name}}
    )
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  
  # Table.
  time <- NULL
  tab <- TabulateKM(df, alpha = alpha) %>%
    dplyr::filter(time <= tau)
  auc <- RMST(status = df$status, time = df$time, tau = tau)
  
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
