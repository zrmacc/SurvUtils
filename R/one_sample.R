# Purpose: One-sample analysis.
# Updated: 2022-03-27

# -----------------------------------------------------------------------------

#' One Sample Rates
#' 
#' @param data Data.frame.
#' @param tau Truncation time.
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
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

  km <- GetCurves(df, alpha = alpha)
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