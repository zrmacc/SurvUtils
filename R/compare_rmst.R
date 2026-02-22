# Purpose: Compare RMSTs.
# Updated: 2022-08-07

#' Compare RMSTs
#'
#' Compares the restricted mean survival times (RMST) of two treatment arms
#' up to time \code{tau}, with difference and ratio contrasts.
#'
#' @param data Data.frame with arm, time, and status (0 = censored, 1 = event).
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @param arm_name Name of arm column (values 0 and 1).
#' @param status_name Name of status column.
#' @param tau Restriction time. Defaults to the maximum observation time if NULL.
#' @param time_name Name of time column.
#' @return Object of class \code{TwoSample} with slots \code{Marginal} (per-arm
#'   RMST and se) and \code{Contrasts} (difference and ratio with est, se, lower, upper, p).
#' @importFrom dplyr "%>%"
#' @export

CompareRMSTs <- function(
    data,
    alpha = 0.05,
    arm_name = "arm",
    status_name = "status",
    tau = NULL,
    time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      arm = dplyr::all_of(arm_name),
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
    )
  
  # Default truncation.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  
  # Arm 0.
  arm <- NULL
  arm0 <- df %>%
    dplyr::filter(arm == 0) %>%
    OneSampleRMST(tau = tau, alpha = alpha) %>%
    dplyr::mutate(arm = 0)
  
  # Arm 1.
  arm1 <- df %>%
    dplyr::filter(arm == 1) %>%
    OneSampleRMST(tau = tau, alpha = alpha) %>%
    dplyr::mutate(arm = 1) 
  
  # Rates.
  rates <- rbind(arm0, arm1)
  
  # Difference.
  rd <- RateDiff(rates, rate_name = "auc")
  
  # Ratio.
  rr <- RateRatio(rates, rate_name = "auc")
  
  # Output.
  out <- methods::new(
    "TwoSample",
    Marginal = rates,
    Contrasts = rbind(rd, rr)
  )
  return(out)
}
