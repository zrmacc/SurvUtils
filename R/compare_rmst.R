# Purpose: Compare RMSTs.
# Updated: 2022-08-07


#' Compare RMSTs
#' 
#' Compare the restricted mean survival times of two treatment arms.
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param tau Truncation time.
#' @param time_name Name of time column.
#' @return Data.frame.
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
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
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
