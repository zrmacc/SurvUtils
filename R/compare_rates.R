# Purpose: Compare survival rates.
# Updated: 2022-03-23

# -----------------------------------------------------------------------------
# Risk difference, ratio, odds ratio.
# -----------------------------------------------------------------------------

#' Calculate Rate Difference
#'
#' Computes the difference in event rates between two groups (arm 1 minus arm 0)
#' with a confidence interval and p-value. Input data.frame should contain arm
#' (0 and 1), event rate, and standard error.
#'
#' @param rates Data.frame with one row per arm.
#' @param alpha Type I error level for the confidence interval (default 0.05).
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame with one row: \code{stat} = "rd", \code{est}, \code{se}, \code{lower}, \code{upper}, \code{p}.
#' @export
RateDiff <- function(
  rates, 
  alpha = 0.05,
  arm_name = "arm",
  rate_name = "rate",
  se_name = "se"
) {
  
  # Initialize.
  arm <- est <- log_se <- rate <- se <- NULL
  z <- stats::qnorm(p = 1 - alpha / 2)
  
  # Rate difference calculation.
  rd <- rates %>%
    dplyr::rename(
      arm = dplyr::all_of(arm_name),
      rate = dplyr::all_of(rate_name),
      se = dplyr::all_of(se_name)
    ) %>%
    dplyr::summarise(
      stat = "rd",
      est = rate[arm == 1] - rate[arm == 0],
      se = sqrt(se[arm == 1]^2 + se[arm == 0]^2)
    ) %>%
    dplyr::mutate(
      lower = est - z * se,
      upper = est + z * se,
      p = 2 * stats::pnorm(
        q = abs(est) / se,
        lower.tail = FALSE
      )
    )
  return(rd)
}


#' Calculate Rate Ratio
#'
#' Calculate the rate ratio comparing two groups. Input data.frame should
#' contain arm (0 and 1), event rate, and standard error.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame with one row: \code{stat} = "rr", \code{est}, \code{se}, \code{lower}, \code{upper}, \code{p}.
#' @export
RateRatio <- function(
  rates, 
  alpha = 0.05,
  arm_name = "arm",
  rate_name = "rate",
  se_name = "se"
) {
  
  # Initialize.
  arm <- est <- log_se <- rate <- se <- NULL
  z <- stats::qnorm(p = 1 - alpha / 2)
  
  # Rate ratio calculation.
  rr <- rates %>%
    dplyr::rename(
      arm = dplyr::all_of(arm_name),
      rate = dplyr::all_of(rate_name),
      se = dplyr::all_of(se_name)
    ) %>%
    dplyr::summarise(
      stat = "rr",
      est = rate[arm == 1] / rate[arm == 0],
      log_se = sqrt(se[arm == 1]^2 / rate[arm == 1]^2 + se[arm == 0]^2 / rate[arm == 0]^2)
    ) %>% 
    dplyr::mutate(
      lower = est * exp(-z * log_se),
      upper = est * exp(+z * log_se),
      se = est * log_se,
      p = 2 * stats::pnorm(
        q = abs(log(est)) / log_se,
        lower.tail = FALSE
      )
    ) %>% 
    dplyr::select(-log_se)
  return(rr)
}


#' Calculate Odds Ratio
#'
#' Computes the odds ratio comparing two treatment arms (event odds in arm 1
#' vs arm 0) with confidence interval and p-value. Input data.frame should
#' contain arm (0 and 1), event rate, and standard error.
#'
#' @param rates Data.frame with one row per arm.
#' @param alpha Type I error level for the confidence interval (default 0.05).
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame with one row: \code{stat} = "or", \code{est}, \code{se}, \code{lower}, \code{upper}, \code{p}.
#' @export
OddsRatio <- function(
  rates, 
  alpha = 0.05,
  arm_name = "arm",
  rate_name = "rate",
  se_name = "se"
) {
  
  # Initialize.
  arm <- est <- log_se <- rate <- se <- NULL
  z <- stats::qnorm(p = 1 - alpha / 2)
  
  # Odds ratio calculation.
  or <- rates %>% 
    dplyr::rename(
      arm = dplyr::all_of(arm_name),
      rate = dplyr::all_of(rate_name),
      se = dplyr::all_of(se_name)
    ) %>%
    dplyr::summarise(
      stat = "or",
      est = rate[arm == 1] * (1 - rate[arm == 0]) / rate[arm == 0] / (1 - rate[arm == 1]),
      log_se = sqrt(
        se[arm == 1]^2 / (rate[arm == 1]^2 * (1 - rate[arm == 1])^2) +
          se[arm == 0]^2 / (rate[arm == 0]^2 * (1 - rate[arm == 0])^2)
      )
    ) %>% 
    dplyr::mutate(
      lower = est * exp(-z * log_se),
      upper = est * exp(+z * log_se),
      se = est * log_se,
      p = 2 * stats::pnorm(
        q = abs(log(est)) / log_se,
        lower.tail = FALSE
      )
    ) %>% 
    dplyr::select(-log_se)
  return(or)
}


# -----------------------------------------------------------------------------
# Compare KM Rates
# -----------------------------------------------------------------------------

#' Compare Rates
#'
#' Compares the Kaplan-Meier survival (or event) rates of two treatment arms
#' at time \code{tau}, with risk difference, risk ratio, and odds ratio.
#'
#' @param data Data.frame with arm, time, and status.
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @param arm_name Name of arm column (values 0 and 1).
#' @param return_surv If TRUE, compare survival probabilities; if FALSE, compare
#'   event (failure) probabilities (1 - survival).
#' @param status_name Name of status column.
#' @param tau Time at which to evaluate the rate. Defaults to the maximum
#'   observation time if NULL.
#' @param time_name Name of time column.
#' @return Object of class \code{TwoSample} with slots \code{Marginal} (per-arm
#'   rate and se) and \code{Contrasts} (rd, rr, or with est, se, lower, upper, p).
#' @export
CompareRates <- function(
  data,
  alpha = 0.05,
  arm_name = "arm",
  return_surv = TRUE,
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
  
  # Get KM curves.
  arm <- NULL
  km0 <- df %>% 
    dplyr::filter(arm == 0) %>%
    SurvCurves()
  km1 <- df %>% 
    dplyr::filter(arm == 1) %>%
    SurvCurves()
  
  # Calculate per-arm statistics.
  p0 <- km0@Surv(tau)
  p1 <- km1@Surv(tau)
  se0 <- sqrt(km0@SurvVar(tau))
  se1 <- sqrt(km1@SurvVar(tau))
  
  if (!return_surv) {
    p0 <- 1 - p0
    p1 <- 1 - p1
  }
  
  rates <- data.frame(
    arm = c(0, 1),
    tau = tau,
    rate = c(p0, p1),
    se = c(se0, se1)
  )
  
  # Compare arms.
  rd <- rates %>% RateDiff(alpha = alpha)
  rr <- rates %>% RateRatio(alpha = alpha)
  or <- rates %>% OddsRatio(alpha = alpha)
  
  # Output.
  out <- methods::new(
    "TwoSample",
    Marginal = rates,
    Contrasts = rbind(rd, rr, or)
  )
  return(out)
}
