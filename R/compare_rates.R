# Purpose: Compare survival rates.
# Updated: 2021-05-01

# -----------------------------------------------------------------------------

#' Get Kaplan-Meier Standard Error Curve
#' 
#' Return a function that calculates the standard error of the survival probability
#' for a single treatment arm.
#' 
#' @param data Data.frame.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Step function.
#' 
#' @importFrom dplyr "%>%" rename
#' @importFrom stats stepfun
#' @importFrom survival survfit Surv
#' @export 

GetKMSECurve <- function(
  data, 
  return_surv = TRUE, 
  status_name = "status", 
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  fit <- survival::survfit(Surv(time, status) ~ 1, data = df)
  if (return_surv) {
    g <- stats::stepfun(x = fit$time, y = c(1, fit$std.err))
  } else {
    g <- stats::stepfun(x = fit$time, y = c(0, 1 - fit$std.err))
  }
  return(g)
}


# -----------------------------------------------------------------------------
# Risk difference, ratio, odds ratio.
# -----------------------------------------------------------------------------

#' Calculate Rate Difference.
#' 
#' Calculate the rate difference comparing two groups.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame.
#' 
#' @importFrom dplyr "%>%" mutate rename summarise
#' @importFrom stats pnorm qnorm
#' @export

RateDiff <- function(
  rates, 
  alpha = 0.05,
  arm_name = "arm",
  rate_name = "rate",
  se_name = "se"
) {
  
  # Initialize.
  arm <- NULL
  est <- NULL
  log_se <- NULL
  rate <- NULL
  se <- NULL
  z <- qnorm(p = 1 - alpha / 2)
  
  # Rate difference calculation.
  rd <- rates %>%
    dplyr::rename(
      arm = {{arm_name}},
      rate = {{rate_name}},
      se = {{se_name}}
    ) %>%
    dplyr::summarise(
      stat = "rd",
      est = rate[arm == 1] - rate[arm == 0],
      se = sqrt(se[arm == 1]^2 + se[arm == 0]^2)
    ) %>%
    dplyr::mutate(
      lower = est - z * se,
      upper = est + z * se,
      p = 2 * pnorm(q = abs(est) / se, lower.tail = FALSE)
    )
  return(rd)
}


#' Calculate Rate Ratio
#' 
#' Calculate the rate difference comparing two groups.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame
#' 
#' @importFrom dplyr "%>%" mutate select summarise
#' @importFrom stats pnorm qnorm
#' @export

RateRatio <- function(
  rates, 
  alpha = 0.05,
  arm_name = "arm",
  rate_name = "rate",
  se_name = "se"
) {
  
  # Initialize.
  arm <- NULL
  est <- NULL
  log_se <- NULL
  rate <- NULL
  se <- NULL
  z <- qnorm(p = 1 - alpha / 2)
  
  # Rate ratio calculation.
  rr <- rates %>%
    dplyr::rename(
      arm = {{arm_name}},
      rate = {{rate_name}},
      se = {{se_name}}
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
      p = 2 * pnorm(q = abs(log(est)) / log_se, lower.tail = FALSE)
    ) %>% 
    dplyr::select(
      - log_se
    )
  return(rr)
}


#' Calculate Rate Odds Ratio
#' 
#' Calculate the odds ratio comparing two treatment arms.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame
#' @importFrom dplyr "%>%" mutate select summarise
#' @importFrom stats pnorm qnorm
#' @export

OddsRatio <- function(
  rates, 
  alpha = 0.05,
  arm_name = "arm",
  rate_name = "rate",
  se_name = "se"
) {
  
  # Initialize.
  arm <- NULL
  est <- NULL
  log_se <- NULL
  rate <- NULL
  se <- NULL
  z <- qnorm(p = 1 - alpha / 2)
  
  # Odds ratio calculation.
  or <- rates %>% 
    dplyr::rename(
      arm = {{arm_name}},
      rate = {{rate_name}},
      se = {{se_name}}
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
      p = 2 * pnorm(q = abs(log(est)) / log_se, lower.tail = FALSE)
    ) %>% 
    dplyr::select(
      - log_se
    )
  return(or)
}


# -----------------------------------------------------------------------------
# Compare KM Rates
# -----------------------------------------------------------------------------

#' Compare Rates
#' 
#' Compare the Kaplan-Meier survival or incidence rates of two treatment arms.
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Truncation time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' 
#' @importFrom dplyr "%>%" filter
#' @importFrom methods new
#' @export

CompareKMRates <- function(
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
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Default truncation.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  
  # Get KM curves.
  arm <- NULL
  km0 <- df %>% 
    dplyr::filter(arm == 0) %>%
    GetKMCurve(return_surv = return_surv)
  km1 <- df %>% 
    dplyr::filter(arm == 1) %>%
    GetKMCurve(return_surv = return_surv)
  
  # Get variance calculations.
  se0 <- df %>% 
    dplyr::filter(arm == 0) %>%
    GetKMSECurve(return_surv = return_surv)
  se1 <- df %>% 
    dplyr::filter(arm == 1) %>%
    GetKMSECurve(return_surv = return_surv)
  
  # Calculate per-arm statistics.
  rates <- data.frame(
    arm = c(0, 1),
    tau = tau,
    rate = c(km0(tau), km1(tau)),
    se = c(se0(tau), se1(tau))
  )
  
  # Compare arms.
  rd <- rates %>% RateDiff(alpha = alpha)
  rr <- rates %>% RateRatio(alpha = alpha)
  or <- rates %>% OddsRatio(alpha = alpha)
  
  # Output.
  out <- methods::new(
    "SurvSummary",
    Marginal = rates,
    Contrasts = rbind(rd, rr, or)
  )
  return(out)
}


