# Purpose: Compare survival rates.
# Updated: 2022-03-23

# -----------------------------------------------------------------------------
# Risk difference, ratio, odds ratio.
# -----------------------------------------------------------------------------

#' Calculate Rate Difference.
#' 
#' Calculate the rate difference comparing two groups.
#' Input data.frame should contain arm, taking values 0 and 1,
#' the event rate, and the standard error.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame.
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
      p = 2 * stats::pnorm(
        q = abs(est) / se,
        lower.tail = FALSE
      )
    )
  return(rd)
}


#' Calculate Rate Ratio
#' 
#' Calculate the rate difference comparing two groups.
#' Input data.frame should contain arm, taking values 0 and 1,
#' the event rate, and the standard error.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame
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
      p = 2 * stats::pnorm(
        q = abs(log(est)) / log_se,
        lower.tail = FALSE
      )
    ) %>% 
    dplyr::select(-log_se)
  return(rr)
}


#' Calculate Rate Odds Ratio
#' 
#' Calculate the odds ratio comparing two treatment arms.
#' Input data.frame should contain arm, taking values 0 and 1,
#' the event rate, and the standard error.
#' 
#' @param rates Data.frame.
#' @param alpha Type I error.
#' @param arm_name Name of arm column.
#' @param rate_name Name of rate column.
#' @param se_name Name of standard error column.
#' @return Data.frame
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
