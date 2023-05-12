# Purpose: Functions for plotting Kaplan-Meier curves.
# Updated: 2021-05-01

#' One Sample Survival Plotting Frame
#' 
#' @param data Data.frame
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @noRd

OneSampleSurvFrame <- function(
  data,
  eval_points = 1000,
  return_surv = TRUE, 
  status_name = "status", 
  tau = NULL,
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  km <- GetCurves(df)
  
  # Time grid.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(
    time = times,
    prob = km@Surv(times),
    lower = km@SurvLower(times),
    upper = km@SurvUpper(times)
  )
  if (!return_surv) {
    out$prob <- 1 - out$prob
    lower <- pmax(1 - out$upper, 0)
    upper <- pmin(1 - out$lower, 1)
    out$lower <- lower
    out$upper <- upper
  }
  return(out)
}


# -----------------------------------------------------------------------------


#' Two Sample Survival Plotting Frame
#' 
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @noRd

TwoSampleSurvFrame <- function(
  data,
  tau,
  arm_name = "arm",
  eval_points = 1000,
  return_surv = TRUE, 
  status_name = "status", 
  time_name = "time"
) {
  
  # Data.frame.
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Prepare data.
  arm <- NULL
  df0 <- data %>%
    dplyr::filter(arm == 0) %>%
    OneSampleSurvFrame(
      eval_points = eval_points,
      return_surv = return_surv,
      tau = tau
    ) %>%
    dplyr::mutate(arm = 0)
  
  df1 <- data %>%
    dplyr::filter(arm == 1) %>%
    OneSampleSurvFrame(
      eval_points = eval_points,
      return_surv = return_surv,
      tau = tau
    ) %>%
    dplyr::mutate(arm = 1)
  
  out <- rbind(df0, df1)
  out$arm <- factor(out$arm, levels = c(0, 1), ordered = TRUE)
  return(out)
}


# -----------------------------------------------------------------------------

#' Plotting Frame for One Sample Parametric Model
#'
#' @param fit Object of class "fit" from \code{Temporal}.
#' @param tau Upper limit of observation period.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @return Data.frame.
#' @noRd

OneSampleModelFrame <- function(
    fit, 
    tau,
    eval_points = 1000,
    return_surv = TRUE
) {
  
  if (!methods::is(fit, "fit")) {
    stop("Parametric model should be an object of class fit from Temporal.")
  }
  
  # Time grid.
  times <- seq(from = 0, to = tau, length.out = eval_points)
  
  # Survival functions.
  Surv <- fit@S
  
  # Plotting frame
  out <- data.frame(
    time = times,
    prob = Surv(times)
  )
  
  if (!return_surv) {
    out$prob <- 1 - out$prob
  }
  
  return(out)
}

# -----------------------------------------------------------------------------

#' Plotting Frame for Parametric Models
#'
#' @param contrast Object of class "contrast" from \code{Temporal}.
#' @param tau Upper limit of observation period.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @return Data.frame.
#' @noRd

TwoSampleModelFrame <- function(
  contrast, 
  tau,
  eval_points = 1000,
  return_surv = TRUE
) {
  
  # Arm 0.
  df0 <- OneSampleModelFrame(
    fit = contrast@Model0, 
    tau = tau,
    eval_points = eval_points,
    return_surv = return_surv
  ) %>% dplyr::mutate(arm = 0)
  
  # Arm 1.
  df1 <- OneSampleModelFrame(
    fit = contrast@Model1, 
    tau = tau,
    eval_points = eval_points,
    return_surv = return_surv
  ) %>% dplyr::mutate(arm = 1)
  
  # Plotting frame
  out <- data.frame(rbind(df1, df0))
  out$arm <- factor(out$arm, levels = c(0, 1), ordered = TRUE)
  return(out)
}


# -----------------------------------------------------------------------------

#' One Sample Cumulative Hazard Frame
#' 
#' @param data Data.frame
#' @param eval_points Number of points at which to evaluate the curve.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @noRd

OneSampleCHFrame <- function(
    data,
    eval_points = 1000,
    status_name = "status", 
    tau = NULL,
    time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  ch <- GetCurves(df)
  
  # Time grid.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(
    time = times,
    cumhaz = ch@CumHaz(times),
    lower = ch@CumHazLower(times),
    upper = ch@CumHazUpper(times)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Two Sample Cumulative Hazard Frame
#' 
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @noRd

TwoSampleCHFrame <- function(
    data,
    tau,
    arm_name = "arm",
    eval_points = 1000,
    status_name = "status", 
    time_name = "time"
) {
  
  # Data.frame.
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Prepare data.
  arm <- NULL
  df0 <- data %>%
    dplyr::filter(arm == 0) %>%
    OneSampleCHFrame(
      eval_points = eval_points,
      tau = tau
    ) %>%
    dplyr::mutate(arm = 0)
  
  df1 <- data %>%
    dplyr::filter(arm == 1) %>%
    OneSampleCHFrame(
      eval_points = eval_points,
      tau = tau
    ) %>%
    dplyr::mutate(arm = 1)
  
  out <- rbind(df0, df1)
  out$arm <- factor(out$arm, levels = c(0, 1), ordered = TRUE)
  return(out)
}
