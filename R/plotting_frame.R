# Purpose: Functions for plotting Kaplan-Meier curves.
# Updated: 2021-05-01

#' Prepare Kaplan-Meier Plotting Frame
#' 
#' @param data Data.frame
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' 
#' @importFrom dplyr "%>%" rename

GetKMPlotFrame <- function(
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
  g <- GetKMCurve(df, return_surv)
  
  # Time grid.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(
    time = times,
    prob = g(times)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Get 2 Arm Plotting Frame
#' 
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' 
#' @importFrom dplyr "%>%" filter mutate rename

Get2ArmPlotFrame <- function(
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
    GetKMPlotFrame(
      eval_points = eval_points,
      return_surv = return_surv,
      tau = tau
    ) %>%
    dplyr::mutate(arm = 0)
  
  df1 <- data %>%
    dplyr::filter(arm == 1) %>%
    GetKMPlotFrame(
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

#' Plotting Frame for Parametric Models
#'
#' @param para_model Fitted parametric model.
#' @param tau Upper limit of observation period.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @return Data.frame containing `time`, `surv`, `arm`.

Get2ArmParaPlotFrame <- function(
  para_model, 
  tau,
  eval_points = 1000,
  return_surv = TRUE
) {
  
  # Time grid.
  times <- seq(from = 0, to = tau, length.out = eval_points)
  
  # Survival functions.
  g0 <- para_model@Model0@S
  g1 <- para_model@Model1@S
  
  # Plotting frame
  df0 <- cbind("time" = times, "prob" = g0(times), "arm" = 0)
  df1 <- cbind("time" = times, "prob" = g1(times), "arm" = 1)
  df <- data.frame(rbind(df1, df0))
  df$arm <- factor(df$arm, levels = c(0, 1), ordered = TRUE)
  if (!return_surv) {
    df$prob <- 1 - df$prob
  }
  return(df)
}
