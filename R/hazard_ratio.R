# Purpose: Functions for hazard ratio analysis.
# Updated: 2022-03-23

#' Calculate Hazard Ratio
#'
#' Fits a Cox model with a single arm covariate and returns the hazard ratio
#' (arm 1 vs 0) with confidence interval and p-value. Also returns the
#' Schoenfeld residual test for the proportional hazards assumption.
#'
#' @param data Data.frame with arm (0 and 1), time, and status (0 = censored, 1 = event).
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with one row: \code{hr}, \code{se_log_hr}, \code{lower}, \code{upper},
#'   \code{p}, and \code{schoenfeld_test} (p-value for proportional hazards).
#' @export

CalcHR <- function(
  data,
  arm_name = "arm",
  status_name = "status",
  time_name = "time"
) {
  
  # Format data.
  df <- data %>%
    dplyr::rename(
      arm = dplyr::all_of(arm_name),
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
    )
  
  # Fit Cox model.
  cox_model <- survival::coxph(
    survival::Surv(time, status) ~ arm, data = df)
  cox_model_summary <- summary(cox_model, conf.int = 0.95)
  
  # Extract hazard ratio.
  hazard_ratio <- data.frame(
    "hr" = cox_model_summary$conf.int[1],
    "se_log_hr" = cox_model_summary$coefficients[3],
    "lower" = cox_model_summary$conf.int[3],
    "upper" = cox_model_summary$conf.int[4],
    "p" = cox_model_summary$coefficients[5]
  )
  
  # Schoenfeld residual test.
  schoenfeld_test <- survival::cox.zph(cox_model)
  hazard_ratio$schoenfeld_test <- schoenfeld_test$table[1, 3]
  return(hazard_ratio)
}
