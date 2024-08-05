# Purpose: Functions for hazard ratio analysis.
# Updated: 2022-03-23

#' Calculate Hazard Ratio
#'
#' Calculate the hazard ratio, comparing two treatment arms, its confidence
#' interval and p-value. Also provides the Schoenfeld residual test of the
#' proportional hazards assumption.
#'
#' @param data Reconstructed data.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
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
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Fit Cox model.
  cox_model <- survival::coxph(
    survival::Surv(time, status) ~ arm, data = df)
  cox_model_summary <- summary(cox_model, conf.int = 0.95)
  
  # Extract hazard ratio.
  hazard_ratio <- data.frame(
    "hr" = cox_model_summary$conf.int[1],
    "lower" = cox_model_summary$conf.int[3],
    "upper" = cox_model_summary$conf.int[4],
    "p" = cox_model_summary$coefficients[5]
  )
  
  # Schoefeld residual test.
  schoefeld_test <- survival::cox.zph(cox_model)
  hazard_ratio$schoenfeld_test <- schoefeld_test$table[1, 3]
  return(hazard_ratio)
}
