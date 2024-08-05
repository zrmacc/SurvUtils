# Purpose: RMST-related functions.
# Updated: 2024-01-08

#' RMST Influence
#' 
#' Calculate RMST influence function for each observation.
#' 
#' @param data Data.frame.
#' @param tau Truncation time at which to calculate the influence.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Numeric vector.
#' @export
RMSTInfluence <- function(
    data,
    tau = NULL,
    status_name = "status",
    time_name = "time"
) {
  
  # Format data.
  time <- status <- NULL
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Evaluation time.
  tmax <- max(df$time)
  if (is.null(tau) || tau > tmax) {
    if (!is.null(tau) && tau > tmax) {
      warning("tau cannot exceed the maximum observation time.")
    }
    tau <- max(df$time)
  }
  
  # Calculate influence function.
  influence <- InfluenceRMST(
    status = df$status,
    time = df$time,
    trunc_time = tau
  )
  
  data$influence <- as.numeric(influence)
  return(data)
}
