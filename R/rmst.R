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
    tau,
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
  
  inf <- InfluenceRMST(
    status = df$status,
    time = df$time,
    trunc_time = tau
  )
  return(as.numeric(inf))
}
