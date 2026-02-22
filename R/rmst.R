# Purpose: RMST-related functions.
# Updated: 2024-01-08

#' RMST Influence
#'
#' Calculates the influence function contribution for the restricted mean
#' survival time (RMST) at \code{tau} for each observation.
#'
#' @param data Data.frame with time and status (0 = censored, 1 = event).
#' @param tau Truncation time at which the RMST is evaluated. Defaults to the
#'   maximum observation time if NULL or if tau exceeds it.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with the same rows as \code{data}, plus column \code{influence}.
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
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
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
