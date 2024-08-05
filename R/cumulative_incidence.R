#' Purpose: Cumulative incidence estimator.
#' Updated: 2024-08-04


#' Tabulate Cumulative Incidence
#' 
#' Status should be coded as 0 for censoring, 1 for the event of interest, 
#' and 2 for the competing risk (e.g. death).
#' 
#' @param data Data.frame.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param alpha Type I error, for confidence intervals.
#' @return Data.frame.
#' @export
TabulateCIC <- function(
    data,
    status_name = "status",
    time_name = "time",
    alpha = 0.05
) {

  # Format data.
  time <- status <- NULL
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Table.
  out <- CalcCIC(status = df$status, time = df$time)

  # Critical value.
  z <- stats::qnorm(p = 1 - alpha / 2)
  
  # Confidence intervals.
  out$cic_event_lower <- out$cic_event * exp(-z * out$se_cic_event / out$cic_event)
  out$cic_event_upper <- out$cic_event * exp(+z * out$se_cic_event / out$cic_event)
  
  out$cic_event_lower <- ifelse(out$cic_event == 0, 0, out$cic_event_lower)
  out$cic_event_upper <- ifelse(out$cic_event == 0, 0, out$cic_event_upper)
  
  # Output.
  return(out)
}


# -----------------------------------------------------------------------------

#' Generate Cumulative Incidence Curves
#' 
#' Intended for data from a single sample. Status should be coded as 0 for
#' censoring, 1 for the event of interest, and 2 for the competing risk (e.g.
#' death).
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @importFrom stats stepfun
#' @export
#' @return OneSampleSurv object.
CICurves <- function(
    data,
    alpha = 0.05,
    status_name = "status",
    time_name = "time"
) {
  
  # Format data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Tabulate cumulative incidence.
  cic <- TabulateCIC(df, alpha = alpha)
  
  # Define step functions.
  fn <- list()
  fn$cic <- stepfun(x = cic$time, y = c(0, cic$cic_event))
  fn$cic_var <- stepfun(x = cic$time, y = c(0, cic$var_cic_event))
  fn$cic_lower <- stepfun(x = cic$time, y = c(0, cic$cic_event_lower))
  fn$cic_upper <- stepfun(x = cic$time, y = c(0, cic$cic_event_upper))
  fn$nar <- stepfun(x = cic$time, y = c(cic$nar, 0), right = TRUE)
  for (i in 1:length(fn)) {
    f <- fn[[i]]
    class(f) <- "function"
    fn[[i]] <- f
  }
  
  # Output.
  out <- methods::new(
    "OneSampleCIC",
    CIC = fn[["cic"]],
    CICVar = fn[["cic_var"]],
    CICLower = fn[["cic_lower"]],
    CICUpper = fn[["cic_upper"]],
    NAR = fn[["nar"]],
    tmax = max(cic$time)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Cumulative Incidence Influence
#' 
#' Calculate cumulative incidence influence function for each observation.
#' 
#' @param data Data.frame.
#' @param tau Truncation time at which to calculate the influence.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with an additional column, `influence`.
#' @export
CICInfluence <- function(
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
  influence <- InfluenceCIC(
    status = df$status,
    time = df$time,
    trunc_time = tau
  )
  
  data$influence <- as.numeric(influence)
  return(data)
}
