#' Purpose: Cumulative incidence estimator.
#' Updated: 2024-08-04


#' Tabulate Cumulative Incidence
#'
#' Builds a table of the cumulative incidence curve for the event of interest
#' in the presence of a competing risk. Status: 0 = censoring, 1 = event of
#' interest, 2 = competing risk (e.g. death).
#'
#' @param data Data.frame with time and status.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @return Data.frame with columns including \code{time}, \code{nar}, \code{cic_event},
#'   \code{var_cic_event}, \code{se_cic_event}, \code{cic_event_lower}, \code{cic_event_upper}.
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
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
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
#' Fits the cumulative incidence curve for the event of interest in the presence
#' of a competing risk and returns step functions for the CIC, variance, and
#' confidence bounds. Status: 0 = censoring, 1 = event of interest, 2 = competing risk.
#'
#' @param data Data.frame with time and status.
#' @param alpha Type I error level for confidence intervals (default 0.05).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @importFrom stats stepfun
#' @return Object of class \code{OneSampleCIC} with slots \code{CIC}, \code{CICVar},
#'   \code{CICLower}, \code{CICUpper}, \code{NAR}, and \code{tmax}.
#' @export
CICurves <- function(
    data,
    alpha = 0.05,
    status_name = "status",
    time_name = "time"
) {
  
  # Format data.
  df <- data %>%
    dplyr::rename(
      status = dplyr::all_of(status_name),
      time = dplyr::all_of(time_name)
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
#' Calculates the influence function contribution for the cumulative incidence
#' of the event of interest at \code{tau} for each observation. Status: 0 =
#' censored, 1 = event of interest, 2 = competing risk.
#'
#' @param data Data.frame with time and status.
#' @param tau Time at which the cumulative incidence is evaluated. Defaults to
#'   the maximum observation time if NULL or if tau exceeds it.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with the same rows as \code{data}, plus column \code{influence}.
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
  influence <- InfluenceCIC(
    status = df$status,
    time = df$time,
    trunc_time = tau
  )
  
  data$influence <- as.numeric(influence)
  return(data)
}
