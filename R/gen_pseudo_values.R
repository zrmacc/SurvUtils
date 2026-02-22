# Purpose: Generate pseudo-values.
# Updated: 2024-08-04

#' Generate Pseudo-values
#'
#' Generates pseudo-values for each observation as the target parameter estimate
#' plus the estimated influence function. Useful for regression modeling and
#' the perturbation bootstrap.
#'
#' @param data Data.frame with survival data (time and status).
#' @param tau Truncation time at which to evaluate the parameter. Defaults to
#'   the maximum observation time if NULL or if tau exceeds it.
#' @param type Type of pseudo-value: \code{"prob"} for event (failure) probability
#'   at tau (1 - S(tau)); \code{"rmst"} for restricted mean survival time up to
#'   tau; \code{"cic"} for cumulative incidence at tau (competing risks; status
#'   0 = censored, 1 = event, 2 = competing risk).
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame with the same rows as \code{data}, plus columns
#'   \code{influence} (estimated influence contribution) and \code{pseudo}
#'   (pseudo-value = parameter estimate + influence).
#' @export 
GenPseudo <- function(
  data,
  tau = NULL,
  type = "prob",
  status_name = "status",
  time_name = "time"
) {
  
  # Format data.
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
  
  # Calculate parameter.
  if (type == "cic") {
    
    param <- OneSampleCIC(df, tau = tau)
    param <- param$rate
    out <- CICInfluence(df, tau = tau)

  } else if (type == "prob") {
    
    param <- OneSampleRates(df, tau = tau)
    param <- param$rate
    out <- KMInfluence(df, tau = tau)
    
  } else if (type == "rmst") {
    
    param <- OneSampleRMST(df, tau = tau)
    param <- param$auc
    out <- RMSTInfluence(df, tau = tau)
    
  } else {
    
    stop("type must be selected from among: cic, prob, rmst")
    
  }
  
  # Output.
  out$pseudo <- param + out$influence
  return(out)
}
