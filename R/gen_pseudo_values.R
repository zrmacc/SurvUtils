# Purpose: Generate pseudo-values.
# Updated: 2024-08-04

#' Generate Pseudo-values
#' 
#' Generates pseudo-values for each patient as the target parameter plus the
#' influence function.
#'
#' @param data Data.frame.
#' @param tau Truncation time.
#' @param type Type of pseudo-value, select "prob" for probability and "rmst"
#'   for restricted mean survival time.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame including the influence function and pseudo-value for each
#'   observation.
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
  
  # Calculate parameter.
  if (type == "cic") {
    
    param <- OneSampleCIC(df, tau = tau)
    param <- param$rate
    out <- CICInfluence(df, tau = tau)

  } else if (type == "prob") {
    
    param <- OneSampleRates(df, tau = tau)
    param <- param$rate
    out <- KMInfluence(df, tau = tau)
    
  } else if(type == "rmst") {
    
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
