# Purpose: Generate pseudo-values.

#' Generate Pseudo-values
#' 
#' Generates pseudo-values for each patient as the target parameter plus
#' the influence function. 
#' 
#' @param data Data.frame.
#' @param tau Truncation time.
#' @param type Type of pseudo-value, select "prob" for probability and "rmst" for 
#'   restricted mean survival time.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
GenPseudo <- function(
  data,
  tau = NULL,
  type = "prob",
  status_name = "status",
  time_name = "time"
) {
  
  # Check type.
  if (! (type %in% c("prob", "rmst"))) {
    stop("Select type from among: 'prob', 'rmst'.")
  }
  
  # Format data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- max(df$time)
  }
  
  # Calculate parameter.
  if (type == "prob") {
    
    param <- OneSampleRates(df, tau = tau)
    param <- param$rate
    
  } else if(type == "rmst") {
    
    param <- OneSampleRMST(df, tau = tau)
    param <- param$auc
    
  }
  
  # Calculate influence function.
  if (type == "prob") {
    
    inf <- KMInfluence(df, tau = tau)
    
  } else {
    
    inf <- RMSTInfluence(df, tau = tau)
    
  }
  
  # Return pseudo-values.
  pseudo <- param + inf
  
}
