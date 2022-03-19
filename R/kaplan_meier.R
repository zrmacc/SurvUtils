#' Purpose: Kaplan-Meier estimator.
#' Updated: 2022-03-19

#' Tabulate Kaplan-Meier 
#' 
#' @param data Data.frame.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param alpha Type I error, for confidence intervals.
#' @return Data.frame.
#' @importFrom dplyr "%>%"
#' @export

TabulateKM <- function(
    data,
    status_name = "status",
    time_name = "time",
    alpha = 0.05
) {
  
  # Format data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Events table.
  out <- df %>%
    dplyr::arrange(time) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      censor = sum(status == 0),
      events = sum(status == 1)
    ) %>% dplyr::ungroup()
  
  # Add initial row.
  out <- rbind(c(0.00, 0, 0), out)
  
  # Number at risk.
  n <- nrow(data)
  n_removed <- cumsum(out$censor + out$events)
  out$nar <- n - c(0, n_removed[1:(nrow(out) - 1)])
  
  # Critical value.
  z <- stats::qnorm(p = 1 - alpha / 2)
  
  # Hazard.
  out$haz <- out$events / out$nar
  out$cum_haz <- cumsum(out$haz)
  out$var_cum_haz <- cumsum(1 / out$nar^2)
  out$cum_haz_lower <- out$cum_haz * 
    exp(-z * sqrt(out$var_cum_haz) / out$cum_haz)
  out$cum_haz_upper <- out$cum_haz *
    exp(+z * sqrt(out$var_cum_haz) / out$cum_haz)
  
  # Survival.
  out$surv <- cumprod(1 - out$haz)
  out$var_surv <- (out$surv^2) * out$var_cum_haz
  out$surv_lower <- out$surv ^ exp(
    -z * sqrt(out$var_surv) / (out$surv * log(out$surv))
  )
  out$surv_upper <- out$surv ^ exp(
    +z * sqrt(out$var_surv) / (out$surv * log(out$surv))
  )
  return(out)
}
