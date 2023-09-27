#' Purpose: Kaplan-Meier estimator.
#' Updated: 2023-09-26


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
  time <- status <- NULL
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
  out$cum_haz_var <- cumsum(out$events / out$nar^2)
  
  out$cum_haz_lower <- out$cum_haz * 
    exp(-z * sqrt(out$cum_haz_var) / out$cum_haz)
  out$cum_haz_upper <- out$cum_haz *
    exp(+z * sqrt(out$cum_haz_var) / out$cum_haz)
  
  out$cum_haz_lower <- ifelse(out$cum_haz == 0.0, 0, out$cum_haz_lower)
  out$cum_haz_upper <- ifelse(out$cum_haz == 0.0, 0, out$cum_haz_upper)
  
  # Survival.
  out$surv <- cumprod(1 - out$haz)
  out$surv_var <- (out$surv^2) * out$cum_haz_var
  
  out$surv_lower <- out$surv ^ exp(
    -z * sqrt(out$surv_var) / (out$surv * log(out$surv))
  )
  out$surv_upper <- out$surv ^ exp(
    +z * sqrt(out$surv_var) / (out$surv * log(out$surv))
  )
  
  out$surv_lower <- ifelse(out$surv == 1.0, 1, out$surv_lower)
  out$surv_upper <- ifelse(out$surv == 1.0, 1, out$surv_upper)
  
  out$surv_lower <- ifelse(out$surv == 0.0, 0, out$surv_lower)
  out$surv_upper <- ifelse(out$surv == 0.0, 0, out$surv_upper)
  
  return(out)
}


# -----------------------------------------------------------------------------


#' Generate Survival and Hazard Curves
#' 
#' Intended for data from a single sample.
#' 
#' @param data Data.frame.
#' @param alpha Type I error.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @importFrom stats stepfun
#' @export
#' @return stepfun.

GetCurves <- function(
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
  
  # Tabulate Kaplan-Meier.
  km <- TabulateKM(df, alpha = alpha)
  
  # Define step functions.
  fn <- list()
  fn$cum_haz <- stepfun(x = km$time, y = c(0, km$cum_haz))
  fn$cum_haz_var <- stepfun(x = km$time, y = c(0, km$cum_haz_var))
  fn$cum_haz_lower <- stepfun(x = km$time, y = c(0, km$cum_haz_lower))
  fn$cum_haz_upper <- stepfun(x = km$time, y = c(0, km$cum_haz_upper))
  fn$nar <- stepfun(x = km$time, y = c(km$nar, 0), right = TRUE)
  fn$surv <- stepfun(x = km$time, y = c(1, km$surv))
  fn$surv_var <- stepfun(x = km$time, y = c(0, km$surv_var))
  fn$surv_lower <- stepfun(x = km$time, y = c(1, km$surv_lower))
  fn$surv_upper <- stepfun(x = km$time, y = c(1, km$surv_upper))
  for (i in 1:length(fn)) {
    f <- fn[[i]]
    class(f) <- "function"
    fn[[i]] <- f
  }
  
  # Output.
  out <- methods::new(
    "OneSample",
    CumHaz = fn[["cum_haz"]],
    CumHazVar = fn[["cum_haz_var"]],
    CumHazLower = fn[["cum_haz_lower"]],
    CumHazUpper = fn[["cum_haz_upper"]],
    NAR = fn[["nar"]],
    Surv = fn[["surv"]],
    SurvVar = fn[["surv_var"]],
    SurvLower = fn[["surv_lower"]],
    SurvUpper = fn[["surv_upper"]],
    tmax = max(km$time)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Kaplan-Meier Influence
#' 
#' Calculate Kaplan-Meier influence function for each observation.
#' 
#' @param data Data.frame.
#' @param tau Truncation time at which to calculate the influence.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @export

KMInfluence <- function(
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
  
  influence <- InfluenceKM(
    status = df$status,
    time = df$time,
    trunc_time = tau
  )
  return(as.numeric(influence))
}

