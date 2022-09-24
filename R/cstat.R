# Purpose: C-statistic calculation.
# Updated: 2022-09-24

#' Calculate Concordance
#'
#' Calculate the C-statistic in the presence of censoring using the method of:
#' <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3079915/>.
#'
#' @param train_data Data for estimating the censoring distribution. Should
#'   include status and time.
#' @param test_data Data for evaluating the risk score. Should include risk,
#'   stats, and time.
#' @param risk_name Name of the risk column.
#' @param status_name Name of status column.
#' @param tau Truncation time.
#' @param time_name Name of time column.
#' @return Data.frame.
#' 
#' @importFrom dplyr "%>%"
#' @export
Cstat <- function(
  train_data,
  test_data,
  risk_name = "risk",
  status_name = "status",
  tau = NULL,
  time_name = "time"
) {
  if (is.null(tau)) {
    tau <- stats::quantile(train_data$time, 0.95)
  }
  
  # Estimate censoring distribution.
  status <- time <- NULL
  df0 <- train_data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    dplyr::select(status, time) %>%
    dplyr::mutate(status = 1 - status)
  df0_curves <- GetCurves(df0)
  km <- df0_curves@Surv
  
  # Evaluate risk score.
  risk <- NULL
  df1 <- test_data %>%
    dplyr::rename(
      risk = {{risk_name}},
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    dplyr::select(risk, status, time)
  
  df1_event <- df1 %>% dplyr::filter((status == 1) & (time < tau))
  df1_event_split <- split(df1_event, seq_len(nrow(df1_event)))
  results <- lapply(df1_event_split, function(x) {
    
    ti <- as.numeric(x$time)
    ri <- as.numeric(x$risk)
    
    p_cens <- km(x$time)
    weight <- 1 / (p_cens^2)
    
    # Indicator Ti < Tj.
    t1_lt_t2 <- 1 * (ti < df1$time)
    
    # Indicator Ri > Rj.
    r1_gt_r2 <- 1 * (ri > df1$risk)
    
    # Calculation.
    numerator <- weight * sum(t1_lt_t2 * r1_gt_r2)
    denominator <- weight * sum(t1_lt_t2)
    return(c(numerator, denominator))
    
  })
  results <- do.call(rbind, results)
  out <- sum(results[, 1]) / sum(results[, 2])
  return(out)
}
