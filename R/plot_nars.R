# Purpose: Functions for plotting numbers at risk.
# Updated: 2021-05-01

# -----------------------------------------------------------------------------
# Numbers at risk.
# -----------------------------------------------------------------------------

#' Get Number at Risk Curve
#' 
#' Return a function that calculates the number at risk for a single treatment
#' arm.
#' 
#' @param data Data.frame.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Step function.
#' 
#' @importFrom dplyr "%>%" rename
#' @importFrom stats stepfun
#' @importFrom survival survfit Surv
#' @export 

GetNARCurve <- function(
  data, 
  status_name = "status", 
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  fit <- survival::survfit(Surv(time, status) ~ 1, data = df)
  g <- stats::stepfun(x = fit$time, y = c(nrow(df), fit$n.risk))
  return(g)
}


#' Get Numbers at Risk
#' 
#' Numbers at risk for competing risks data.
#' 
#' @param data Data.frame.
#' @param x_breaks Time points at which to determine the NARs.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame containing `time`, `nar_ctrl`, `nar_trt`.
#' 
#' @importFrom dplyr "%>%" filter rename

GetKMNARs <- function(
  data, 
  x_breaks, 
  arm_name = "arm",
  status_name = "status",
  time_name = "time"
) {
  
  # Data prep.
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # NAR functions.
  arm <- NULL
  g0 <- data %>% 
    dplyr::filter(arm == 0) %>%
    GetNARCurve()
  g1 <- data %>%
    dplyr::filter(arm == 1) %>%
    GetNARCurve()
  
  # Output.
  out <- data.frame(
    time = x_breaks,
    nar_ctrl = g0(x_breaks),
    nar_trt = g1(x_breaks)
  )
  return(out)
}


# -----------------------------------------------------------------------------
# Plot numbers at risk.
# -----------------------------------------------------------------------------

#' Plot Numbers at Risk
#' 
#' @param data Data.frame.
#' @param x_breaks X-axis breaks.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param x_labs X-axis labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis name.
#' @param y_labs Y-axis labels.
#' @return ggplot.
#' 
#' @importFrom dplyr "%>%" rename mutate
#' @importFrom ggplot2 aes geom_text ggplot scale_x_continuous 
#'   scale_y_discrete theme theme_bw 
#' @importFrom tidyr pivot_longer
#' @export

PlotNARs <- function(
  data,
  x_breaks,
  arm_name = "arm",
  status_name = "status",
  time_name = "time",
  x_labs = NULL,
  x_max = NULL,
  x_name = NULL,
  y_labs = c("Ctrl", "Trt")
) {
  
  # Defaults.
  if (is.null(x_labs)) {
    x_labs = x_breaks
  }
  if (is.null(x_max)) {
    x_max = max(x_breaks)
  }
  
  # Data prep.
  nar_ctrl <- NULL
  nar_trt <- NULL
  df <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    GetKMNARs(x_breaks) %>%
    tidyr::pivot_longer(
      cols = c(nar_ctrl, nar_trt),
      names_to = "arm",
      values_to = "nar"
    ) %>%
    dplyr::mutate(
      arm = factor(arm, c("nar_ctrl", "nar_trt"), y_labs)
    )
  
  # Plotting.
  arm <- NULL
  nar <- NULL
  time <- NULL
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ggplot2::geom_text(
      aes(x = time, y = arm, label = nar)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) + 
    ggplot2::scale_y_discrete(
      name = NULL,
      labels = y_labs
    )
  return(q)
}