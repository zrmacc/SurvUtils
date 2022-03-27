# Purpose: Plot area under the curve.
# Updated: 2022-03-27

# -----------------------------------------------------------------------------
# Plot area under the curve.
# -----------------------------------------------------------------------------

#' Plot Area Under the Kaplan-Meier Curve.
#'
#' @param data Data including time, status, arm.
#' @param color Color.
#' @param label Label for the arm.
#' @param legend_pos Legend position.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Truncation time for shading.
#' @param time_name Name of time column.
#' @param title Plot title.
#' @param which_arm 0 for control, 1 for treatment.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_max X-axis upper limit; may differ from tau.
#' @param x_name X-axis name.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot.
#' 
#' @importFrom dplyr "%>%"
#' @export

PlotOneSampleAUC <- function(
  data,
  color = "#0F9D58",
  label = "Ctrl",
  legend_pos = "top",
  return_surv = TRUE,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  which_arm = 0,
  x_breaks = NULL,
  x_labs = NULL,
  x_name = "Time",
  x_max = NULL,
  y_name = "Survival",
  y_lim = c(0, 1)
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(time_name))
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  if (is.null(x_breaks)) {
    x_breaks <- seq(from = 0.0, to = x_max, length.out = 10)
  }
  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }

  # Data prep.
  df <- OneSamplePlotFrame(
    data,
    return_surv = return_surv,
    status_name = status_name,
    time_name = time_name
  )
  
  # Plotting.
  arm <- NULL
  prob <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::geom_ribbon(
      data = df %>% dplyr::filter(time <= tau),
      ggplot2::aes(x = time, ymin = 0, ymax = prob, fill = label),
      alpha = 0.5
    ) + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_pos
    ) + 
    ggplot2::scale_fill_manual(
      name = NULL, 
      labels = label, 
      values = color
    ) + 
    ggplot2::geom_step(
      data = df, 
      ggplot2::aes(x = time, y = prob), 
      size = 1, 
      color = color
    ) +
    ggplot2::scale_x_continuous(
      name = x_name,
      breaks = x_breaks,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_continuous(
      name = y_name,
      limits = y_lim
    ) + 
    ggplot2::ggtitle(
      label = title
    )
  return(q)
}
