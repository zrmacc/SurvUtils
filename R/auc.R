# Purpose: Plot area under the curve.
# Updated: 2021-05-01

# -----------------------------------------------------------------------------
# Area under the curve.
# -----------------------------------------------------------------------------

#' Plot Area Under the Kaplan-Meier Curve.
#'
#' @param data Data including time, status, arm.
#' @param arm_label Label for the arm.
#' @param arm_name Name of arm column.
#' @param color Color.
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
#' @importFrom dplyr "%>%" filter rename
#' @importFrom ggplot2 aes ggplot geom_ribbon geom_step 
#'   scale_fill_manual theme

PlotAUC <- function(
  data,
  arm_label = "Placebo",
  arm_name = "arm",
  color = "#C65842",
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
  
  # Data prep.
  df <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    dplyr::filter(arm == which_arm) %>%
    GetKMPlotFrame(
      return_surv = return_surv,
      tau = tau
    )
  
  # Plotting.
  arm <- NULL
  prob <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::geom_ribbon(
      data = df,
      aes(x = time, ymin = 0, ymax = prob, fill = arm_label),
      alpha = 0.5
    ) + 
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = legend_pos
    ) + 
    ggplot2::scale_fill_manual(
      name = NULL, 
      labels = arm_label, 
      values = color
    ) + 
    ggplot2::geom_step(
      data = df, 
      aes(x = time, y = prob), 
      size = 1, 
      color = color
    ) +
    ggplot2::scale_x_continuous(
      name = x_name,
      breaks = x_breaks,
      labels = x_labs
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
