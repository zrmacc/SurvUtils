# Purpose: Functions for plotting Kaplan-Meier curves.
# Updated: 2021-05-01

# -----------------------------------------------------------------------------
# Kaplan-Meier Curve
# -----------------------------------------------------------------------------

#' Get Kaplan-Meier Curve
#' 
#' Return a function that calculates the survival probability for a single
#' treatment arm.
#' 
#' @param data Data.frame.
#' @param return_surv Logical, TRUE for survival, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Step function.
#' 
#' @importFrom dplyr "%>%" rename
#' @importFrom stats stepfun
#' @importFrom survival survfit Surv
#' @export 

GetKMCurve <- function(
  data, 
  return_surv = TRUE, 
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
  if (return_surv) {
    g <- stats::stepfun(x = fit$time, y = c(1, fit$surv))
  } else {
    g <- stats::stepfun(x = fit$time, y = c(0, 1 - fit$surv))
  }
  return(g)
}


# -----------------------------------------------------------------------------
# Main KM plotting function.
# -----------------------------------------------------------------------------

#' Plot Survival Data.
#'
#' @param data Reconstructed data.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param color_labs Color labels.
#' @param ctrl_color Color for control arm.
#' @param para_model Fitted parametric contrast from the Temporal package.
#' @param plot_surv Logical, TRUE for survival curves, FALSE for cumulative incidence.
#' @param tau Truncation time.
#' @param title Plot title.
#' @param trt_color Color for treatment arm.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_name X-axis name.
#' @param x_max X-axis upper limit.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot
#' 
#' @importFrom dplyr "%>%"
#' @importFrom ggplot2 aes element_blank geom_step geom_line ggplot ggtitle
#'   scale_color_manual scale_x_continuous scale_y_continuous theme theme_bw 
#' @export

PlotKMCurves <- function(
  data,
  arm_name = "arm",
  color_labs = c("Ctrl", "Trt"),
  ctrl_color = "#EFC000FF",
  para_model = NULL,
  plot_surv = TRUE,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  trt_color = "#6385B8",
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
  
  # Prepare data.
  df_km <- data %>% Get2ArmPlotFrame(
    tau = tau,
    arm_name = arm_name,
    return_surv = plot_surv,
    status_name = status_name,
    time_name = time_name
  )
  
  # Prepare plotting frame for parametric model.
  if (!is.null(para_model)) {
    df_para <- Get2ArmParaPlotFrame(
      para_model = para_model,
      tau = tau,
      return_surv = plot_surv
    )
  }
  
  # Plotting.
  arm <- NULL
  prob <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) + 
    ggplot2::geom_step(
      data = df_km, 
      aes(x = time, y = prob, color = arm),
      size = 1
    )
  
  # Add parametric curves.
  if (!is.null(para_model)) {
    q <- q + 
      ggplot2::geom_line(
        data = df_para,
        aes(x = time, y = prob, color = arm),
        linetype = "dashed",
        size = 1
      ) 
  }
  
  # Plot adjustments.
  q <- q + 
    ggplot2::scale_color_manual(
      name = NULL, 
      labels = color_labs, 
      values = c(ctrl_color, trt_color)
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