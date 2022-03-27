# Purpose: Functions for plotting Kaplan-Meier curves.
# Updated: 2022-03-25

# -----------------------------------------------------------------------------
# One Sample KM.
# -----------------------------------------------------------------------------

#' Plot Survival Data.
#'
#' @param data Data.frame containing time and status.
#' @param color Color for KM curve.
#' @param fit Parametric fit from the Temporal package.
#' @param plot_surv Logical, TRUE for survival curves, FALSE for cumulative incidence.
#' @param status_name Name of status column.
#' @param tau Truncation time.
#' @param time_name Name of time column.
#' @param title Plot title.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_name X-axis name.
#' @param x_max X-axis upper limit.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot.
#' 
#' @importFrom dplyr "%>%"
#' @export

PlotOneSampleKM <- function(
    data,
    color = "#0F9D58",
    fit = NULL,
    plot_surv = TRUE,
    status_name = "status",
    tau = NULL,
    time_name = "time",
    title = NULL,
    x_breaks = NULL,
    x_labs = NULL,
    x_name = "Time",
    x_max = NULL,
    y_name = "Survival",
    y_lim = c(0, 1)
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(dplyr::all_of(time_name)))
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
  
  # Prepare data.
  df_km <- data %>% OneSamplePlotFrame(
    tau = tau,
    return_surv = plot_surv,
    status_name = status_name,
    time_name = time_name
  )
  
  # Prepare plotting frame for parametric model.
  if (!is.null(fit)) {
    df_para <- OneSampleModelFrame(
      fit = fit,
      tau = tau,
      return_surv = plot_surv
    )
  }
  
  # Plotting.
  lower <- prob <- time <- upper <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    ) + 
    ggplot2::geom_ribbon(
      data = df_km,
      ggplot2::aes(x = time, ymin = lower, ymax = upper),
      fill = color,
      alpha = 0.2
    ) + 
    ggplot2::geom_step(
      data = df_km, 
      ggplot2::aes(x = time, y = prob),
      color = color,
      size = 1
    )
  
  # Add parametric curves.
  if (!is.null(fit)) {
    curve <- NULL
    df0 <- df_km %>% 
      dplyr::select("time", "prob") %>%
      dplyr::mutate(curve = "1km")
    df1 <- df_para %>%
      dplyr::mutate(curve = "2fit")
    df <- rbind(df0, df1)
    
    q <- q + 
      ggplot2::geom_line(
        data = df,
        ggplot2::aes(x = time, y = prob, linetype = curve),
        color = color,
        size = 1
      ) +
      ggplot2::scale_linetype_manual(
        name = "Curve",
        values = c("solid", "dotted"),
        labels = c("KM", "Fit")
      )
  }
  
  # Plot adjustments.
  q <- q + 
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

# -----------------------------------------------------------------------------
# Two Sample KM.
# -----------------------------------------------------------------------------

#' Plot Survival Data.
#'
#' @param data Data.frame containing time, status, and arm.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param color_labs Color labels.
#' @param color_ctrl Color for control arm.
#' @param color_trt Color for treatment arm.
#' @param contrast Fitted parametric contrast from the Temporal package.
#' @param plot_surv Logical, TRUE for survival curves, FALSE for cumulative incidence.
#' @param tau Truncation time.
#' @param title Plot title.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_name X-axis name.
#' @param x_max X-axis upper limit.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot.
#' 
#' @importFrom dplyr "%>%"
#' @export

PlotTwoSampleKM <- function(
  data,
  arm_name = "arm",
  color_labs = c("Ctrl", "Trt"),
  color_ctrl = "#EFC000FF",
  color_trt = "#6385B8",
  contrast = NULL,
  plot_surv = TRUE,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  x_breaks = NULL,
  x_labs = NULL,
  x_name = "Time",
  x_max = NULL,
  y_name = "Survival",
  y_lim = c(0, 1)
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(dplyr::all_of(time_name)))
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
  
  # Prepare data.
  df_km <- data %>% TwoSamplePlotFrame(
    tau = tau,
    arm_name = arm_name,
    return_surv = plot_surv,
    status_name = status_name,
    time_name = time_name
  )
  
  # Prepare plotting frame for parametric model.
  if (!is.null(contrast)) {
    df_para <- TwoSampleModelFrame(
      contrast = contrast,
      tau = tau,
      return_surv = plot_surv
    )
  }
  
  # Plotting.
  arm <- lower <- prob <- time <- upper <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    ) + 
    ggplot2::geom_ribbon(
      data = df_km,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = arm),
      alpha = 0.2
    ) + 
    ggplot2::geom_step(
      data = df_km, 
      ggplot2::aes(x = time, y = prob, color = arm),
      size = 1
    )
  
  # Add parametric curves.
  if (!is.null(contrast)) {
    
    curve <- NULL
    df0 <- df_km %>% 
      dplyr::select("time", "prob", "arm") %>%
      dplyr::mutate(curve = "1km")
    df1 <- df_para %>%
      dplyr::mutate(curve = "2fit")
    df <- rbind(df0, df1)
    
    q <- q + 
      ggplot2::geom_line(
        data = df,
        ggplot2::aes(x = time, y = prob, color = arm, linetype = curve),
        size = 1
      ) +
      ggplot2::scale_linetype_manual(
        name = NULL,
        values = c("solid", "dotted"),
        labels = c("KM", "Fit")
      )
  }
  
  # Plot adjustments.
  q <- q + 
    ggplot2::scale_color_manual(
      name = NULL, 
      labels = color_labs, 
      values = c(color_ctrl, color_trt)
    ) + ggplot2::scale_fill_manual(
      name = NULL, 
      labels = color_labs, 
      values = c(color_ctrl, color_trt)
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
