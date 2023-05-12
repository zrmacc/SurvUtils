# Purpose: Functions for plotting Kaplan-Meier curves.
# Updated: 2023-05-11


# -----------------------------------------------------------------------------
# One Sample Cumulative Hazard.
# -----------------------------------------------------------------------------

#' Plot One Sample Cumulative Hazard
#'
#' @param data Data.frame containing time and status.
#' @param color Color for KM curve.
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

PlotOneSampleCH <- function(
    data,
    color = "#0F9D58",
    status_name = "status",
    tau = NULL,
    time_name = "time",
    title = NULL,
    x_breaks = NULL,
    x_labs = NULL,
    x_name = "Time",
    x_max = NULL,
    y_name = "Cumulactive Hazard",
    y_lim = NULL
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(dplyr::all_of(time_name)))
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  if (is.null(x_breaks)) {
    x_breaks <- round(seq(from = 0.0, to = x_max, length.out = 10))
  }
  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  
  # Prepare data.
  df_ch <- data %>% OneSampleCHFrame(
    tau = tau,
    status_name = status_name,
    time_name = time_name
  )
  
  # Plotting.
  lower <- cumhaz <- time <- upper <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    ) + 
    ggplot2::geom_ribbon(
      data = df_ch,
      ggplot2::aes(x = time, ymin = lower, ymax = upper),
      fill = color,
      alpha = 0.2
    ) + 
    ggplot2::geom_step(
      data = df_ch, 
      ggplot2::aes(x = time, y = cumhaz),
      color = color,
      size = 1
    )
  
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

#' Plot Two Sample Cumulative Hazard
#'
#' @param data Data.frame containing time, status, and arm.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param color_labs Color labels.
#' @param color_ctrl Color for control arm.
#' @param color_trt Color for treatment arm.
#' @param plot_haz_haz Logical, plot hazard vs. hazard? If FALSE, plots hazards vs. time.
#' @param tau Truncation time.
#' @param title Plot title.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_lim X-axis limits.
#' @param x_name X-axis name.
#' @param x_max X-axis upper limit.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot.
#' 
#' @importFrom dplyr "%>%"
#' @export

PlotTwoSampleCH <- function(
  data,
  arm_name = "arm",
  color_labs = c("Ctrl", "Trt"),
  color_ctrl = "#EFC000FF",
  color_trt = "#6385B8",
  plot_haz_haz = TRUE,
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  x_breaks = NULL,
  x_labs = NULL,
  x_lim = NULL,
  x_name = NULL,
  x_max = NULL,
  y_name = NULL,
  y_lim = NULL
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(dplyr::all_of(time_name)))
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  if (is.null(x_breaks)) {
    x_breaks <- round(seq(from = 0.0, to = x_max, length.out = 10))
  }
  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  
  # Prepare data.
  df_ch <- data %>% TwoSampleCHFrame(
    tau = tau,
    arm_name = arm_name,
    status_name = status_name,
    time_name = time_name
  )
  
  cumhaz <- NULL
  if (plot_haz_haz) {
    df_ch_0 <- df_ch %>% 
      dplyr::filter(arm == 0) %>% 
      dplyr::select(time, cumhaz) %>%
      dplyr::rename(x = cumhaz)
    df_ch_1 <- df_ch %>%
      dplyr::filter(arm == 1) %>%
      dplyr::select(cumhaz) %>%
      dplyr::rename(y = cumhaz)
    df_ch <- cbind(df_ch_0, df_ch_1)
  }
  
  # Plotting.
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    )
  
  # Cumulative hazards vs. cumulative hazard version.
  if (plot_haz_haz) {
    
    hr <- CalcHR(
      data = data,
      arm_name = arm_name,
      status_name = status_name,
      time_name = time_name
    )
    
    if (is.null(x_name)) {
      x_name = expression(Lambda[0](t))
    }
    if (is.null(y_name)) {
      y_name = expression(Lambda[1](t))
    }
    
    x <- y <- NULL
    q <- q + 
      ggplot2::geom_abline(
        intercept = 0,
        slope = hr$hr,
        linetype = "dashed",
        color = "gray"
      ) +
      ggplot2::geom_step(
        data = df_ch,
        ggplot2::aes(x = x, y = y),
        color = color_trt
      ) + 
      ggplot2::annotate(
        geom = "text",
        x = -Inf,
        y = Inf,
        label = bquote(HR==.(round(hr$hr, digits = 3))),
        vjust = 2,
        hjust = -0.5
      ) +
      ggplot2::scale_x_continuous(
        name = x_name,
        limits = x_lim
      ) +
      ggplot2::scale_y_continuous(
        name = y_name,
        limits = y_lim
      )
    
  } else {
    
    # Cumulative hazards vs. time version.
    
    if (is.null(x_name)) {
      x_name <- "Time"
    }
    if (is.null(y_name)) {
      y_name <- "Cumulative Hazard"
    }
    if (is.null(x_lim)) {
      x_lim <- c(0, x_max)
    }
    
    arm <- lower <- prob <- time <- upper <- NULL
    q <- q + 
      ggplot2::geom_ribbon(
        data = df_ch,
        ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = arm),
        alpha = 0.2
      ) + 
      ggplot2::geom_step(
        data = df_ch, 
        ggplot2::aes(x = time, y = cumhaz, color = arm),
        size = 1
      ) + 
      ggplot2::scale_color_manual(
        name = NULL, 
        labels = color_labs, 
        values = c(color_ctrl, color_trt)
      ) + 
      ggplot2::scale_fill_manual(
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
      )
  }
  
  # Plot adjustments.
  q <- q +
    ggplot2::ggtitle(
      label = title
    )
  return(q)
}
