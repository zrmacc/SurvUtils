# Purpose: Functions for plotting numbers at risk.
# Updated: 2022-03-27


# -----------------------------------------------------------------------------

#' One Sample NAR Frame
#' 
#' @param data Data.frame.
#' @param x_breaks X-axis breaks.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame.
#' @noRd
OneSampleNARs <- function(
    data,
    x_breaks,
    status_name = "status",
    time_name = "time"
) {
  km <- SurvCurves(
    data,
    status_name = status_name,
    time_name = time_name
  )
  out <- data.frame(
    time = x_breaks,
    nar = km@NAR(x_breaks)
  )
  return(out)
}


# -----------------------------------------------------------------------------
# One sample.
# -----------------------------------------------------------------------------


#' Plot Two Sample Numbers at Risk
#' 
#' @param data Data.frame.
#' @param specified_nars Optional data.frame of specified NARs. Should include
#'   "time" and "nar". Overrides `data` if present.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis name.
#' @param y_lab Y-axis label.
#' @return ggplot.
#' @importFrom dplyr "%>%"
#' @export
PlotOneSampleNARs <- function(
    data,
    specified_nars = NULL,
    status_name = "status",
    time_name = "time",
    x_breaks = NULL,
    x_labs = NULL,
    x_max = NULL,
    x_name = NULL,
    y_lab = "NAR"
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(dplyr::all_of(time_name)))
  }
  if (is.null(x_breaks)) {
    x_breaks <- round(seq(from = 0.0, to = x_max, length.out = 10))
  }
  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  
  # Data prep.
  if (is.null(specified_nars)) {
    data <- data %>%
      dplyr::rename(
        status = {{status_name}},
        time = {{time_name}}
      )
    df <- OneSampleNARs(data, x_breaks)
  } else {
    stopifnot(
      "specified_nars must contain time and nar." = 
        all(c("time", "nar") %in% colnames(specified_nars))
    )
    df <- specified_nars
  }
  
  # Placeholder for arm.
  df$arm <- 0
  df$arm <- factor(
    df$arm,
    levels = 0,
    labels = y_lab
  )
  
  # Plotting.
  arm <- NULL
  nar <- NULL
  time <- NULL
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = time, y = arm, label = nar)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) + 
    ggplot2::scale_y_discrete(
      name = NULL
    )
  return(q)
}


# -----------------------------------------------------------------------------
# Two sample.
# -----------------------------------------------------------------------------

#' Plot Two Sample Numbers at Risk
#' 
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param specified_nars Optional data.frame of specified NARs. Should include
#'   "arm", "time", and "nar". Overrides `data` if present.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis name.
#' @param y_labs Y-axis labels.
#' @return ggplot.
#' @importFrom dplyr "%>%"
#' @export
PlotTwoSampleNARs <- function(
  data,
  arm_name = "arm",
  specified_nars = NULL,
  status_name = "status",
  time_name = "time",
  x_breaks = NULL,
  x_labs = NULL,
  x_max = NULL,
  x_name = NULL,
  y_labs = c("Ctrl", "Trt")
) {
  
  # Defaults.
  if (is.null(x_max)) {
    x_max <- max(data %>% dplyr::select(dplyr::all_of(time_name)))
  }
  if (is.null(x_breaks)) {
    x_breaks <- round(seq(from = 0.0, to = x_max, length.out = 10))
  }
  if (is.null(x_labs)) {
    x_labs <- x_breaks
  }
  
  # Data prep.
  if (is.null(specified_nars)) {
    data <- data %>%
      dplyr::rename(
        arm = {{arm_name}},
        status = {{status_name}},
        time = {{time_name}}
      )
    df0 <- data %>%
      dplyr::filter(arm == 0) %>%
      OneSampleNARs(x_breaks) %>%
      dplyr::mutate(arm = 0)
    df1 <- data %>%
      dplyr::filter(arm == 1) %>%
      OneSampleNARs(x_breaks) %>%
      dplyr::mutate(arm = 1)
    df <- rbind(df0, df1)
  } else {
    df <- specified_nars
    stopifnot(
      "specified_nars must contain arm, time, and nar." = 
        all(c("arm", "time", "nar") %in% colnames(specified_nars))
    )
    stopifnot(
      "The levels of arm in specified_nars must be c(0, 1)." = 
        all(sort(unique(specified_nars$arm)) == c(0, 1))
    )
  }
  df$arm <- factor(
    df$arm,
    levels = c(0, 1),
    labels = y_labs,
    ordered = TRUE
  )
  
  # Plotting.
  arm <- NULL
  nar <- NULL
  time <- NULL
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = time, y = arm, label = nar)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) + 
    ggplot2::scale_y_discrete(
      name = NULL
    )
  return(q)
}
