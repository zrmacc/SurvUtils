% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_nars.R
\name{PlotTwoSampleNARs}
\alias{PlotTwoSampleNARs}
\title{Plot Two Sample Numbers at Risk}
\usage{
PlotTwoSampleNARs(
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
)
}
\arguments{
\item{data}{Data.frame.}

\item{arm_name}{Name of arm column.}

\item{specified_nars}{Optional data.frame of specified NARs. Should include
"arm", "time", and "nar". Overrides \code{data} if present.}

\item{status_name}{Name of status column.}

\item{time_name}{Name of time column.}

\item{x_breaks}{X-axis breaks.}

\item{x_labs}{X-axis labels.}

\item{x_max}{X-axis upper limit.}

\item{x_name}{X-axis name.}

\item{y_labs}{Y-axis labels.}
}
\value{
ggplot.
}
\description{
Plot Two Sample Numbers at Risk
}
