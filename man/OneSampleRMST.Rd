% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/one_sample.R
\name{OneSampleRMST}
\alias{OneSampleRMST}
\title{One Sample RMST}
\usage{
OneSampleRMST(
  data,
  tau = NULL,
  alpha = 0.05,
  status_name = "status",
  time_name = "time"
)
}
\arguments{
\item{data}{Data.frame.}

\item{tau}{Truncation time.}

\item{alpha}{Type I error.}

\item{status_name}{Name of status column.}

\item{time_name}{Name of time column.}
}
\value{
Data.frame.
}
\description{
One Sample RMST
}
