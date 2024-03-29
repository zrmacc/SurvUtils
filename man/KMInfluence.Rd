% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kaplan_meier.R
\name{KMInfluence}
\alias{KMInfluence}
\title{Kaplan-Meier Influence}
\usage{
KMInfluence(data, tau, status_name = "status", time_name = "time")
}
\arguments{
\item{data}{Data.frame.}

\item{tau}{Truncation time at which to calculate the influence.}

\item{status_name}{Name of status column.}

\item{time_name}{Name of time column.}
}
\value{
Numeric vector.
}
\description{
Calculate Kaplan-Meier influence function for each observation.
}
