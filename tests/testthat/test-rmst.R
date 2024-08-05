library(testthat)

test_that("Test restricted mean survival time calculation.", {

  # Test RMST calculation.
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 1, 1, 1, 0)
  
  # Truncation time of 5 (implicitly).
  obs <- RMST(status = status, time = time)
  exp <- 1.0 + 0.8 + 0.6 + 0.4 + 0.2
  expect_equal(obs, exp)
  
  # Truncation time of 2.
  obs <- RMST(status = status, time = time, tau = 2)
  exp <- 1.0 + 0.8
  expect_equal(obs, exp)
  
  # Truncation time of 7, do not extend.
  obs <- RMST(status = status, time = time, tau = 7)
  exp <- 1.0 + 0.8 + 0.6 + 0.4 + 0.2
  expect_equal(obs, exp)

  # Truncation time of 7, extend past last observation time.
  obs <- RMST(status = status, time = time, tau = 7, extend = TRUE)
  exp <- 1.0 + 0.8 + 0.6 + 0.4 + 0.2 + (7 - 5) * 0.2
  expect_equal(obs, exp)
  
})


test_that("Test influence function calculation.", {
  
  # Data.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 0)
  )
  n <- nrow(data)
  tau <- 4
  
  # Manual calculation.
  tab <- TabulateKM(data) %>%
    dplyr::filter(time <= tau)  
  
  event_times <- tab$time
  surv <- tab$surv
  nar <- tab$nar
  
  delta_t <- diff(c(event_times, tau))
  mu_t <- rev(cumsum(rev(delta_t * surv)))
  
  mart_t <- rbind(
    c(0, 1 - 0.2, 0, 0, 0),
    c(0, -0.2, 1 - 0.25, 0, 0),
    c(0, -0.2, -0.25, 1 - 1/3, 0),
    c(0, -0.2, -0.25, -1/3, 1 - 0.5),
    c(0, -0.2, -0.25, -1/3, -0.5)
  )
  
  exp_psi <- rep(0, n)
  for (i in 1:n) {
    exp_psi[i] <- -1 * sum(mu_t / (nar / n) * mart_t[i, ])
  }
  
  # Observed.
  obs <- RMSTInfluence(data, tau = tau)
  expect_equal(obs$influence, exp_psi)
  
})

