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