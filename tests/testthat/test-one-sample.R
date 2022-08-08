library(testthat)

test_that("Test percentile calculation", {
  
  withr::local_seed(101)
  
  data <- GenData(base_event_rate = 1.0, n = 1e3)
  pct <- OneSamplePercentiles(data, p = c(0.75, 0.50, 0.25))
  expect_equal(pct$time[1], -log(0.75), tolerance = 0.05)
  expect_equal(pct$time[2], -log(0.50), tolerance = 0.05)
  expect_equal(pct$time[3], -log(0.25), tolerance = 0.05)

})


test_that("Test RMST calculation", {
  
  data <- data.frame(
    time = c(1, 1, 2, 2, 3, 3, 4, 4),
    status = rep(1, 8)
  )
  
  # Observed.
  obs <- OneSampleRMST(data, tau = 4)
  
  # Variance calculation.
  mu_t <- c(0.25 + 0.50 + 0.75, 0.25 + 0.50, 0.25, 0)
  dn_t <- c(2, 2, 2, 2)
  y_t <- c(8, 6, 4, 2)
  var <- sum(mu_t^2 * dn_t / y_t^2)
  
  # Expected.
  expect_equal(obs$auc, 2.5)
  expect_equal(obs$se, sqrt(var))
  
  
})
