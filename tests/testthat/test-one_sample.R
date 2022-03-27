library(testthat)

test_that("Test percentile calculation", {
  
  withr::local_seed(101)
  
  data <- GenData(base_event_rate = 1.0, n = 1e3)
  pct <- OneSamplePercentiles(data, p = c(0.75, 0.50, 0.25))
  expect_equal(pct$time[1], -log(0.75), tolerance = 0.05)
  expect_equal(pct$time[2], -log(0.50), tolerance = 0.05)
  expect_equal(pct$time[3], -log(0.25), tolerance = 0.05)

})