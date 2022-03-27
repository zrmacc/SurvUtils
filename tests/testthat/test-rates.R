library(testthat)

test_that("Test rate comparison.", {
  
  withr::local_seed(104)
  
  arm1 <- GenData(base_event_rate = 0.5, n = 1e3)
  arm1$arm <- 1
  arm0 <- GenData(base_event_rate = 1.0, n = 1e3)
  arm0$arm <- 0
  data <- rbind(arm1, arm0)
  
  # Expected values:
  # rate difference = 0.50 - 0.25 = 0.25
  # rate ratio = 0.50 / 0.25 = 2.0
  # odds ratio = 0.50 / (1 - 0.50) * (1 - 0.25) / 0.25 = 3.0
  compare_rates <- CompareRates(data, tau = log(4.0))
  p0 <- compare_rates@Marginal$rate[1]
  p1 <- compare_rates@Marginal$rate[2]
  
  rd <- compare_rates@Contrasts$est[1]
  rr <- compare_rates@Contrasts$est[2]
  or <- compare_rates@Contrasts$est[3]
  
  expect_equal(p0, 0.25, tolerance = 0.05)
  expect_equal(p1, 0.50, tolerance = 0.05)
  
  expect_equal(rd, 0.25, tolerance = 0.05)
  expect_equal(rr, 2.0, tolerance = 0.2)
  expect_equal(or, 3.0, tolerance = 0.2)

})