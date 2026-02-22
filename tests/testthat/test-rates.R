library(testthat)

test_that("RateDiff returns correct structure and calculation.", {
  rates <- data.frame(arm = c(0, 1), rate = c(0.2, 0.5), se = c(0.04, 0.05))
  rd <- RateDiff(rates)
  expect_equal(rd$stat, "rd")
  expect_equal(rd$est, 0.3)
  expect_equal(rd$se, sqrt(0.04^2 + 0.05^2))
  expect_true(rd$lower < rd$est)
  expect_true(rd$upper > rd$est)
  expect_true(rd$p >= 0 && rd$p <= 1)
})

test_that("RateRatio returns correct structure and calculation.", {
  rates <- data.frame(arm = c(0, 1), rate = c(0.25, 0.50), se = c(0.03, 0.04))
  rr <- RateRatio(rates)
  expect_equal(rr$stat, "rr")
  expect_equal(rr$est, 2.0)
  expect_true(rr$lower < rr$est)
  expect_true(rr$upper > rr$est)
  expect_true(rr$p >= 0 && rr$p <= 1)
})

test_that("OddsRatio returns correct structure and calculation.", {
  # odds0 = 0.25/0.75 = 1/3, odds1 = 0.5/0.5 = 1, OR = 3
  rates <- data.frame(arm = c(0, 1), rate = c(0.25, 0.50), se = c(0.03, 0.04))
  or <- OddsRatio(rates)
  expect_equal(or$stat, "or")
  expect_equal(or$est, (0.5 / 0.5) / (0.25 / 0.75), tolerance = 0.01)
  expect_true(or$lower < or$est)
  expect_true(or$upper > or$est)
  expect_true(or$p >= 0 && or$p <= 1)
})

test_that("CompareRates with return_surv = FALSE compares event rates.", {
  withr::local_seed(110)
  arm0 <- GenData(base_event_rate = 1.0, n = 500)
  arm0$arm <- 0
  arm1 <- GenData(base_event_rate = 0.5, n = 500)
  arm1$arm <- 1
  data <- rbind(arm0, arm1)
  comp <- CompareRates(data, tau = 1.0, return_surv = FALSE)
  # Event rate = 1 - survival; arm1 has lower event rate
  expect_lt(comp@Marginal$rate[2], comp@Marginal$rate[1])
  expect_equal(comp@Contrasts$stat, c("rd", "rr", "or"))
})

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