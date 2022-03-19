library(testthat)

test_that("Test Kaplan-Meier estimator.", {
  
  # Case without ties.
  data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 1, 0)
  )
  
  exp_haz <- c(0, 1 / 4, 0 / 3, 1 / 2, 0 / 1)
  km <- TabulateKM(data)
  expect_equal(km$haz, exp_haz)
  expect_equal(km$surv, cumprod(1 - exp_haz))
  
  # Case with ties.
  data <- data.frame(
    time = c(1, 1, 2, 2, 3),
    status = c(1, 1, 0, 0, 1)
  )
  
  exp_haz <- c(0, 2 / 5, 0, 1/1)
  km <- TabulateKM(data)
  expect_equal(km$haz, exp_haz)
  expect_equal(km$surv, cumprod(1 - exp_haz))
  
})