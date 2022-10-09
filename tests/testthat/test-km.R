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


test_that("Influence function calculation.", {
  
  withr::local_seed(101)
  n <- 1000
  data <- GenData(n = n)
  km <- GetCurves(data)
  
  tau <- 1.0
  ref_var <- km@SurvVar(tau)
  influence <- KMInfluence(data, tau = tau)
  inf_var <- sum(influence^2) / (n^2)
  expect_equal(n * ref_var, n * inf_var, tolerance = 0.05)
  
})
