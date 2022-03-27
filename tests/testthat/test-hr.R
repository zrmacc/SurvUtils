library(testthat)

test_that("Test hazard ratio.", {
  
  withr::local_seed(102)
  
 # Case: HR = 0.8
  arm1 <- GenData(base_event_rate = 0.8, n = 1e3)
  arm1$arm <- 1
  arm0 <- GenData(base_event_rate = 1.0, n = 1e3)
  arm0$arm <- 0
  data <- rbind(arm1, arm0)
  
  hr <- CalcHR(data)
  expect_equal(hr$hr, 0.8, tolerance = 0.01)
  
})