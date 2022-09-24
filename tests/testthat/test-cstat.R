library(testthat)

test_that("Test C-stat calculation.", {
  
  # No censoring, all weights are 1.
  train_data <- data.frame(
    time = c(0, 1, 2),
    status = c(1, 1, 1)
  )
  
  # Risk is properly ordered, expecting C = 1.
  test_data <- data.frame(
    time = c(0, 1, 2),
    status = c(1, 1, 1),
    risk = c(3, 2, 1)
  )
  
  exp <- 1
  obs <- Cstat(train_data, test_data, tau = 4)
  expect_equal(obs, exp)
  
  # Risk ordering is partially correct. 
  test_data <- data.frame(
    time = c(0, 1, 2),
    status = c(1, 1, 1),
    risk = c(3, 1, 2)
  )
  exp <- (2 + 0 + 0) / (2 + 1 + 0)
  obs <- Cstat(train_data, test_data, tau = 4)
  expect_equal(obs, exp)
  
  # Risk ordering is anti-correct.
  test_data <- data.frame(
    time = c(0, 1, 2),
    status = c(1, 1, 1),
    risk = c(1, 2, 3)
  )
  exp <- 0
  obs <- Cstat(train_data, test_data, tau = 4)
  expect_equal(obs, exp)
  
  
})