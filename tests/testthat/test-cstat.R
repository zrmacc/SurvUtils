library(testthat)

test_that("Test simple c-stat calcualtion.", {
  
  data <- data.frame(
    time = c(0, 1, 2, 3),
    status = c(0, 1, 0, 1),
    risk = c(4, 3, 2, 1)
  )
  
  # With pseudo counts.
  obs <- Cstat(data, pseudo_counts = TRUE)
  exp <- (0.5 + 2) / (1.0 + 2)
  expect_equal(obs, exp)
  
  # Without pseudo counts.
  obs <- Cstat(data, pseudo_counts = FALSE)
  exp <- 2 / 2
  expect_equal(obs, exp)
  
})


test_that("Test weighted C-stat calculation.", {
  
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
  obs <- WeightedCstat(
    train_data = train_data, 
    test_data = test_data, 
    tau = 4
  )
  expect_equal(obs, exp)
  
  # Risk ordering is partially correct. 
  test_data <- data.frame(
    time = c(0, 1, 2),
    status = c(1, 1, 1),
    risk = c(3, 1, 2)
  )
  exp <- (2 + 0 + 0) / (2 + 1 + 0)
  obs <- WeightedCstat(
    train_data = train_data,
    test_data = test_data, 
    tau = 4
  )
  expect_equal(obs, exp)
  
  # Risk ordering is anti-correct.
  test_data <- data.frame(
    time = c(0, 1, 2),
    status = c(1, 1, 1),
    risk = c(1, 2, 3)
  )
  exp <- 0
  obs <- WeightedCstat(
    train_data = train_data, 
    test_data = test_data, 
    tau = 4
  )
  expect_equal(obs, exp)
  
})


test_that("Test simple c-stat calcualtion.", {
  
  withr::local_seed(101)
  
  # Test Cox model comparison.
  n <- 1000
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  data <- SurvUtils::GenData(
    covariates = cbind(x1, x2),
    beta_event = c(1.0, -1.0)
  )
  
  eval <- CompreCoxCstat(
    status = data$status,
    time = data$time,
    x1 = data %>% dplyr::select(x1, x2),
    x2 = data %>% dplyr::select(x1)
  )
  
  obs_diff <- mean(eval$diff)
  obs_ratio <- mean(eval$ratio)
  expect_gt(obs_diff, 0.0)
  expect_gt(obs_ratio, 1.0)
  
})
