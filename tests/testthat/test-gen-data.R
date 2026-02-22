library(testthat)

test_that("GenData returns correct structure and columns.", {
  withr::local_seed(101)

  # Default simple = TRUE: idx, time, status only
  data <- GenData(n = 50, simple = TRUE)
  expect_true(is.data.frame(data))
  expect_equal(colnames(data), c("idx", "time", "status"))
  expect_equal(nrow(data), 50)
  expect_equal(data$idx, 1:50)
  expect_true(all(data$status %in% c(0, 1)))
  expect_true(all(data$time >= 0))

  # simple = FALSE: includes additional columns
  data_full <- GenData(n = 50, simple = FALSE)
  expect_true(all(c("idx", "time", "status", "true_event_rate", "event_time", "censor_time") %in% colnames(data_full)))
  expect_equal(nrow(data_full), 50)
  # Status 1 iff event occurred before censoring
  expect_equal(data_full$status, 1 * (data_full$time == data_full$event_time))
})

test_that("GenData respects n and tau.", {
  withr::local_seed(102)
  data <- GenData(n = 100, tau = 2.0, simple = TRUE)
  expect_equal(nrow(data), 100)
  expect_true(all(data$time <= 2.0))
})

test_that("GenData with covariates uses nrow(covariates).", {
  withr::local_seed(103)
  x <- matrix(rnorm(40), ncol = 2)
  data <- GenData(covariates = x, base_event_rate = 1.0, simple = TRUE)
  expect_equal(nrow(data), 20)
})

test_that("GenCRData returns correct structure and status coding.", {
  withr::local_seed(104)

  data <- GenCRData(n = 50, simple = TRUE)
  expect_true(is.data.frame(data))
  expect_equal(colnames(data), c("idx", "time", "status"))
  expect_true(all(data$status %in% c(0, 1, 2)))
  expect_true(all(data$time >= 0))

  data_full <- GenCRData(n = 50, simple = FALSE)
  expect_true(all(c("idx", "time", "status", "event_time", "event_type") %in% colnames(data_full)))
  expect_equal(nrow(data_full), 50)
})

test_that("GenCRData respects tau.", {
  withr::local_seed(105)
  data <- GenCRData(n = 100, tau = 1.5, simple = TRUE)
  expect_true(all(data$time <= 1.5))
})
