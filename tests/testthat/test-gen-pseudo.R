library(testthat)

test_that("GenPseudo type = 'prob' returns param + influence = pseudo.", {
  withr::local_seed(107)
  data <- GenData(n = 200, tau = 4)
  out <- GenPseudo(data, tau = 1.0, type = "prob")
  expect_true("influence" %in% colnames(out))
  expect_true("pseudo" %in% colnames(out))
  param <- OneSampleRates(data, tau = 1.0)$rate
  expect_equal(out$pseudo, param + out$influence)
  expect_equal(nrow(out), nrow(data))
})

test_that("GenPseudo type = 'rmst' returns param + influence = pseudo.", {
  withr::local_seed(108)
  data <- GenData(n = 200, tau = 4)
  out <- GenPseudo(data, tau = 2.0, type = "rmst")
  expect_true("influence" %in% colnames(out))
  expect_true("pseudo" %in% colnames(out))
  param <- OneSampleRMST(data, tau = 2.0)$auc
  expect_equal(out$pseudo, param + out$influence)
  expect_equal(nrow(out), nrow(data))
})

test_that("GenPseudo type = 'cic' returns param + influence = pseudo.", {
  withr::local_seed(109)
  data <- GenCRData(n = 200, tau = 4)
  out <- GenPseudo(data, tau = 1.0, type = "cic")
  expect_true("influence" %in% colnames(out))
  expect_true("pseudo" %in% colnames(out))
  param <- OneSampleCIC(data, tau = 1.0)$rate
  expect_equal(out$pseudo, param + out$influence)
  expect_equal(nrow(out), nrow(data))
})

test_that("GenPseudo errors on invalid type.", {
  data <- data.frame(time = c(1, 2), status = c(1, 0))
  expect_error(GenPseudo(data, type = "invalid"), "type must be selected from among")
})
