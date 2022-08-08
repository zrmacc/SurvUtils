test_that("Test two sample RMST comparison.", {
  
  withr::local_seed(101)
  
  # True difference.
  arm1 <- GenData(base_event_rate = 0.5, n = 1e3)
  arm1$arm <- 1
  arm0 <- GenData(base_event_rate = 1.0, n = 1e3)
  arm0$arm <- 0
  data <- rbind(arm1, arm0)
  
  # Observed.
  obs <- CompareRMSTs(data, tau = 4)
  marg <- obs@Marginal
  cont <- obs@Contrasts

  # Expected value.
  auc0 <- 1 - exp(-4)
  auc1 <- 2 * (1 - exp(-2))
  expect_equal(marg$auc[1], auc0, tolerance = 0.05)
  expect_equal(marg$auc[2], auc1, tolerance = 0.05)
  expect_equal(cont$est[1], auc1 - auc0, tolerance = 0.05)
  expect_equal(cont$est[2], auc1 / auc0, tolerance = 0.05)
  expect_equal(cont$p[1], 0)
  expect_equal(cont$p[2], 0)
  
  
  # No true difference.
  arm1 <- GenData(base_event_rate = 1.0, n = 1e3)
  arm1$arm <- 1
  arm0 <- GenData(base_event_rate = 1.0, n = 1e3)
  arm0$arm <- 0
  data <- rbind(arm1, arm0)
  
  # Observed.
  obs <- CompareRMSTs(data, tau = 4)
  marg <- obs@Marginal
  cont <- obs@Contrasts
  
  # Expected value.
  auc0 <- 1 - exp(-4)
  expect_equal(marg$auc[1], auc0, tolerance = 0.05)
  expect_equal(marg$auc[2], auc0, tolerance = 0.05)
  expect_equal(cont$est[1], 0, tolerance = 0.05)
  expect_equal(cont$est[2], 1, tolerance = 0.05)
  expect_gt(cont$p[1], 0.05)
  expect_gt(cont$p[2], 0.05)
  
  
})