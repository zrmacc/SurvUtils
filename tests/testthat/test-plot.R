library(testthat)

test_that("Plot functions return ggplot objects.", {
  withr::local_seed(111)
  data <- GenData(n = 50)
  data_two <- rbind(
    cbind(GenData(n = 25), arm = 0),
    cbind(GenData(n = 25), arm = 1)
  )

  p1 <- PlotOneSampleKM(data)
  expect_s3_class(p1, "ggplot")

  p2 <- PlotOneSampleNARs(data)
  expect_s3_class(p2, "ggplot")

  p3 <- PlotOneSampleAUC(data, tau = 2)
  expect_s3_class(p3, "ggplot")

  p4 <- PlotOneSampleCH(data)
  expect_s3_class(p4, "ggplot")

  p5 <- PlotTwoSampleKM(data_two)
  expect_s3_class(p5, "ggplot")

  p6 <- PlotTwoSampleNARs(data_two)
  expect_s3_class(p6, "ggplot")

  p7 <- PlotTwoSampleCH(data_two, tau = 2)
  expect_s3_class(p7, "ggplot")
})
