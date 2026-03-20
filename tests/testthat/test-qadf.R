test_that("qadf returns correct class and structure", {
  set.seed(1)
  y <- cumsum(rnorm(60))
  res <- qadf(y, tau = 0.5, model = "c", max_lags = 4, ic = "aic")
  expect_s3_class(res, "qadf")
  expect_named(res, c("statistic", "coef_stat", "rho_tau", "rho_ols",
                       "alpha_tau", "delta2", "half_life", "opt_lags",
                       "nobs", "critical_values", "tau", "model", "ic",
                       "varname"))
})

test_that("qadf validates tau", {
  expect_error(qadf(rnorm(50), tau = 0), "strictly between 0 and 1")
  expect_error(qadf(rnorm(50), tau = 1), "strictly between 0 and 1")
  expect_error(qadf(rnorm(50), tau = -0.1), "strictly between 0 and 1")
})

test_that("qadf validates model", {
  expect_error(qadf(rnorm(50), model = "none"), "\"c\" \\(constant\\)")
})

test_that("qadf validates ic", {
  expect_error(qadf(rnorm(50), ic = "hqic"), "\"aic\", \"bic\", or \"tstat\"")
})

test_that("qadf rejects short series", {
  expect_error(qadf(rnorm(15)), "at least 20")
})

test_that("qadf constant model produces negative t-stat for unit root", {
  set.seed(42)
  y <- cumsum(rnorm(80))
  res <- qadf(y, tau = 0.5, model = "c", max_lags = 4)
  expect_true(is.finite(res$statistic))
  expect_equal(res$tau, 0.5)
  expect_equal(res$model, "c")
})

test_that("qadf trend model runs without error", {
  set.seed(7)
  y <- cumsum(rnorm(80))
  res <- qadf(y, tau = 0.5, model = "ct", max_lags = 3, ic = "bic")
  expect_s3_class(res, "qadf")
})

test_that("critical values are negative and ordered", {
  set.seed(3)
  y <- cumsum(rnorm(60))
  res <- qadf(y, tau = 0.5, model = "c")
  cv  <- res$critical_values
  expect_true(cv["cv1"] < cv["cv5"])
  expect_true(cv["cv5"] < cv["cv10"])
  expect_true(all(cv < 0))
})

test_that("print.qadf is silent (uses message, not cat/print)", {
  set.seed(1)
  y <- cumsum(rnorm(60))
  res <- qadf(y, tau = 0.5)
  # print.qadf uses message(), not cat() — output goes to stderr
  expect_message(print(res))
  # No output to stdout
  out <- capture.output(print(res))
  expect_equal(length(out), 0L)
})
