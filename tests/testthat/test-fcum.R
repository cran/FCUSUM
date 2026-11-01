make_data <- function(n = 150, p = 1, k = 1.5, amp = 2, coint = TRUE, seed = 1) {
  set.seed(seed)
  x <- apply(matrix(rnorm(n * p), n, p), 2, cumsum)
  x <- matrix(x, n, p, dimnames = list(NULL, paste0("x", seq_len(p))))
  u <- if (coint) rnorm(n) else cumsum(rnorm(n))
  brk <- amp * sin(2 * pi * k * seq_len(n) / n)
  y <- 1 + brk + drop(x %*% rep(1, p)) + u
  list(y = y, x = x)
}

test_that("fcum returns a well-formed object", {
  d <- make_data()
  r <- fcum(d$y, d$x[, 1], nboot = 49, seed = 1)
  expect_s3_class(r, "fcum")
  expect_true(is.finite(r$statistic) && r$statistic > 0)
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_named(r$critical_values, c("10%", "5%", "1%"))
  expect_true(all(diff(r$critical_values) > 0))   # right-tail ordering
  expect_equal(r$p, 1)
  expect_true(r$test %in% c("Fourier CUSUM", "Xiao-Phillips CUSUM"))
  expect_output(print(r), "sup-F")
  expect_output(summary(r), "DOLS coefficients")
})

test_that("critical values are upper-tail quantiles (bug in 1.0.0)", {
  d <- make_data()
  r <- fcum(d$y, d$x, cv = "table")
  expect_gt(unname(r$critical_values["5%"]), 0.5)
  r0 <- fcum(d$y, d$x, fourier = "never")
  expect_gt(unname(r0$critical_values["5%"]), 1)
})

test_that("a clear break is detected and handled by the Fourier stage", {
  d <- make_data(n = 300, k = 1.5, amp = 3)
  r <- fcum(d$y, d$x, kstar = 3, nboot = 49, seed = 2)
  expect_true(r$pretest$fourier_retained)
  expect_equal(r$test, "Fourier CUSUM")
  expect_lt(abs(r$best_frequency - 1.5), 0.1)
})

test_that("no break leads to the Xiao-Phillips stage most of the time", {
  kept <- vapply(1:20, function(i) {
    d <- make_data(n = 200, amp = 0, seed = 100 + i)
    fcum(d$y, d$x, cv = "table")$pretest$fourier_retained
  }, logical(1))
  expect_lt(mean(kept), 0.3)
})

test_that("multiple regressors are handled (bug in 1.0.0: flattened x)", {
  d <- make_data(p = 3)
  r <- fcum(d$y, d$x, kstar = 2, cv = "table")
  expect_equal(r$p, 3)
  r2 <- fcum(d$y, as.data.frame(d$x), kstar = 2, cv = "table")
  expect_equal(unname(r$statistic), unname(r2$statistic))
})

test_that("integer grid returns an integer frequency", {
  d <- make_data(k = 2)
  r <- fcum(d$y, d$x, kstar = 3, freq = "integer", cv = "table")
  expect_true(r$best_frequency %in% 1:3)
})

test_that("bootstrap is reproducible and both schemes run", {
  d <- make_data(n = 120, amp = 3)
  r1 <- fcum(d$y, d$x, nboot = 49, seed = 7, fourier = "always")
  r2 <- fcum(d$y, d$x, nboot = 49, seed = 7, fourier = "always")
  expect_identical(r1$critical_values, r2$critical_values)
  r3 <- fcum(d$y, d$x, nboot = 49, seed = 7, fourier = "always", boot = "sieve")
  expect_match(r3$cv_source, "sieve")
})

test_that("configurations outside the tables are simulated", {
  d <- make_data(n = 150, p = 5, amp = 0)
  expect_message(r <- fcum(d$y, d$x, kstar = 1, nboot = 49, seed = 1),
                 "simulating")
  expect_true(is.finite(r$statistic))
})

test_that("all bandwidth rules run", {
  d <- make_data()
  for (b in c("nw", "kurozumi", "andrews")) {
    expect_s3_class(fcum(d$y, d$x, bw = b, cv = "table"), "fcum")
  }
})

test_that("input validation", {
  d <- make_data()
  expect_error(fcum(d$y, d$x, kstar = 2.5), "positive integer")
  expect_error(fcum(d$y[-1], d$x), "same number")
  y_na <- d$y; y_na[3] <- NA
  expect_error(fcum(y_na, d$x), "missing")
  expect_error(fcum(d$y, d$x, q = -1), "non-negative")
  expect_error(fcum(d$y, d$x, nboot = 5), "at least 19")
})

test_that("fcum_cv_sim returns ordered quantiles", {
  s <- fcum_cv_sim(n = 100, p = 1, type = "plain", nsim = 200, seed = 1)
  expect_true(all(diff(s$quantiles) >= 0))
})

test_that("empirical size is close to nominal (no break, table route)", {
  skip_on_cran()
  rej <- vapply(1:200, function(i) {
    d <- make_data(n = 200, amp = 0, seed = 1000 + i)
    r <- fcum(d$y, d$x, kstar = 2, nboot = 49, seed = i)
    unname(r$statistic > r$critical_values["5%"])
  }, logical(1))
  expect_lt(mean(rej), 0.12)
})

test_that("power against no cointegration", {
  skip_on_cran()
  rej <- vapply(1:40, function(i) {
    d <- make_data(n = 300, coint = FALSE, seed = 2000 + i)
    r <- fcum(d$y, d$x, kstar = 2, nboot = 99, seed = i)
    unname(r$statistic > r$critical_values["5%"])
  }, logical(1))
  expect_gt(mean(rej), 0.35)
})
