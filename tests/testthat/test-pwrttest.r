test_that("pwrttest computes power (cal_power) for two-sample, two.sided", {
  # Two-sample, equal allocation: df = N - 2, ncp = d * sqrt(N) / 2
  N <- 128L
  d <- 0.50
  alpha <- 0.05

  out <- pwrttest(
    paired = FALSE, onesample = FALSE, alternative = "two.sided",
    n_total = N, delta = d, alpha = alpha, power = NULL
  )

  expect_s3_class(out, "cal_power")
  expect_equal(out$df, N - 2)
  expect_equal(out$n_total, N)
  # critical value under central t
  expect_equal(out$t_critical, qt(1 - alpha/2, df = N - 2), tolerance = 1e-12)
  # ncp formula
  expect_equal(out$ncp, d * sqrt(N) / 2, tolerance = 1e-12)

  # recompute power independently (two-sided)
  crit <- out$t_critical
  ncp  <- out$ncp
  pow_expected <- (1 - pt(crit, df = out$df, ncp = ncp)) + pt(-crit, df = out$df, ncp = ncp)
  expect_equal(out$power, pow_expected, tolerance = 1e-10)
})

test_that("pwrttest solves N (cal_n) for two-sample using cohensf (d = 2f)", {
  # target: two-sample, two.sided, specify f and solve N
  f <- 0.25
  alpha <- 0.05
  target_power <- 0.80

  out <- pwrttest(
    paired = FALSE, onesample = FALSE, alternative = "two.sided",
    n_total = NULL, cohensf = f, alpha = alpha, power = target_power
  )

  expect_s3_class(out, "cal_n")
  expect_true(out$n_total %% 2 == 0)           # equal allocation enforced
  expect_equal(out$df, out$n_total - 2)
  # d = 2f により内部 delta=0.5 相当
  expect_equal(out$delta, 2 * f, tolerance = 1e-12)
  # achieved power should be >= target (first N that reaches target)
  expect_true(out$power >= target_power - 1e-12)
})

test_that("pwrttest paired design accepts f (dz = f) and computes power", {
  # paired: df = N - 1, ncp = dz * sqrt(N)
  N <- 40L
  f <- 0.30  # here f should be equal to dz
  alpha <- 0.05

  out <- pwrttest(
    paired = TRUE, onesample = FALSE, alternative = "one.sided",
    n_total = N, cohensf = f, alpha = alpha, power = NULL
  )

  expect_s3_class(out, "cal_power")
  expect_equal(out$df, N - 1)
  expect_equal(out$ncp, f * sqrt(N), tolerance = 1e-12)
  expect_equal(out$t_critical, qt(1 - alpha, df = N - 1), tolerance = 1e-12)

  # one-sided power check
  crit <- out$t_critical
  pow_expected <- 1 - pt(crit, df = out$df, ncp = out$ncp)
  expect_equal(out$power, pow_expected, tolerance = 1e-10)
})

test_that("pwrttest one-sample: solve alpha (cal_alpha) given power", {
  N <- 30L
  d <- 0.40
  target_power <- 0.90

  out <- pwrttest(
    onesample = TRUE, paired = FALSE, alternative = "two.sided",
    n_total = N, delta = d, alpha = NULL, power = target_power
  )

  expect_s3_class(out, "cal_alpha")
  expect_equal(out$df, N - 1)

  # back-check: with found alpha -> recompute power equals target (within tol)
  crit <- out$t_critical
  pow_expected <- (1 - pt(crit, df = out$df, ncp = d * sqrt(N))) +
    pt(-crit, df = out$df, ncp = d * sqrt(N))
  expect_equal(out$power, pow_expected, tolerance = 1e-8)
  # alpha from central t two-sided
  expect_equal(out$alpha, 2 * (1 - pt(crit, df = out$df)), tolerance = 1e-12)
})

test_that("pwrttest solves minimal detectable effect (cal_es), two-sample two.sided", {
  N <- 100L
  alpha <- 0.05
  target_power <- 0.80

  out <- pwrttest(
    paired = FALSE, onesample = FALSE, alternative = "two.sided",
    n_total = N, delta = NULL, cohensf = NULL, peta2 = NULL,
    alpha = alpha, power = target_power
  )

  expect_s3_class(out, "cal_es")
  # recompute power using solved effect
  d_solved <- out$delta
  df <- out$df
  crit <- qt(1 - alpha/2, df = df)
  ncp <- d_solved * sqrt(N) / 2
  pow_expected <- (1 - pt(crit, df = df, ncp = ncp)) + pt(-crit, df = df, ncp = ncp)
  expect_equal(out$t_critical, crit, tolerance = 1e-12)
  expect_equal(out$ncp, ncp, tolerance = 1e-10)
  expect_equal(out$power, pow_expected, tolerance = 1e-8)
  # f と peta2 の一貫性（two-sample: d = 2f）
  expect_equal(out$cohensf, out$delta / 2, tolerance = 1e-12)
  expect_true(out$peta2 > 0 && out$peta2 < 1)
})

test_that("pwrttest input validation: n_total must be scalar integer", {
  expect_error(pwrttest(onesample = TRUE, n_total = c(20, 30), delta = 0.4, alpha = 0.05, power = NULL))
  expect_error(pwrttest(onesample = TRUE, n_total = 20.5,     delta = 0.4, alpha = 0.05, power = NULL))
})

test_that("pwrttest enforces exactly one NULL among (n_total, delta/f/eta2, alpha, power)", {
  # Here two NULLs (n_total and power) -> error
  expect_error(pwrttest(onesample = TRUE, n_total = NULL, delta = 0.4, alpha = 0.05, power = NULL))
})

test_that("pwrttest ignores 'paired' when onesample = TRUE (with warning)", {
  expect_warning(
    pwrttest(onesample = TRUE, paired = TRUE, n_total = 16, delta = 0.5, alpha = 0.05, power = NULL)
  )
})

test_that("pwrttest two-sample requires even n_total", {
  expect_error(
    pwrttest(paired = FALSE, onesample = FALSE, n_total = 21, delta = 0.5, alpha = 0.05, power = NULL)
  )
})

test_that("pwrttest: effect-size inputs precedence and design-dependent conversions", {
  # precedence: delta > f > peta2

  expect_warning(
    out1 <- pwrttest(paired = FALSE, onesample = FALSE,
                     n_total = 80, delta = 0.6, cohensf = 0.25, peta2 = 0.06,
                     alpha = 0.05, power = NULL),
    "ignored"
  )
  expect_equal(out1$delta, 0.6, tolerance = 1e-12)

  # paired: f allowed (dz = f)
  out2 <- pwrttest(paired = TRUE, onesample = FALSE,
                   n_total = 30, cohensf = 0.35, alpha = 0.05, power = NULL)
  expect_equal(out2$delta, 0.35, tolerance = 1e-12)

  # one-sample: f / peta2 not supported -> should error
  expect_error(
    pwrttest(onesample = TRUE, n_total = 20, cohensf = 0.25, alpha = 0.05, power = NULL)
  )
  expect_error(
    pwrttest(onesample = TRUE, n_total = 20, peta2 = 0.06, alpha = 0.05, power = NULL)
  )
})
