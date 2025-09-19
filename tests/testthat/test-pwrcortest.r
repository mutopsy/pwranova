test_that("pwrcortest: returns expected columns and classes (t-method)", {
  out <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = 80, alpha = 0.05, power = NULL, rho = 0.30
  )
  expect_s3_class(out, "cal_power")
  expect_true(is.data.frame(out))
  expect_true(all(c("df","n_total","alpha","power","rho","t_critical","ncp") %in% names(out)))
  expect_equal(out$df, 78)
  expect_equal(out$n_total, 80)
  expect_gt(out$power, 0)
  expect_lt(out$power, 1)
})

test_that("pwrcortest: power(N) -> solve N -> re-evaluate gives target power (t-method)", {
  # Step 1: compute power at a chosen N
  pwr <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = 80, alpha = 0.05, rho = 0.30
  )
  target <- as.numeric(pwr$power)

  # Step 2: solve N for that power
  solN <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = NULL, alpha = 0.05, power = target, rho = 0.30
  )
  expect_s3_class(solN, "cal_n")
  expect_true(solN$n_total >= 3)

  # Step 3: plug back and check power
  pwr2 <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = solN$n_total, alpha = 0.05, rho = 0.30
  )
  expect_equal(as.numeric(pwr2$power), target, tolerance = 1e-8)
})

test_that("pwrcortest: solve alpha then reproduce target power (t-method)", {
  target_power <- 0.85
  solA <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = 120, alpha = NULL, power = target_power, rho = 0.25
  )
  expect_s3_class(solA, "cal_alpha")
  expect_true(solA$alpha > 0 && solA$alpha < 1)

  # Recompute power using solved alpha
  pwr <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = 120, alpha = solA$alpha, rho = 0.25
  )
  expect_equal(as.numeric(pwr$power), target_power, tolerance = 1e-10)
})

test_that("pwrcortest: solve rho then reproduce target power (t-method)", {
  target_power <- 0.80
  solR <- pwrcortest(
    alternative = "one.sided", method = "t",
    n_total = 90, alpha = 0.05, power = target_power, rho = NULL
  )
  expect_s3_class(solR, "cal_es")
  expect_true(abs(solR$rho) > 0 && abs(solR$rho) < 1)

  # Recompute power using solved rho
  pwr <- pwrcortest(
    alternative = "one.sided", method = "t",
    n_total = 90, alpha = 0.05, rho = solR$rho
  )
  expect_equal(as.numeric(pwr$power), target_power, tolerance = 1e-10)
})

test_that("pwrcortest: sign of rho is irrelevant for two-sided", {
  p_pos <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = 60, alpha = 0.05, rho = 0.3
  )$power
  p_neg <- pwrcortest(
    alternative = "two.sided", method = "t",
    n_total = 60, alpha = 0.05, rho = -0.3
  )$power
  expect_equal(as.numeric(p_pos), as.numeric(p_neg), tolerance = 0)
})

test_that("pwrcortest: z-method columns and roundtrip (power -> N -> power)", {
  out <- pwrcortest(
    alternative = "two.sided", method = "z",
    n_total = 80, alpha = 0.05, rho = 0.20
  )
  expect_s3_class(out, "cal_power")
  expect_true(all(c("n_total","alpha","power","rho","z_critical","ncp") %in% names(out)))
  expect_false("df" %in% names(out))
  target <- as.numeric(out$power)

  solN <- pwrcortest(
    alternative = "two.sided", method = "z",
    n_total = NULL, alpha = 0.05, power = target, rho = 0.20
  )
  expect_s3_class(solN, "cal_n")
  pwr2 <- pwrcortest(
    alternative = "two.sided", method = "z",
    n_total = solN$n_total, alpha = 0.05, rho = 0.20
  )
  expect_equal(as.numeric(pwr2$power), target, tolerance = 1e-8)
})

test_that("pwrcortest: z-method with bias correction behaves consistently", {
  # Solve alpha then reproduce power (bias_correction = TRUE)
  target_power <- 0.9
  solA <- pwrcortest(
    alternative = "two.sided", method = "z", bias_correction = TRUE,
    n_total = 70, alpha = NULL, power = target_power, rho = 0.25
  )
  expect_s3_class(solA, "cal_alpha")
  expect_true(solA$alpha > 0 && solA$alpha < 1)

  pwr <- pwrcortest(
    alternative = "two.sided", method = "z", bias_correction = TRUE,
    n_total = 70, alpha = solA$alpha, rho = 0.25
  )
  expect_equal(as.numeric(pwr$power), target_power, tolerance = 1e-10)
})

test_that("pwrcortest: minimal detectable rho (solve es) reproduces power (z-method)", {
  target_power <- 0.80
  solR <- pwrcortest(
    alternative = "two.sided", method = "z", bias_correction = FALSE,
    n_total = 100, alpha = 0.05, power = target_power, rho = NULL
  )
  expect_s3_class(solR, "cal_es")
  expect_true(abs(solR$rho) > 0 && abs(solR$rho) < 1)

  pwr <- pwrcortest(
    alternative = "two.sided", method = "z",
    n_total = 100, alpha = 0.05, rho = solR$rho
  )
  expect_equal(as.numeric(pwr$power), target_power, tolerance = 1e-10)
})

test_that("pwrcortest: method='t' ignores bias_correction but stays valid", {
  expect_warning(
    out <- pwrcortest(
      alternative = "one.sided", method = "t", bias_correction = TRUE,
      n_total = 40, alpha = 0.05, rho = 0.35
    ),
    "ignored"
  )
  expect_s3_class(out, "cal_power")
  expect_true(all(c("df","t_critical","ncp") %in% names(out)))
})
