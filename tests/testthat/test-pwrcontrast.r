test_that("pwrcontrast calculates power (cal_power) for between-subjects correctly", {
  w <- c(1, -1)
  n_total <- 40
  f <- 0.25
  alpha <- 0.05

  out <- pwrcontrast(weight = w, paired = FALSE,
                     n_total = n_total, cohensf = f, alpha = alpha, power = NULL)

  expect_s3_class(out, "cal_power")
  expect_equal(out$df_num, 1)
  expect_equal(out$df_denom, n_total - length(w))
  expect_equal(out$ncp, f^2 * n_total, tolerance = 1e-12)

  # Check critical F and power numerically
  expected_crit <- qf(1 - alpha, df1 = 1, df2 = n_total - length(w))
  expected_pow  <- 1 - pf(expected_crit, df1 = 1, df2 = n_total - length(w), ncp = f^2 * n_total)

  expect_equal(out$F_critical, expected_crit, tolerance = 1e-12)
  expect_equal(out$power, expected_pow, tolerance = 1e-10)

})

test_that("pwrcontrast solves N (cal_n) for between-subjects and respects multiples of K", {
  w <- c(3, -1, -1, -1)
  alpha <- 0.05
  target_power <- 0.80

  out <- pwrcontrast(weight = w, paired = FALSE,
                     n_total = NULL, peta2 = 0.06, alpha = alpha, power = target_power)

  expect_s3_class(out, "cal_n")
  expect_equal(out$df_num, 1)
  # denominator df = n - K (between)
  expect_equal(out$df_denom, out$n_total - length(w))
  # N must be multiple of K
  expect_equal(out$n_total %% length(w), 0)
  # Achieved power should be >= target (by construction)
  expect_true(out$power >= target_power - 1e-12)
})

test_that("pwrcontrast solves alpha (cal_alpha) for paired design", {
  w <- c(1, 0, -1)
  n_total <- 30
  f <- 0.2
  target_power <- 0.9

  out <- pwrcontrast(weight = w, paired = TRUE,
                     n_total = n_total, cohensf = f, alpha = NULL, power = target_power)

  expect_s3_class(out, "cal_alpha")
  # denominator df = (n - 1)*(K - 1) in paired design
  expect_equal(out$df_denom, (n_total - 1) * (length(w) - 1))
  expect_true(out$alpha > 0 && out$alpha < 1)
})

test_that("pwrcontrast solves minimal detectable effect size (cal_es)", {
  w <- c(1, -1, 0)
  n_total <- 60
  alpha <- 0.05
  target_power <- 0.80

  out <- pwrcontrast(weight = w, paired = FALSE,
                     n_total = n_total, cohensf = NULL, peta2 = NULL,
                     alpha = alpha, power = target_power)

  expect_s3_class(out, "cal_es")
  expect_true(out$cohensf >= 0)
  expect_true(out$peta2 > 0 && out$peta2 < 1)

  # Back-check: with the returned effect, recompute power (should match out$power)
  df2 <- n_total - length(w)
  crit <- qf(1 - alpha, 1, df2)
  ncp  <- out$cohensf^2 * n_total
  expected_power <- 1 - pf(crit, 1, df2, ncp = ncp)
  expect_equal(out$power, expected_power, tolerance = 1e-5)
})

test_that("pwrcontrast input validation errors", {
  # weight checks
  expect_error(pwrcontrast(weight = NULL))
  expect_error(pwrcontrast(weight = 1))                      # length < 2
  expect_error(pwrcontrast(weight = c(1, NA)))               # NA in weight
  expect_error(pwrcontrast(weight = c(1, 1)))                # not sum-to-zero

  # paired must be length-1 logical
  expect_error(pwrcontrast(weight = c(1, -1), paired = c(TRUE, FALSE),
                           n_total = 20, cohensf = 0.2, alpha = 0.05, power = NULL))

  # nlim validity
  expect_error(pwrcontrast(weight = c(1, -1), nlim = 5,
                           n_total = NULL, cohensf = 0.2, alpha = 0.05, power = 0.8))
  expect_error(pwrcontrast(weight = c(1, -1), nlim = c(10, 5),
                           n_total = NULL, cohensf = 0.2, alpha = 0.05, power = 0.8))

  # exactly one of n_total, cohensf/peta2, alpha, power must be NULL
  # (two NULLs here: n_total and power)
  expect_error(pwrcontrast(weight = c(1, -1), paired = FALSE,
                           n_total = NULL, cohensf = 0.2, alpha = 0.05, power = NULL))

  # n_total multiple of K (between-subjects)
  expect_error(pwrcontrast(weight = c(1, 0, -1), paired = FALSE,
                           n_total = 20, cohensf = 0.2, alpha = 0.05, power = NULL))
})

test_that("pwrcontrast warns when both cohensf and peta2 are supplied and ignores peta2", {
  # cal_power branch; peta2 should be ignored and remain NA in output
  expect_warning({
    out <- pwrcontrast(weight = c(1, -1), paired = FALSE,
                       n_total = 40, cohensf = 0.25, peta2 = 0.06,
                       alpha = 0.05, power = NULL)
    expect_equal(out$cohensf, 0.25)
  })
})

test_that("pwrcontrast solve-N warns when target power not reached within nlim", {
  w <- c(1, -1)  # K=2
  # Very small effect + tight upper bound to force warning
  expect_warning({
    out <- pwrcontrast(weight = w, paired = FALSE,
                       n_total = NULL, cohensf = 0.01,
                       alpha = 0.05, power = 0.95,
                       nlim = c(2, 10))
    # For between-subjects, n candidates go by K steps; max candidate should be <= nlim[2]
    expect_true(out$n_total <= 10)
    expect_s3_class(out, "cal_n")
  })
})

test_that("pwrcontrast matches expected results", {
  path <- testthat::test_path("expected", "expected_pwrcontrast.csv")
  expected <- utils::read.csv(path)

  for(i in 1:nrow(expected)){
    e <- expected[i,]

    weight <- scan(text = e$weight, what = integer(), sep = ",", quiet = TRUE)
    paired <- e$paired
    n_total <- e$n_total
    alpha <- e$alpha
    power <- e$power
    cohensf <- e$cohensf

    # N

    res_n <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = NULL,
      alpha = alpha,
      power = power,
      cohensf = cohensf
    )[,c("n_total", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_n),
      as.numeric(e[,c("n_total", "ncp", "df_num", "df_denom")])
    )

    # Alpha

    res_alpha <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = n_total,
      alpha = NULL,
      power = power,
      cohensf = cohensf,
    )[,c("alpha", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_alpha),
      as.numeric(e[,c("alpha_est", "ncp", "df_num", "df_denom")]),
      tolerance = 10e-5
    )

    # Power

    res_power <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = n_total,
      alpha = alpha,
      power = NULL,
      cohensf = cohensf,
    )[,c("power", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_power),
      as.numeric(e[,c("power_est", "ncp", "df_num", "df_denom")]),
      tolerance = 10e-5
    )

    # Cohen's f

    res_cohensf <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = n_total,
      alpha = alpha,
      power = power,
      cohensf = NULL
    )[,c("cohensf", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_cohensf),
      as.numeric(e[,c("cohensf_est", "ncp_sen", "df_num", "df_denom")]),
      tolerance = 10e-5
    )
  }

})

test_that("pwrcontrast matches pwranova results", {
  set.seed(610)

  nrep <- 50

  for(i in 1:nrep){
    weight <- sample(c(-1,1),2, replace = FALSE) * sample(1:5, 1)
    paired <- sample(c(FALSE, TRUE), 1)
    n_total <- sample(5:100, 1) * 2
    alpha <- runif(1, 0.01,0.10)
    power <- runif(1,0.10,0.99)
    cohensf <- runif(1, 0.1,1)

    if(paired){
      nlevels_b <- NULL
      nlevels_w <- 2
    } else{
      nlevels_b <- 2
      nlevels_w <- NULL
    }

    # N

    res_n <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = NULL,
      alpha = alpha,
      power = power,
      cohensf = cohensf
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    res_n_anova <- pwranova(
      nlevels_b = nlevels_b,
      nlevels_w = nlevels_w,
      n_total = NULL,
      alpha = alpha,
      power = power,
      cohensf = cohensf
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    testthat::expect_equal(
      as.numeric(res_n),
      as.numeric(res_n_anova)
    )

    # Alpha

    res_alpha <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = n_total,
      alpha = NULL,
      power = power,
      cohensf = cohensf,
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    res_alpha_anova <- pwranova(
      nlevels_b = nlevels_b,
      nlevels_w = nlevels_w,
      n_total = n_total,
      alpha = NULL,
      power = power,
      cohensf = cohensf
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    testthat::expect_equal(
      as.numeric(res_alpha),
      as.numeric(res_alpha_anova)
    )

    # Power

    res_power <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = n_total,
      alpha = alpha,
      power = NULL,
      cohensf = cohensf,
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    res_power_anova <- pwranova(
      nlevels_b = nlevels_b,
      nlevels_w = nlevels_w,
      n_total = n_total,
      alpha = alpha,
      power = NULL,
      cohensf = cohensf
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    testthat::expect_equal(
      as.numeric(res_power),
      as.numeric(res_power_anova)
    )

    # Cohen's f

    res_cohensf <- pwrcontrast(
      weight = weight,
      paired = paired,
      n_total = n_total,
      alpha = alpha,
      power = power,
      cohensf = NULL,
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    res_cohensf_anova <- pwranova(
      nlevels_b = nlevels_b,
      nlevels_w = nlevels_w,
      n_total = n_total,
      alpha = alpha,
      power = power,
      cohensf = NULL
    )[,c("df_num", "df_denom", "n_total", "alpha", "power", "cohensf", "peta2", "F_critical", "ncp")]

    testthat::expect_equal(
      as.numeric(res_cohensf),
      as.numeric(res_cohensf_anova)
    )
  }

})

# ## Make expexted data
# path <- testthat::test_path("expected", "expected_pwrcontrast.csv")
# expected <- utils::read.csv(path)
# expected_out <- expected
#
# for(i in 1:nrow(expected)){
#   e <- expected[i,]
#
#   weight <- scan(text = e$weight, what = integer(), sep = ",", quiet = TRUE)
#
#   paired <- e$paired
#   alpha <- e$alpha
#   power <- e$power
#   cohensf <- e$cohensf
#
#   # solve N
#   res_n <- pwrcontrast(weight = weight, paired = paired, alpha = alpha, cohensf = cohensf, power = power)
#
#   n_total <- res_n$n_total
#   ncp <- res_n$ncp
#   df_num <- res_n$df_num
#   df_denom <- res_n$df_denom
#
#   # solve alpha
#   res_alpha <- pwrcontrast(weight = weight, paired = paired, n_total = n_total, cohensf = cohensf, power = power)
#   alpha_est <- res_alpha$alpha
#
#   # solve f
#   res_f <- pwrcontrast(weight = weight, paired = paired, n_total = n_total, alpha = alpha, power = power)
#   cohensf_est <- res_f$cohensf
#   ncp_sen <- res_f$ncp
#
#   # solve power
#   res_power <- pwrcontrast(weight = weight, paired = paired, n_total = n_total, alpha = alpha, cohensf = cohensf)
#   power_est <- res_power$power
#
#   expected_out[i,"n_total"] <- n_total
#   expected_out[i,"alpha_est"] <- alpha_est
#   expected_out[i,"power_est"] <- power_est
#   expected_out[i,"cohensf_est"] <- cohensf_est
#   expected_out[i,"ncp"] <- ncp
#   expected_out[i,"ncp_sen"] <- ncp_sen
#   expected_out[i,"df_num"] <- df_num
#   expected_out[i,"df_denom"] <- df_denom
#
# }
#
# clipr::write_clip(expected_out)
#
