# tests/testthat/test-pwranova.R
test_that("basic power calculation (between-subjects, 1 factor) works", {
  res <- pwranova(
    nlevels_b = 3, nlevels_w = NULL,
    n_total = 60, cohensf = 0.25, alpha = 0.05, power = NULL
  )
  expect_s3_class(res, "cal_power")
  expect_equal(nrow(res), 1L)
  expect_equal(res$term, "B1")
  expect_equal(res$df_num, 2)
  expect_equal(res$df_denom, 60 - 3)
  expect_true(res$power > 0 && res$power < 1)
})

test_that("pwranova works with two between-subject factors (2 x 3 design)", {
  res <- pwranova(
    nlevels_b = c(2, 3),
    nlevels_w = NULL,
    n_total   = 60,
    cohensf   = 0.25,
    alpha     = 0.05
  )

  expect_s3_class(res, "cal_power")
  expect_equal(nrow(res), 3)  # B1, B2, B1:B2 の3つ

  expect_equal(res$df_num[res$term == "B1"], 1) # (2-1)
  expect_equal(res$df_num[res$term == "B2"], 2) # (3-1)
  expect_equal(res$df_num[res$term == "B1:B2"], 2) # (2-1)*(3-1)

  expect_true(all(res$df_denom == 54))

  expect_true(all(res$power > 0 & res$power < 1))
})



test_that("mixed design returns terms and epsilon behaves", {
  res <- pwranova(
    nlevels_b = 2, nlevels_w = 3,
    n_total = 40, cohensf = 0.25, alpha = 0.05,
    target = c("B1", "W1", "B1:W1"), epsilon = 0.8
  )
  expect_setequal(res$term, c("B1", "W1", "B1:W1"))
  # B1 is between-only -> epsilon should be 1 for that term
  expect_equal(res$epsilon[res$term == "B1"], 1)
  # W1 has df1 = 2; epsilon should apply
  expect_equal(res$df_num[res$term == "W1"], (3 - 1) * 0.8, tolerance = 1e-12)
  # interaction W1:B1 also includes within, so epsilon applies to numerator df
  expect_equal(res$df_num[res$term == "B1:W1"], (1 * (3 - 1)) * 0.8, tolerance = 1e-12)
  # powers are numeric and in bounds
  expect_true(all(res$power > 0 & res$power < 1))
})

test_that("epsilon is ignored/reset when no/2-level within factors", {
  # No within factor -> epsilon ignored, defaulting to 1
  res1 <- pwranova(
    nlevels_b = 2, nlevels_w = NULL,
    n_total = 20, cohensf = 0.25, alpha = 0.05,
    epsilon = 0.7
  )
  expect_true(all(res1$epsilon == 1))

  # Within factor has only 2 levels -> epsilon reset to 1
  res2 <- pwranova(
    nlevels_b = 2, nlevels_w = 2,
    n_total = 20, cohensf = 0.25, alpha = 0.05,
    epsilon = 0.6
  )
  expect_true(all(res2$epsilon == 1))
})

test_that("solve required N (between-subjects, 2 groups) works and is multiple of groups", {
  res <- pwranova(
    nlevels_b = 2, nlevels_w = NULL,
    n_total = NULL, cohensf = 0.25, alpha = 0.05, power = 0.80,
    nlim = c(2, 2000)
  )
  expect_s3_class(res, "cal_n")
  expect_true(res$n_total %% 2 == 0)
  expect_true(res$power >= 0.80)
})

test_that("solve alpha given N and power works", {
  res <- pwranova(
    nlevels_b = 3, nlevels_w = NULL,
    n_total = 90, cohensf = 0.2, alpha = NULL, power = 0.90
  )
  expect_s3_class(res, "cal_alpha")
  expect_true(res$alpha > 0 && res$alpha < 1)
})

test_that("solve minimal detectable effect size (cohensf, peta2) works", {
  res <- pwranova(
    nlevels_b = 2, nlevels_w = 3,
    n_total = 48, cohensf = NULL, peta2 = NULL,
    alpha = 0.05, power = 0.80, epsilon = 0.85,
    target = c("W1")
  )
  expect_s3_class(res, "cal_es")
  expect_true(res$cohensf > 0)
  expect_true(res$peta2 > 0 && res$peta2 < 1)
})

test_that("target filtering works and errors on invalid target", {
  res <- pwranova(
    nlevels_b = 2, nlevels_w = 3,
    n_total = 48, cohensf = 0.2, alpha = 0.05,
    target = c("B1", "W1")
  )
  expect_setequal(res$term, c("B1", "W1"))
  # invalid target name should error
  expect_error(
    pwranova(
      nlevels_b = 2, nlevels_w = 3,
      n_total = 48, cohensf = 0.2, alpha = 0.05,
      target = "Z9"
    ),
    "contains no valid"
  )
})

test_that("input validation catches common issues", {
  # neither between nor within specified
  expect_error(
    pwranova(
      nlevels_b = NULL, nlevels_w = NULL,
      n_total = 10, cohensf = 0.2, alpha = 0.05
    ),
    "must be specified"
  )
  # n_total not multiple of groups
  expect_error(
    pwranova(
      nlevels_b = c(2, 2), nlevels_w = NULL,
      n_total = 10, cohensf = 0.2, alpha = 0.05
    ),
    "multiple of the number of groups"
  )
  # epsilon out of range
  expect_error(
    pwranova(
      nlevels_b = 2, nlevels_w = 3,
      n_total = 24, cohensf = 0.2, alpha = 0.05,
      epsilon = 1.2
    ),
    "in \\(0, 1\\]"
  )
  # exactly one of n_total/cohensf/alpha/power must be NULL
  expect_error(
    pwranova(
      nlevels_b = 2, nlevels_w = 3,
      n_total = 24, cohensf = 0.2, alpha = 0.05, power = 0.8
    ),
    "Exactly one"
  )
})


test_that("pwranova matches G*Power results", {
  path <- testthat::test_path("expected", "expected_pwranova.csv")
  expected <- utils::read.csv(path)
  i <- 3
  for(i in 1:nrow(expected)){
    e <- expected[i,]

    if(e$nlevels_b == "NULL"){
      nlevels_b <- NULL
    } else{
      nlevels_b <- scan(text = e$nlevels_b, what = integer(), sep = ",", quiet = TRUE)
    }

    if(e$nlevels_w == "NULL"){
      nlevels_w <- NULL
    } else{
      nlevels_w <- scan(text = e$nlevels_w, what = integer(), sep = ",", quiet = TRUE)
    }

    target <- e$target
    n_total <- e$n_total
    alpha <- e$alpha
    power <- e$power
    cohensf <- e$cohensf

    # N

    res_n <- pwranova(
      nlevels_b = nlevels_b, nlevels_w = nlevels_w,
      n_total = NULL,
      alpha = alpha,
      power = power,
      cohensf = cohensf,
      target = target
    )[,c("n_total", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_n),
      as.numeric(e[,c("n_total", "ncp", "df_num", "df_denom")])
    )

    # Alpha

    res_alpha <- pwranova(
      nlevels_b = nlevels_b, nlevels_w = nlevels_w,
      n_total = n_total,
      alpha = NULL,
      power = power,
      cohensf = cohensf,
      target = target
    )[,c("alpha", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_alpha),
      as.numeric(e[,c("alpha_est", "ncp", "df_num", "df_denom")]),
      tolerance = 10e-5
    )

    # Power

    res_power <- pwranova(
      nlevels_b = nlevels_b, nlevels_w = nlevels_w,
      n_total = n_total,
      alpha = alpha,
      power = NULL,
      cohensf = cohensf,
      target = target
    )[,c("power", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_power),
      as.numeric(e[,c("power_est", "ncp", "df_num", "df_denom")]),
      tolerance = 10e-5
    )

    # Cohen's f

    res_cohensf <- pwranova(
      nlevels_b = nlevels_b, nlevels_w = nlevels_w,
      n_total = n_total,
      alpha = alpha,
      power = power,
      cohensf = NULL,
      target = target
    )[,c("cohensf", "ncp", "df_num", "df_denom")]

    testthat::expect_equal(
      as.numeric(res_cohensf),
      as.numeric(e[,c("cohensf_est", "ncp_sen", "df_num", "df_denom")]),
      tolerance = 10e-5
    )




  }


})
