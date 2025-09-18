test_that("peta2_to_cohensf works for valid numeric inputs", {
  # Single value
  expect_equal(peta2_to_cohensf(0.06), sqrt(0.06 / (1 - 0.06)))

  # Multiple values
  input <- c(0.01, 0.06, 0.14)
  expected <- sqrt(input / (1 - input))
  expect_equal(peta2_to_cohensf(input), expected)

  # Return type should be numeric vector
  expect_type(peta2_to_cohensf(input), "double")
})

test_that("peta2_to_cohensf errors for non-numeric input", {
  # Input is not numeric
  expect_error(peta2_to_cohensf("a"))
  expect_error(peta2_to_cohensf(list(0.1)))
})

test_that("peta2_to_cohensf errors for NA input", {
  # Input contains NA
  expect_error(peta2_to_cohensf(NA))
  expect_error(peta2_to_cohensf(c(0.1, NA)))
})

test_that("peta2_to_cohensf errors for out-of-range input", {
  # Values <= 0
  expect_error(peta2_to_cohensf(0))
  expect_error(peta2_to_cohensf(-0.01))

  # Values >= 1
  expect_error(peta2_to_cohensf(1))
  expect_error(peta2_to_cohensf(1.01))
})
