test_that("cohensf_to_peta2 works for valid numeric inputs", {
  # Single value
  expect_equal(cohensf_to_peta2(0.25), 0.25^2 / (1 + 0.25^2))

  # Multiple values
  input <- c(0.1, 0.25, 0.4)
  expected <- input^2 / (1 + input^2)
  expect_equal(cohensf_to_peta2(input), expected)

  # Zero value (boundary case)
  expect_equal(cohensf_to_peta2(0), 0)

  # Return type should be numeric vector
  expect_type(cohensf_to_peta2(input), "double")
})

test_that("cohensf_to_peta2 throws errors for invalid inputs", {
  # Non-numeric input
  expect_error(cohensf_to_peta2("a"))
  expect_error(cohensf_to_peta2(list(0.1)))

  # NA input
  expect_error(cohensf_to_peta2(NA))
  expect_error(cohensf_to_peta2(c(0.1, NA)))

  # Negative values (invalid)
  expect_error(cohensf_to_peta2(-0.1))
  expect_error(cohensf_to_peta2(c(0.1, -0.2)))
})
