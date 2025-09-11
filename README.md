# pwranova: Power Analysis for Flexible ANOVA Designs and Related Tests

pwranova is an R package for power analysis in ANOVA designs — including
between-, within-, and mixed-factor designs — with full support for both
main effects and interactions across any number of factors.

The package provides functions to compute statistical power, required
total sample size, significance level, and minimal detectable
effect sizes (partial eta squared or Cohen’s f) for ANOVA terms and
planned contrasts. In addition, complementary helpers for commonly used
tests (e.g., t-tests and correlation tests) are planned, making
pwranova a convenient toolkit for power analysis in experimental
psychology and related fields.

## Installation

Install the development version from GitHub:

  # Install devtools if not already installed
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }

  # Install pwranova
  devtools::install_github("mutopsy/pwranova")

Tip: If you build a pkgdown site, use
  devtools::install_github("mutopsy/pwranova", build_vignettes = TRUE)
to include vignettes.

## Dependencies

- R (>= 4.1.0)

No heavy external dependencies are required for the core functionality.

## Quick start

### 1) Compute power (given N, effect size, alpha)

  library(pwranova)

  # One between factor (3 levels), no within factor
  res_power_between <- pwranova(
    nlevels_b = 3,
    n_total   = 60,
    cohensf   = 0.25,
    alpha     = 0.05
  )
  res_power_between

  # One within factor (4 levels), no between factor
  # (Optionally set epsilon for sphericity correction; default is 1)
  res_power_within <- pwranova(
    nlevels_w = 4,
    n_total   = 30,
    cohensf   = 0.50,
    alpha     = 0.05,
    epsilon   = 1.00
  )
  res_power_within

  # Mixed design: one between factor (2 levels) and two within factors (2 and 3 levels)
  # Show only selected terms with `target` if you want a compact output
  res_power_mixed <- pwranova(
    nlevels_b = 2,
    nlevels_w = c(2, 3),
    n_total   = 30,
    cohensf   = 0.50,
    alpha     = 0.05,
    # epsilon applies to within terms with df1 >= 2 (here, W2 and terms including W2)
    epsilon   = 1.00,
    target    = c("B1", "W1", "W2", "B1:W2")  # example subset
  )
  res_power_mixed

### 2) Solve required total N (given target power)

  # One between factor (3 levels), no within factor
  res_n_between <- pwranova(
    nlevels_b = 3,
    cohensf   = 0.25,
    alpha     = 0.05,
    power     = 0.80
  )
  res_n_between  # returns the smallest total N (multiple of groups) meeting the target power

### 3) Planned contrast power

  # Contrast weights must sum to 0
  res_contrast <- pwrcontrast(
    weight  = c(1, -1, 0),  # three conditions, compare 1 vs 2
    paired  = FALSE,
    n_total = 60,
    cohensf = 0.25,
    alpha   = 0.05
  )
  res_contrast

## Functions

Current functions include:

- pwranova() — Power analysis for between-/within-/mixed-factor ANOVA, covering all main effects and interactions.
- pwrcontrast() — Power analysis for a single planned contrast (1 df) in between-subjects or paired/repeated-measures designs.

For full documentation, see the reference site (pkgdown):
https://mutopsy.github.io/pwranova/reference/

## Version history

See the Changelog:
https://mutopsy.github.io/pwranova/news/

## License

GPL-3 (see LICENSE).
