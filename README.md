# pwranova: Power Analysis of Flexible ANOVA Designs and Related Tests

`pwranova` is an R package for power analysis in ANOVA designs, 
including between-, within-, and mixed-factor designs, 
with full support for both **main effects and interactions** across any number of factors.

The package allows calculation of **statistical power**, 
**required total sample size**, **significance level**, and **minimal 
detectable effect sizes** expressed as partial eta squared or Cohen's *f* 
for ANOVA terms and planned contrasts. 
In addition, complementary functions are included for common related 
tests such as *t*-tests and correlation tests, making the package a 
convenient toolkit for power analysis in experimental psychology 
and related fields.

## Links

<!-- - CRAN: <https://CRAN.R-project.org/package=pwranova>  -->
- Documentation: <https://mutopsy.github.io/pwranova>  
- Source code: <https://github.com/mutopsy/pwranova>

## Installation

You can install `pwranova` from GitHub using `devtools`:  

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install pwranova
devtools::install_github("mutopsy/pwranova")
```

## Dependencies

- R (>= 4.1.0)

No heavy external dependencies are required for the core functionality.

## Quick start

```r
library(pwranova)

### 1) Compute power (given N, effect size, alpha)

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
# Show only a selected term with `target` if you want a compact output
res_power_mixed <- pwranova(
  nlevels_b = 2,
  nlevels_w = c(2, 3),
  n_total   = 30,
  cohensf   = 0.50,
  alpha     = 0.05,
  epsilon   = 1.00,
  target    = "B1:W2"  # show only 2x3 interaction of the between factor and the second within factor 
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
```

## Functions

Current functions include:

- `pwranova()` — Power analysis for between-/within-/mixed-factor ANOVA, covering all main effects and interactions.
- `pwrcontrast()` — Power analysis for a single planned contrast (1 df) in between-subjects or paired/repeated-measures designs.
- `pwrttest()` — Power analysis for *t*-tests (one-sample, paired, and two-sample).
- `pwrcortest()` — Power analysis for Pearson’s correlation (using either the *t*-distribution or Fisher’s *z*-transformation approach).

For full documentation, see the reference site (pkgdown):
https://mutopsy.github.io/pwranova/reference/


## Citation

Until it is formally published, please cite the following preprint:

Muto, H. (2025). pwranova: An R package for power analysis of flexible ANOVA designs and related tests. Jxiv.  https://doi.org/10.51094/jxiv.1555

## Version history

See the Changelog:
https://mutopsy.github.io/pwranova/news/

## License

GPL-3
