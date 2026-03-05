# pwranova 1.1.1 (2026-03-05)

* Strengthened input validation across functions.
  Added checks for finite values, parameter bounds, and boundary cases
  in effect-size conversion utilities and power analysis functions.
* Clarified documentation for `pwrcortest()` and `pwrttest()` to state that
  signed effect-size inputs (`rho`, `delta`) are converted to their absolute
  values.

# pwranova 1.1.0 (2026-03-05)

* Added `ncp_scale` argument to `pwrcortest()` for alternative scaling of
  the noncentrality parameter in correlation power analysis.

# pwranova 1.0.3 (2026-01-19)

* Fixed a bug in the denominator df for paired/repeated-measures contrast tests.

# pwranova 1.0.2 (2025-10-01)

* Minor adjustments in DESCRIPTION and documentation links.  
* Small fixes in preparation for CRAN submission.

# pwranova 1.0.1 (2025-09-30)

* Removed obsolete comments in code.
* Updated citation information.
* Minor updates to DESCRIPTION.

# pwranova 1.0.0 (2025-09-28)

* First stable release of the package.

* Provides power analysis functions for:
  - Between-, within-, and mixed-factor ANOVA designs (main effects and interactions).
  - Planned contrasts with user-defined weights.
  - t-tests (one-sample, paired, and two-sample).
  - Pearson correlations (using either the t-distribution or Fisher’s z-transformation).

* Supports solving for:
  - Statistical power
  - Required total sample size
  - Significance level (alpha)
  - Minimal detectable effect size (Cohen’s f, partial eta squared, or Cohen’s d-type indices).

* Validated against G*Power for supported ANOVA designs, t-tests, and correlation tests.
