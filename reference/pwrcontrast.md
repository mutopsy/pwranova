# Power Analysis for Planned Contrast in Between- or Within-Factor ANOVA

Computes power, required total sample size, alpha, or minimal detectable
effect size for a **single planned contrast** (1 df) in
between-participants or paired/repeated-measures settings.

## Usage

``` r
pwrcontrast(
  weight = NULL,
  paired = FALSE,
  n_total = NULL,
  cohensf = NULL,
  peta2 = NULL,
  alpha = NULL,
  power = NULL,
  nlim = c(2, 10000)
)
```

## Arguments

- weight:

  Numeric vector (length \\K \ge 2\\). Contrast weights whose sum must
  be zero.

- paired:

  Logical. `FALSE` for between-subjects (default), `TRUE` for
  paired/repeated-measures.

- n_total:

  Integer scalar. Total sample size. If `NULL`, the function solves for
  `n_total`.

- cohensf:

  Numeric (non-negative). Cohen's \\f\\. If `NULL`, it is derived from
  `peta2` when available.

- peta2:

  Numeric in \\(0,1)\\. Partial eta squared. If `NULL`, it is derived
  from `cohensf` when available.

- alpha:

  Numeric in \\(0,1)\\. If `NULL`, it is solved for.

- power:

  Numeric in \\(0,1)\\. If `NULL`, it is computed; if `n_total` is
  `NULL`, `n_total` is solved to achieve this power.

- nlim:

  Integer length-2. Search range of total `n` when solving sample size.

## Value

A one-row data frame with class:

- `"cal_power"` when power is calculated given `n_total`, `alpha`, and
  effect size;

- `"cal_n"` when `n_total` is solved;

- `"cal_alpha"` when `alpha` is solved;

- `"cal_es"` when minimal detectable effect sizes are solved.

Columns: `term` (always `"contrast"`), `weight` (comma-separated
string), `df_num`, `df_denom`, `n_total`, `alpha`, `power`, `cohensf`,
`peta2`, `F_critical`, `ncp`.

## Details

For a contrast with weights \\w_1, \dots, w_K\\ that sum to zero, the
numerator df is 1. The denominator df is \\n - K\\ for between-subjects
(unpaired) designs and \\n - 1\\ for paired/repeated-measures designs.
Power uses the noncentral *F*-with \\\lambda = f^2 \cdot
n\_{\mathrm{total}}\\.

- Contrast weights (`weight`) are not centered internally; only the
  zero-sum condition is enforced (up to numerical tolerance).

- When `paired = FALSE`, the total sample size `n_total` must be a
  multiple of the number of contrast groups \\K\\.

- Exactly one of `n_total`, an effect-size specification
  (`cohensf`/`peta2`), `alpha`, or `power` must be `NULL`; that quantity
  is then solved.

- Critical values are computed from the central *F*-distribution; power
  is based on the noncentral *F*-distribution with noncentrality
  parameter \\\lambda = f^2 \cdot n\_{\mathrm{total}}\\.

- Effect-size inputs can be given as Cohen’s \\f\\ or partial
  eta-squared \\\eta_p^2\\ (internally converted via \\f =
  \sqrt{\eta_p^2/(1-\eta_p^2)}\\). If both are `NULL`, the minimal
  detectable effect size is solved for given `n_total`, `alpha`, and
  `power`.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

## Examples

``` r
# Two-group contrast (1, -1), between-subjects: compute power
pwrcontrast(weight = c(1, -1), paired = FALSE,
            n_total = 40, cohensf = 0.25, alpha = 0.05)
#>       term weight df_num df_denom n_total alpha    power cohensf      peta2
#> 1 contrast   1,-1      1       38      40  0.05 0.337939    0.25 0.05882353
#>   F_critical ncp
#> 1   4.098172 2.5
#> Power (1 - beta) was calculated based on the total sample size, effect size, and alpha.

# Four-level contrast (e.g., Helmert-like), solve required N for target power
pwrcontrast(weight = c(3, -1, -1, -1), paired = FALSE,
            n_total = NULL, peta2 = 0.06, alpha = 0.05, power = 0.80)
#>       term     weight df_num df_denom n_total alpha     power   cohensf peta2
#> 1 contrast 3,-1,-1,-1      1      124     128  0.05 0.8095375 0.2526456  0.06
#>   F_critical      ncp
#> 1    3.91755 8.170213
#> The required total sample size was calculated based on the effect size, alpha, and power.
#> Note: 'power' indicates the achieved power rather than the target power.

# Paired contrast across K=3 conditions
pwrcontrast(weight = c(1, 0, -1), paired = TRUE,
            n_total = NULL, cohensf = 0.2, alpha = 0.05, power = 0.9)
#>       term weight df_num df_denom n_total alpha     power cohensf      peta2
#> 1 contrast 1,0,-1      1      264     265  0.05 0.9004175     0.2 0.03846154
#>   F_critical  ncp
#> 1   3.876924 10.6
#> The required total sample size was calculated based on the effect size, alpha, and power.
#> Note: 'power' indicates the achieved power rather than the target power.
```
