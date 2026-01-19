# Power Analysis for Between-, Within-, or Mixed-Factor ANOVA

Computes power, required total sample size, alpha, or minimal detectable
effect size for fixed-effects terms in between-/within-/mixed-factor
ANOVA designs.

## Usage

``` r
pwranova(
  nlevels_b = NULL,
  nlevels_w = NULL,
  n_total = NULL,
  alpha = NULL,
  power = NULL,
  cohensf = NULL,
  peta2 = NULL,
  epsilon = 1,
  target = NULL,
  max_nfactor = 6,
  nlim = c(2, 10000)
)
```

## Arguments

- nlevels_b:

  Integer scalar or vector. Numbers of levels for between-subjects
  factors. Omit or set `NULL` if there is no between-subjects factor.

- nlevels_w:

  Integer scalar or vector. Numbers of levels for within-subjects
  factors. Omit or set `NULL` if there is no within-subjects factor.

- n_total:

  Integer scalar or vector. Total sample size across all groups. If
  `NULL`, the function solves for `n_total`.

- alpha:

  Numeric in \\(0,1)\\. If `NULL`, it is solved for.

- power:

  Numeric in \\(0,1)\\. If `NULL`, it is computed; if `n_total` is
  `NULL`, `n_total` is solved to achieve this power.

- cohensf:

  Numeric (non-negative). Cohen's \\f\\. If `NULL`, it is derived from
  `peta2` when available. If both effect-size arguments (`cohensf` and
  `peta2`) are `NULL`, the effect size is treated as unknown and solved
  for given `n_total`, `alpha`, and `power`.

- peta2:

  Numeric in \\(0,1)\\. Partial eta squared. If `NULL`, it is derived
  from `cohensf` when available.

- epsilon:

  Numeric in \\(0,1\]\\. Nonsphericity parameter applied to
  within-subjects terms with \\\mathrm{df}\_1 \ge 2\\. Ignored if there
  is no within-subjects factor or if all within factors have two levels.

- target:

  Character vector of term labels to compute (e.g., `"B1"`, `"W1"`,
  `"B1:W1"`, ...). If `NULL`, all terms are returned.

- max_nfactor:

  Integer. Safety cap for the total number of factors.

- nlim:

  Integer length-2. Search range of total `n` when solving sample size.

## Value

A data frame with S3 class:

- `"cal_power"` when power is calculated given `n_total`, `alpha`, and
  effect size;

- `"cal_n"` or `"cal_ns"` when `n_total` is solved;

- `"cal_alpha"` or `"cal_alphas"` when `alpha` is solved;

- `"cal_es"` when minimal detectable effect sizes are solved.

Columns include `term`, `df_num`, `df_denom`, `n_total`, `alpha`,
`power`, `cohensf`, `peta2`, `F_critical`, `ncp`, `epsilon`.

## Details

- Fixed-effects, balanced designs are assumed. All groups/cells have
  equal cell sizes and effects are tested with standard fixed-effects
  ANOVA models.

- Numerator degrees of freedom for within-subjects terms with
  \\\mathrm{df}\_1 \ge 2\\ are adjusted by the nonsphericity parameter
  `epsilon`.

- Denominator degrees of freedom follow standard mixed-ANOVA formulas
  and are multiplied by the same `epsilon` for within-subjects terms.

- Critical values are computed from the central *F*-distribution; power
  uses the noncentral *F*-distribution with noncentrality parameter
  \\\lambda = f^2 \cdot n\_{\mathrm{total}}\\.

- Effect-size inputs can be given as Cohen’s \\f\\ or partial
  eta-squared \\\eta_p^2\\ (internally converted via \\f =
  \sqrt{\eta_p^2/(1-\eta_p^2)}\\). If both are `NULL`, the minimal
  detectable effect size is solved for given `n_total`, `alpha`, and
  `power`.

- Exactly one of `n_total`, an effect-size specification
  (`cohensf`/`peta2`), `alpha`, or `power` must be `NULL`; that quantity
  is then solved.

- **Validation against GPower:** For the subset of designs supported by
  GPower (between-, within-, and mixed-factor ANOVA with equal cell
  sizes), `pwranova()` was validated to produce results identical to
  those of GPower.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

## Examples

``` r
# One between factor (k = 3), one within factor (m = 4), compute power
pwranova(nlevels_b = 3, nlevels_w = 4, n_total = 60,
         cohensf = 0.25, alpha = 0.05, power = NULL, epsilon = 0.8)
#>    term df_num df_denom n_total alpha     power cohensf      peta2 F_critical
#> 1    B1    2.0     57.0      60  0.05 0.3744311    0.25 0.05882353   3.158843
#> 2    W1    2.4    136.8      60  0.05 0.3592414    0.25 0.05882353   2.876716
#> 3 B1:W1    4.8    136.8      60  0.05 0.2687648    0.25 0.05882353   2.307783
#>    ncp epsilon
#> 1 3.75     1.0
#> 2 3.75     0.8
#> 3 3.75     0.8
#> Power (1 - beta) was calculated based on the total sample size, effect size, and alpha.

# Solve required total N for target power
pwranova(nlevels_b = 2, nlevels_w = NULL, n_total = NULL,
         peta2 = 0.06, alpha = 0.05, power = 0.8)
#>   term df_num df_denom n_total alpha     power   cohensf peta2 F_critical
#> 1   B1      1      124     126  0.05 0.8034337 0.2526456  0.06    3.91755
#>        ncp epsilon
#> 1 8.042553       1
#> The required total sample size was calculated based on the effect size, alpha, and power.
#> Note: 'power' indicates the achieved power rather than the target power.
```
