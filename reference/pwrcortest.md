# Power Analysis for Pearson Correlation

Computes statistical power, required total sample size, \\\alpha\\, or
the minimal detectable correlation coefficient for a Pearson correlation
test. Two computational methods are supported: exact noncentral *t*
(`method = "t"`) and Fisher's *z*-transformation with normal
approximation (`method = "z"`).

## Usage

``` r
pwrcortest(
  alternative = c("two.sided", "one.sided"),
  n_total = NULL,
  alpha = NULL,
  power = NULL,
  rho = NULL,
  method = c("t", "z"),
  bias_correction = FALSE,
  nlim = c(2, 10000)
)
```

## Arguments

- alternative:

  Character. Either `"two.sided"` or `"one.sided"`.

- n_total:

  Integer scalar. Total sample size (\\n\\). Must be \\\ge 3\\ for
  `method = "t"` and \\\ge 4\\ for `method = "z"`. If `NULL`, the
  function solves for `n_total`.

- alpha:

  Numeric in \\(0,1)\\. If `NULL`, it is solved for.

- power:

  Numeric in \\(0,1)\\. If `NULL`, it is computed; if `n_total` is
  `NULL`, `n_total` is solved to attain this power.

- rho:

  Numeric correlation coefficient in \\(-1,1)\\, nonzero. If `NULL`,
  `rho` is solved for given the other inputs.

- method:

  Character. Either `"t"` (noncentral *t*-distribution) or `"z"`
  (Fisher's *z* transformation with normal approximation).

- bias_correction:

  Logical. Applies only to `method = "z"`. If `TRUE`, uses the
  bias-corrected Fisher *z*-transformation \\z_p =
  \operatorname{atanh}(r) + r/(2(n-1))\\.

- nlim:

  Integer vector of length 2. Search range of `n_total` when solving
  sample size.

## Value

A one-row `data.frame` with class `"cal_power"`, `"cal_n"`,
`"cal_alpha"`, or `"cal_es"`, depending on the solved quantity. Columns:

- `df` (only for `method = "t"`)

- `n_total`, `alpha`, `power`

- `rho`, `t_critical` or `z_critical`

- `ncp` (noncentrality parameter or mean under the alternative: see
  Details)

## Details

- Exactly one of `n_total`, `rho`, `alpha`, or `power` must be `NULL`;
  that quantity is then solved.

- For `method = "t"`, computations are based on the noncentral
  *t*-distribution with noncentrality parameter \\\lambda =
  \tfrac{\rho}{\sqrt{1-\rho^2}} \sqrt{n}\\.

- For `method = "z"`, computations use Fisher's *z* transformation of
  the population correlation, \\z\_\rho = \operatorname{atanh}(\rho)\\.
  Let \\W = \sqrt{n-3}\\ z\\. Under the alternative hypothesis, \\W \sim
  \mathrm{Normal}(\mu,\\1)\\ with \\\mu = \sqrt{n-3}\\ z\_\rho\\. If
  `bias_correction = TRUE`, \\\rho\\ is first bias-corrected before
  applying Fisher's transform. Critical values are taken from the
  central normal distribution under \\H_0:\rho=0\\ (i.e., \\W \sim
  \mathrm{Normal}(0,1)\\ under the null). The returned `ncp` equals
  \\\mu\\.

- **Validation against GPower:** Results have been confirmed to match
  those produced by GPower for equivalent correlation tests using the
  noncentral *t*-distribution.

- *Note:* Results from `method = "z"` will not exactly match
  `pwr::pwr.r.test`, because `pwr` uses a hybrid approach combining the
  Fisher-*z* approximation with a *t*-based critical value.

## Examples

``` r
# (1) Compute power for rho = 0.3, N = 50, two-sided test
pwrcortest(alternative = "two.sided", n_total = 50, rho = 0.3, alpha = 0.05)
#>   df n_total alpha     power rho t_critical      ncp
#> 1 48      50  0.05 0.5867505 0.3   2.010635 2.223748
#> Power (1 - beta) was calculated based on the total sample size, effect size, and alpha.

# (2) Solve required N for target power, using Fisher-z method
pwrcortest(method = "z", rho = 0.2, alpha = 0.05, power = 0.8)
#>   n_total alpha     power rho z_critical     ncp
#> 1     194  0.05 0.8000666 0.2   1.959964 2.80182
#> The required total sample size was calculated based on the effect size, alpha, and power.
#> Note: 'power' indicates the achieved power rather than the target power.

# (3) Solve minimal detectable correlation
pwrcortest(n_total = 60, alpha = 0.05, power = 0.9, rho = NULL)
#>   df n_total alpha power       rho t_critical      ncp
#> 1 58      60  0.05   0.9 0.3915971   2.001717 3.296573
#> The minimal detectable Cohen's f and partial eta squared were calculated based on the total sample size, alpha, and power.
```
