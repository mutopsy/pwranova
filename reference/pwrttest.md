# Power Analysis for One-/Two-Sample and Paired t Tests

Computes statistical power, required total sample size, \\\alpha\\, or
the minimal detectable effect size for a *t*-test in one of three
designs: one-sample, two-sample (independent), or paired/repeated
measures.

## Usage

``` r
pwrttest(
  paired = FALSE,
  onesample = FALSE,
  n_total = NULL,
  alpha = NULL,
  power = NULL,
  delta = NULL,
  cohensf = NULL,
  peta2 = NULL,
  alternative = c("two.sided", "one.sided"),
  nlim = c(2, 10000)
)
```

## Arguments

- paired:

  Logical. `FALSE` for two-sample (independent; default), `TRUE` for
  paired/repeated-measures. Ignored when `onesample = TRUE`.

- onesample:

  Logical. `TRUE` for the one-sample *t*-test; if `TRUE`, `paired` is
  ignored.

- n_total:

  Integer scalar. Total sample size. If `NULL`, the function solves for
  `n_total`.

- alpha:

  Numeric in \\(0,1)\\. If `NULL`, it is solved for given the other
  inputs.

- power:

  Numeric in \\(0,1)\\. If `NULL`, it is computed; if `n_total` is
  `NULL`, `n_total` is solved to attain this power.

- delta:

  Numeric. Cohen's \\d\\-type effect size. If negative, it is converted
  to its absolute value. If `NULL`, it is derived from `cohensf` or
  `peta2` when available. If all three effect-size arguments (`delta`,
  `cohensf`, `peta2`) are `NULL`, then the effect size is treated as the
  unknown quantity and is solved for given `n_total`, `alpha`, and
  `power`. The exact definition depends on the design:

  - *One-sample*: Cohen's \\d = (\mu - \mu_0)/\sigma\\.

  - *Paired*: Cohen's \\d_z = \bar{d}/s_d\\, i.e., the mean of the
    difference scores divided by their standard deviation.

  - *Two-sample (equal allocation)*: Cohen's \\d\\ is defined as the
    mean difference divided by the pooled standard deviation; internally
    related to \\f\\ via \\d = 2f\\.

- cohensf:

  Numeric (non-negative). Cohen's \\f\\. If `NULL`, it can be derived
  from `delta`; if `delta` is supplied, `cohensf` is ignored.
  Effect-size relations by design:

  - Two-sample (equal allocation): \\d = 2f\\

  - Paired: \\d_z = f\\

  - One-sample: \\f\\ and \\\eta_p^2\\ are not supported

- peta2:

  Numeric in \\(0,1)\\. Partial eta squared. If `NULL`, it can be
  derived from `cohensf`; if `delta` is supplied, `peta2` is ignored.
  Not defined for one-sample designs.

- alternative:

  Character. Either `"two.sided"` or `"one.sided"`.

- nlim:

  Integer vector of length 2. Search range of total `n` when solving
  sample size.

## Value

A one-row `data.frame` with class `"cal_power"`, `"cal_n"`,
`"cal_alpha"`, or `"cal_es"`, depending on the solved quantity. Columns:
`df`, `n_total`, `alpha`, `power`, `delta`, `cohensf`, `peta2`,
`t_critical`, `ncp`.

## Details

- The sign of `delta` is ignored; its absolute value is used, because
  statistical power depends on the magnitude of the effect rather than
  its direction.

- If multiple effect-size arguments are supplied (`delta`, `cohensf`,
  `peta2`), precedence is `delta` \\\>\\ `cohensf` \\\>\\ `peta2`; the
  rest are ignored with a warning.

- For the two-sample design, equal allocation is assumed; `n_total` must
  be even when provided, and the solved `n_total` will be an even
  number.

- For the paired design, the effect size is interpreted as \\d_z\\.

- Computations use the central and noncentral *t*-distributions
  ([`stats::qt`](https://rdrr.io/r/stats/TDist.html),
  [`stats::pt`](https://rdrr.io/r/stats/TDist.html)); root finding uses
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) where
  needed.

- Results have been validated to match those produced by G\*Power for
  equivalent one-sample, paired, and two-sample *t* tests.

## Examples

``` r
# (1) Two-sample (independent), compute power given N and d
pwrttest(paired = FALSE, onesample = FALSE, alternative = "two.sided",
         n_total = 128, delta = 0.50, alpha = 0.05)
#>    df n_total alpha     power delta cohensf      peta2 t_critical      ncp
#> 1 126     128  0.05 0.8014596   0.5    0.25 0.05882353   1.978971 2.828427
#> Power (1 - beta) was calculated based on the total sample size, effect size, and alpha.

# (2) Paired t-test, solve required N for target power
pwrttest(paired = TRUE, onesample = FALSE, alternative = "one.sided",
         n_total = NULL, delta = 0.40, alpha = 0.05, power = 0.90)
#>   df n_total alpha     power delta_z cohensf    peta2 t_critical      ncp
#> 1 54      55  0.05 0.9004524     0.4     0.4 0.137931   1.673565 2.966479
#> The required total sample size was calculated based on the effect size, alpha, and power.
#> Note: 'power' indicates the achieved power rather than the target power.

# (3) One-sample t-test, solve alpha given N and power
pwrttest(onesample = TRUE, alternative = "two.sided",
         n_total = 40, delta = 0.40, alpha = NULL, power = 0.80)
#>   df n_total     alpha power delta cohensf peta2 t_critical      ncp
#> 1 39      40 0.1001708   0.8   0.4      NA    NA   1.683997 2.529822
#> Alpha was calculated based on the total sample size, effect size, and power.

# (4) Two-sample, specify effect via f or partial eta^2 (converted internally)
pwrttest(paired = FALSE, cohensf = 0.25, n_total = NULL, alpha = 0.05, power = 0.80)
#>    df n_total alpha     power delta cohensf      peta2 t_critical      ncp
#> 1 126     128  0.05 0.8014596   0.5    0.25 0.05882353   1.978971 2.828427
#> The required total sample size was calculated based on the effect size, alpha, and power.
#> Note: 'power' indicates the achieved power rather than the target power.
```
