# Convert Cohen's f to Partial Eta Squared

Converts Cohen's f to partial eta squared (\\\eta_p^2\\) using the
standard definition in Cohen (1988).

## Usage

``` r
cohensf_to_peta2(f)
```

## Arguments

- f:

  A numeric vector of Cohen's f values. Each value must be greater than
  or equal to 0.

## Value

A numeric vector of partial eta squared values.

## Details

The conversion is defined as: \$\$\eta_p^2 = \frac{f^2}{1 + f^2}\$\$

This follows from the relationship: \\f = \sqrt{\eta_p^2 / (1 -
\eta_p^2)}\\

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

## See also

[peta2_to_cohensf](https://mutopsy.github.io/pwranova/reference/peta2_to_cohensf.md)

## Examples

``` r
# Convert a single Cohen's f value
cohensf_to_peta2(0.25)
#> [1] 0.05882353

# Convert multiple values
cohensf_to_peta2(c(0.1, 0.25, 0.4))
#> [1] 0.00990099 0.05882353 0.13793103
```
