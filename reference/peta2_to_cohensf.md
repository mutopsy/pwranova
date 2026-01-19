# Convert Partial Eta Squared to Cohen's f

Converts partial eta squared (\\\eta_p^2\\) to Cohen's f using the
standard definition in Cohen (1988).

## Usage

``` r
peta2_to_cohensf(peta2)
```

## Arguments

- peta2:

  A numeric vector of partial eta squared values. Each value must be
  within the range of 0 to 1.

## Value

A numeric vector of Cohen's f values.

## Details

The conversion is defined as: \$\$f = \sqrt{\eta_p^2 / (1 -
\eta_p^2)}\$\$

This follows from the inverse relationship: \$\$\eta_p^2 =
\frac{f^2}{1 + f^2}\$\$

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

## See also

[cohensf_to_peta2](https://mutopsy.github.io/pwranova/reference/cohensf_to_peta2.md)

## Examples

``` r
# Convert a single partial eta squared value
peta2_to_cohensf(0.06)
#> [1] 0.2526456

# Convert multiple values
peta2_to_cohensf(c(0.01, 0.06, 0.14))
#> [1] 0.1005038 0.2526456 0.4034733
```
