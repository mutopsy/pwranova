## Test environments
* Local: Windows 11, R 4.4.3, Rtools, x86_64
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 note

## Comments
This submission updates the CRAN version from 1.0.3 to 1.1.5.

Notable changes since the previous CRAN release include:
* Added the `ncp_scale` argument to `pwrcortest()` for an alternative scaling of the noncentrality parameter.
* Strengthened input validation across functions.
* Improved documentation, including additional explanations of effect-size definitions, noncentrality parameter calculations, and conventional effect-size benchmarks.
* Fixed an issue in `pwranova()` where epsilon adjustments were not fully incorporated into noncentrality parameter calculations.
* Expanded validation tests, including additional validation cases against G\*Power.
* Updated citation information following publication in JOSS.
