#' Convert Partial Eta Squared to Cohen's f
#'
#' Converts partial eta squared (\eqn{\eta_p^2}) to
#' Cohen's f using the standard definition in Cohen (1988).
#'
#' The conversion is defined as:
#' \deqn{f = \sqrt{\eta_p^2 / (1 - \eta_p^2)}}
#'
#' This follows from the inverse relationship:
#' \deqn{\eta_p^2 = \frac{f^2}{1 + f^2}}
#'
#' @param peta2 A numeric vector of partial eta squared values.
#'   Each value must be within the range of 0 to 1.
#'
#' @return A numeric vector of Cohen's f values.
#'
#' @examples
#' # Convert a single partial eta squared value
#' peta2_to_cohensf(0.06)
#'
#' # Convert multiple values
#' peta2_to_cohensf(c(0.01, 0.06, 0.14))
#'
#' @seealso [cohensf_to_peta2]
#'
#' @references
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.).
#' Hillsdale, NJ: Lawrence Erlbaum Associates.
#'
#' @importFrom stats pf qf uniroot
#' @export

peta2_to_cohensf <- function(peta2) {
  if (!is.numeric(peta2)) {
    stop("'peta2' must be numeric.")
  }
  if (anyNA(peta2)) {
    stop("'peta2' must not contain NA values.")
  }
  if (min(peta2) <= 0 || max(peta2) >= 1) {
    stop("'peta2' must be between 0 and 1.")
  }
  sqrt(peta2 / (1 - peta2))
}
