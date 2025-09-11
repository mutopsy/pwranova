#' Convert Cohen's f to Partial Eta Squared
#'
#' This function converts Cohen's f to partial eta squared (\eqn{\eta_p^2})
#' using the method described in Cohen (1988). Partial eta squared is a
#' commonly reported effect size index in ANOVA.
#'
#' The conversion is defined as:
#' \deqn{\eta_p^2 = \frac{f^2}{1 + f^2}}
#'
#' @param f A numeric vector of Cohen's f values.
#'   Each value must be greater than or equal to 0.
#'
#' @return A numeric vector of partial eta squared values.
#'
#' @examples
#' # Convert a single Cohen's f value
#' cohensf_to_peta2(0.25)
#'
#' # Convert multiple values
#' cohensf_to_peta2(c(0.1, 0.25, 0.4))
#'
#' @seealso [peta2_to_cohensf]
#'
#' @references
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.).
#' Hillsdale, NJ: Lawrence Erlbaum Associates.
#'
#' @export

cohensf_to_peta2 <- function(f) {
  if (!is.numeric(f)) {
    stop("'f' must be numeric.")
  }
  if (anyNA(f)) {
    stop("'f' must not contain NA values.")
  }
  if (min(f) < 0) {
    stop("'f' must be greater than or equal to 0.")
  }
  f^2 / (1 + f^2)
}
