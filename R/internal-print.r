#' @keywords internal
#' @noRd
#' @exportS3Method print cal_power
print.cal_power <- function(x, ...) {
  NextMethod()
  cat("Power (1 - beta) was calculated based on the total sample size, effect size, and alpha.\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
#' @exportS3Method print cal_n
print.cal_n <- function(x, ...) {
  NextMethod()
  cat("The required total sample size was calculated based on the effect size, alpha, and power.\n")
  cat("Note: 'power' indicates the achieved power rather than the target power.\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
#' @exportS3Method print cal_ns
print.cal_ns <- function(x, ...) {
  NextMethod()
  cat("The required total sample size was calculated based on the effect size, alpha, and power.\n")
  cat("Note: 'power' indicates the achieved power rather than the target power.\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
#' @exportS3Method print cal_alpha
print.cal_alpha <- function(x, ...) {
  NextMethod()
  cat("Alpha was calculated based on the total sample size, effect size, and power.\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
#' @exportS3Method print cal_alphas
print.cal_alphas <- function(x, ...) {
  NextMethod()
  cat("Alphas were calculated based on the total sample size, effect size, and power.\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
#' @exportS3Method print cal_es
print.cal_es <- function(x, ...) {
  NextMethod()
  cat("The minimal detectable Cohen's f and partial eta squared were calculated based on the total sample size, alpha, and power.\n")
  invisible(x)
}
