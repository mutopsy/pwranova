#' Power Analysis for Planned Contrast in Between- or Within-Factor ANOVA
#'
#' Computes power, required total sample size, alpha, or minimal detectable
#' effect size for a \strong{single planned contrast} (1 df) in between-participants
#' or paired/repeated-measures settings.
#'
#' For a contrast with weights \eqn{w_1, \dots, w_K} that sum to zero,
#' the numerator df is 1. The denominator df is \eqn{n - K} for
#' between-subjects (unpaired) designs and \eqn{(n - 1)(K - 1)} for
#' paired/repeated-measures designs. Power uses the noncentral F with
#' \eqn{\lambda = f^2 \cdot n_{\mathrm{total}}}.
#'
#' @param weight Numeric vector (length \eqn{K \ge 2}). Contrast weights whose sum must be (approximately) zero.
#' @param paired Logical. \code{FALSE} for between-subjects (default), \code{TRUE} for paired/repeated-measures.
#' @param n_total Integer or integer vector. Total sample size(s). If \code{NULL}, the function solves for \code{n_total}.
#' @param cohensf Numeric (non-negative). Cohen's \eqn{f}. If \code{NULL}, it is derived from \code{peta2} when available.
#' @param peta2 Numeric in \eqn{(0,1)}. Partial eta squared. If \code{NULL}, it is derived from \code{cohensf} when available.
#' @param alpha Numeric in \eqn{(0,1)}. If \code{NULL}, it is solved for.
#' @param power Numeric in \eqn{(0,1)}. If \code{NULL}, it is computed; if \code{n_total} is \code{NULL}, \code{n_total} is solved to achieve this power.
#' @param nlim Integer length-2. Search range of total \code{n} when solving sample size.
#'
#' @details
#' - \code{weight} is not normalized internally; only the zero-sum condition is enforced (up to numerical tolerance).
#' - When \code{paired = FALSE}, \code{n_total} must be a multiple of \eqn{K}.
#' - The function enforces exactly one of \code{n_total}, \code{cohensf}/\code{peta2}, \code{alpha}, or \code{power} to be \code{NULL}.
#'
#' @return A one-row data frame with class:
#'   \itemize{
#'     \item \code{"cal_power"} when power is calculated given \code{n_total}, \code{alpha}, and effect size;
#'     \item \code{"cal_n"} when \code{n_total} is solved;
#'     \item \code{"cal_alpha"} when \code{alpha} is solved;
#'     \item \code{"cal_es"} when minimal detectable effect sizes are solved.
#'   }
#'   Columns: \code{term} (always \code{"contrast"}), \code{weight} (comma-separated string),
#'   \code{df_num}, \code{df_denom}, \code{n_total}, \code{alpha}, \code{power},
#'   \code{cohensf}, \code{peta2}, \code{F_critical}, \code{ncp}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.).
#' Hillsdale, NJ: Lawrence Erlbaum Associates.
#'
#' @examples
#' # Two-group contrast (1, -1), between-subjects: compute power
#' pwrcontrast(weight = c(1, -1), paired = FALSE,
#'             n_total = 40, cohensf = 0.25, alpha = 0.05)
#'
#' # Four-level contrast (e.g., Helmert-like), solve required N for target power
#' pwrcontrast(weight = c(3, -1, -1, -1), paired = FALSE,
#'             n_total = NULL, peta2 = 0.06, alpha = 0.05, power = 0.80)
#'
#' # Paired contrast across K=3 conditions
#' pwrcontrast(weight = c(1, 0, -1), paired = TRUE,
#'             n_total = NULL, cohensf = 0.2, alpha = 0.05, power = 0.9)
#'
#' @importFrom stats pf qf uniroot
#' @export
pwrcontrast <- function(
    weight = NULL, paired = FALSE,
    n_total = NULL, cohensf = NULL, peta2 = NULL, alpha = NULL, power = NULL,
    nlim = c(2, 10000)
) {
  ## -------- Initial checks & conversions --------
  if (is.null(weight)) stop("'weight' must be specified.")
  if (!is.numeric(weight) || length(weight) < 2L) {
    stop("'weight' must be a numeric vector of length >= 2.")
  }
  if (anyNA(weight)) stop("'weight' must not contain NA values.")

  # sum-to-zero with tolerance
  if (abs(sum(weight)) > 1e-10) {
    stop("The sum of 'weight' must be (approximately) zero.")
  }

  if (!is.logical(paired) || length(paired) != 1L) {
    stop("'paired' must be a single logical value.")
  }

  if (length(nlim) != 2L || any(nlim %% 1 != 0)) stop("'nlim' must be two integers.")
  if (nlim[1] < 2) stop("'nlim[1]' must be 2 or larger.")
  if (nlim[1] >= nlim[2]) stop("'nlim[1]' must be smaller than 'nlim[2]'.")

  if (!is.null(cohensf) && !is.null(peta2)) {
    peta2 <- NULL
    warning("Both 'cohensf' and 'peta2' were supplied; 'peta2' was ignored.")
  }

  if (is.null(cohensf) && !is.null(peta2)) {
    cohensf <- peta2_to_cohensf(peta2)
  }

  if (is.null(peta2) && !is.null(cohensf)) {
    if (anyNA(cohensf) || any(cohensf < 0)) {
      stop("'cohensf' must be non-missing and non-negative.")
    }
    peta2 <- cohensf_to_peta2(cohensf)
  }

  if (!is.null(alpha)) {
    if (length(alpha) != 1L) stop("'alpha' must be length 1.")
    if (alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1).")
  }
  if (!is.null(power)) {
    if (length(power) != 1L) stop("'power' must be length 1.")
    if (power <= 0 || power >= 1) stop("'power' must be in (0, 1).")
  }

  if ((is.null(n_total) + is.null(cohensf) + is.null(alpha) + is.null(power)) != 1) {
    stop("Exactly one of 'n_total', 'cohensf' (or 'peta2'), 'alpha', or 'power' must be NULL.")
  }

  K <- length(weight)
  df_num <- 1

  if (!is.null(n_total) && !paired) {
    if (any(n_total %% K != 0)) {
      stop(paste0("'n_total' must be a multiple of the number of groups = ", K, "."))
    }
  }

  ## -------- Result scaffold --------
  res <- data.frame(
    term      = "contrast",
    weight    = paste0(weight, collapse = ","),
    df_num    = df_num,
    df_denom  = NA_real_,
    n_total   = NA_real_,
    alpha     = NA_real_,
    power     = NA_real_,
    cohensf   = if (is.null(cohensf)) NA_real_ else cohensf,
    peta2     = if (is.null(peta2))   NA_real_ else peta2,
    F_critical = NA_real_,
    ncp       = NA_real_
  )

  if (!is.null(n_total)) res$n_total <- n_total
  if (!is.null(alpha))   res$alpha   <- alpha
  if (!is.null(power))   res$power   <- power

  ## -------- Helper: denominator df --------
  denom_df_fun <- function(n) {
    if (paired) {
      (n - 1) * (K - 1)
    } else {
      n - K
    }
  }

  ## -------- Power (given N) --------
  if (is.null(power)) {
    res$df_denom  <- denom_df_fun(res$n_total)
    res$F_critical <- qf(1 - res$alpha, res$df_num, res$df_denom)
    res$ncp       <- res$cohensf^2 * res$n_total
    res$power     <- 1 - pf(res$F_critical, res$df_num, res$df_denom, ncp = res$ncp)
    return(structure(res, class = c("cal_power", "data.frame")))
  }

  ## -------- Solve N (given target power) --------
  if (is.null(n_total)) {
    if (paired) {
      nmin <- ceiling(nlim[1])
      n_candi <- seq.int(nmin, nlim[2], by = 1L)
    } else {
      nmin <- ceiling(nlim[1] / K) * K
      if (nmin <= K) nmin <- K * 2L
      n_candi <- seq.int(nmin, nlim[2], by = K)
    }

    df_denom_candi  <- denom_df_fun(n_candi)
    F_critical_candi <- qf(1 - res$alpha, df_num, df_denom_candi)
    ncp_candi       <- res$cohensf^2 * n_candi
    power_candi     <- 1 - pf(F_critical_candi, df_num, df_denom_candi, ncp = ncp_candi)

    idx <- which(power_candi >= power)[1]
    if (is.na(idx)) {
      warning(paste0("Power did not reach ", power,
                     " within N <= ", max(n_candi), "; the maximal N was returned."))
      idx <- length(n_candi)
    }

    res$n_total   <- n_candi[idx]
    res$df_denom  <- df_denom_candi[idx]
    res$power     <- power_candi[idx]
    res$F_critical <- F_critical_candi[idx]
    res$ncp       <- ncp_candi[idx]

    return(structure(res, class = c("cal_n", "data.frame")))
  }

  ## -------- Solve alpha (given power) --------
  if (is.null(alpha)) {
    res$df_denom  <- denom_df_fun(res$n_total)
    res$ncp       <- res$cohensf^2 * res$n_total
    res$F_critical <- qf(1 - res$power, res$df_num, res$df_denom, ncp = res$ncp)
    res$alpha     <- 1 - pf(res$F_critical, res$df_num, res$df_denom)
    return(structure(res, class = c("cal_alpha", "data.frame")))
  }

  ## -------- Solve effect size (given N, alpha, power) --------
  if (is.null(cohensf) && is.null(peta2)) {
    res$df_denom  <- denom_df_fun(res$n_total)
    res$F_critical <- qf(1 - res$alpha, res$df_num, res$df_denom)

    froot <- function(x) 1 - pf(res$F_critical, res$df_num, res$df_denom, ncp = x) - res$power
    upper <- 100
    val_u <- froot(upper)
    while (val_u < 0 && upper < 1e6) {  # increase upper until achievable
      upper <- upper * 2
      val_u <- froot(upper)
    }
    res$ncp <- uniroot(froot, lower = 0, upper = upper)$root

    res$cohensf <- sqrt(res$ncp / res$n_total)
    res$peta2   <- cohensf_to_peta2(res$cohensf)
    return(structure(res, class = c("cal_es", "data.frame")))
  }
}
