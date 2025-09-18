#' Power Analysis for One-/Two-Sample and Paired t Tests
#'
#' Computes statistical \strong{power}, required total sample size, \eqn{\alpha},
#' or the minimal detectable effect size for a \emph{t} test in one of three designs:
#' one-sample, two-sample (independent), or paired/repeated measures.
#'
#' @param paired Logical. \code{FALSE} for two-sample (independent; default),
#'   \code{TRUE} for paired/repeated-measures. Ignored when \code{onesample = TRUE}.
#' @param onesample Logical. \code{TRUE} for the one-sample t test; if \code{TRUE},
#'   \code{paired} is ignored.
#' @param alternative Character. Either \code{"two.sided"} or \code{"one.sided"}.
#' @param n_total Integer or integer vector. Total sample size.
#'   If \code{NULL}, the function solves for \code{n_total}.
#' @param delta Numeric (non-negative). Cohen's \eqn{d}-type effect size.
#'   If \code{NULL}, it is derived from \code{cohensf} or \code{peta2} when available.
#'   (Two-sample, equal allocation relation: \eqn{\delta = 2f}.)
#' @param cohensf Numeric (non-negative). Cohen's \eqn{f}.
#'   If \code{NULL}, it can be derived from \code{delta}; if \code{delta} is supplied,
#'   \code{cohensf} is ignored.
#' @param peta2 Numeric in \eqn{(0,1)}. Partial eta squared.
#'   If \code{NULL}, it can be derived from \code{cohensf}; if \code{delta} is supplied,
#'   \code{peta2} is ignored.
#' @param alpha Numeric in \eqn{(0,1)}. If \code{NULL}, it is solved for given the other inputs.
#' @param power Numeric in \eqn{(0,1)}. If \code{NULL}, it is computed; if \code{n_total} is \code{NULL},
#'   \code{n_total} is solved to attain this power.
#' @param nlim Integer vector of length 2. Search range of total \code{n} when solving sample size.
#'
#' @details
#' \itemize{
#'   \item If multiple effect-size arguments are supplied (\code{delta}, \code{cohensf}, \code{peta2}),
#'         precedence is \code{delta} \eqn{>} \code{cohensf} \eqn{>} \code{peta2}; the rest are ignored with a warning.
#'   \item For the two-sample design, equal allocation is assumed; \code{n_total} must be even when provided,
#'         and the solved \code{n_total} will be an even number.
#'   \item Computations use the central and noncentral \emph{t} distributions (\code{stats::qt}, \code{stats::pt});
#'         root finding uses \code{stats::uniroot()} where needed.
#' }
#'
#' @return A one-row \code{data.frame} with class
#'   \code{"cal_power"}, \code{"cal_n"}, \code{"cal_alpha"}, or \code{"cal_es"},
#'   depending on the solved quantity. Columns:
#'   \code{df}, \code{n_total}, \code{alpha}, \code{power},
#'   \code{delta}, \code{cohensf}, \code{peta2},
#'   \code{t_critical}, \code{ncp}.
#'
#' @examples
#' # (1) Two-sample (independent), compute power given N and d
#' pwrttest(paired = FALSE, onesample = FALSE, alternative = "two.sided",
#'          n_total = 128, delta = 0.50, alpha = 0.05)
#'
#' # (2) Paired t test, solve required N for target power
#' pwrttest(paired = TRUE, onesample = FALSE, alternative = "one.sided",
#'          n_total = NULL, delta = 0.40, alpha = 0.05, power = 0.90)
#'
#' # (3) One-sample t test, solve alpha given N and power
#' pwrttest(onesample = TRUE, alternative = "two.sided",
#'          n_total = 40, delta = 0.40, alpha = NULL, power = 0.80)
#'
#' # (4) Two-sample, specify effect via f or partial eta^2 (converted internally)
#' pwrttest(paired = FALSE, cohensf = 0.25, n_total = NULL, alpha = 0.05, power = 0.80)
#'
#' @importFrom stats pt qt uniroot
#' @export
pwrttest <- function(
    paired = FALSE, onesample = FALSE, alternative = c("two.sided", "one.sided"),
    n_total = NULL, delta = NULL, cohensf = NULL, peta2 = NULL, alpha = NULL, power = NULL,
    nlim = c(2, 10000)
) {
  ## -------- Initial checks & conversions --------

  alternative <- alternative[1]

  if (!is.logical(paired) || length(paired) != 1L) {
    stop("'paired' must be a single logical value.")
  }

  if(!alternative %in% c("two.sided", "one.sided")){
    stop("'alternative' must be 'two.sided' or 'one.sided'.")
  }

  if (length(nlim) != 2L || any(nlim %% 1 != 0)) stop("'nlim' must be two integers.")
  if (nlim[1] < 2) stop("'nlim[1]' must be 2 or larger.")
  if (nlim[1] >= nlim[2]) stop("'nlim[1]' must be smaller than 'nlim[2]'.")

  if (!is.null(delta) && !is.null(cohensf) && !is.null(peta2)) {
    cohensf <- NULL
    peta2 <- NULL
    warning("All of 'delta', 'cohensf' and 'peta2' were supplied; cohensf and 'peta2' were ignored.")
  } else if(!is.null(cohensf) && !is.null(peta2)){
    peta2 <- NULL
    warning("Both 'cohensf' and 'peta2' were supplied; 'peta2' was ignored.")
  } else if(!is.null(delta) && !is.null(cohensf)){
    cohensf <- NULL
    warning("Both 'delta' and 'cohensf' were supplied; 'cohensf' was ignored.")
  } else if(!is.null(delta) && !is.null(peta2)){
    peta2 <- NULL
    warning("Both 'delta' and 'peta2' were supplied; 'peat2' was ignored.")
  }

  if(!is.null(delta)){
    delta <- abs(delta)
    cohensf <- abs(delta) / 2
    peta2 <- cohensf_to_peta2(cohensf)
  } else if(!is.null(cohensf)){
    peta2 <- cohensf_to_peta2(cohensf)
    delta <- cohensf * 2
  } else if(!is.null(peta2)){
    cohensf <- peta2_to_cohensf(peta2)
    delta <- cohensf * 2
  }

  if (!is.null(alpha)) {
    if (length(alpha) != 1L) stop("'alpha' must be length 1.")
    if (alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1).")
  }
  if (!is.null(power)) {
    if (length(power) != 1L) stop("'power' must be length 1.")
    if (power <= 0 || power >= 1) stop("'power' must be in (0, 1).")
  }

  if ((is.null(n_total) + is.null(delta) + is.null(alpha) + is.null(power)) != 1) {
    stop("Exactly one of 'n_total', 'delta' (or 'cohensf' or 'peta2'), 'alpha', or 'power' must be NULL.")
  }

  df <- NA_real_

  if(onesample){
    if(paired){
      warning("Because 'onesample' was true, 'paired' was ignored.")
      paired <- FALSE
    }
    design <- "one.sample"
    if(!is.null(n_total)){
      df <- n_total - 1
    }
  } else if(paired){
    design <- "paired"
    if(!is.null(n_total)){
      df <- n_total - 1
    }
  } else{
    design <- "two.sample"
    if(!is.null(n_total)){
      df <- n_total - 2
    }
  }

  if(alternative == "two.sided"){
    divider <- 2
  } else{
    divider <- 1
  }

  if (!is.null(n_total) && design == "two.sample") {
    if (any(n_total %% 2 != 0)) {
      stop(paste0("'n_total' must be a multiple of the number of groups = 2."))
    }
  }

  ## -------- Result scaffold --------
  res <- data.frame(
    df        = df,
    n_total   = NA_real_,
    alpha     = NA_real_,
    power     = NA_real_,
    delta     = if (is.null(delta)) NA_real_ else delta,
    cohensf   = if (is.null(cohensf)) NA_real_ else cohensf,
    peta2     = if (is.null(peta2))   NA_real_ else peta2,
    t_critical = NA_real_,
    ncp       = NA_real_
  )

  if (!is.null(n_total)) res$n_total <- n_total
  if (!is.null(alpha))   res$alpha   <- alpha
  if (!is.null(power))   res$power   <- power

  ## -------- Power (given N) --------
  if (is.null(power)) {
    res$t_critical <- qt(1 - res$alpha / divider, res$df)
    res$ncp <- res$delta * sqrt(res$n_total)
    if(design == "two.sample"){
      res$ncp <- res$ncp / 2
    }
    res$power <- 1 - pt(res$t_critical, res$df, ncp = res$ncp)
    if(alternative == "two.sided"){
      res$power <- res$power + pt(-res$t_critical, res$df, ncp = res$ncp)
    }

    return(structure(res, class = c("cal_power", "data.frame")))
  }

  ## -------- Solve N (given target power) --------
  if (is.null(n_total)) {
    if (design != "two.sample") {
      nmin <- ceiling(nlim[1])
      n_candi <- seq.int(nmin, nlim[2], by = 1L)
      df_candi <- n_candi - 1
    } else {
      nmin <- ceiling(nlim[1] / 2) * 2
      if (nmin <= 2) nmin <- 4
      n_candi <- seq.int(nmin, nlim[2], by = 2)
      df_candi <- n_candi - 2
    }

    t_critical_candi <- qt(1 - res$alpha/divider, df_candi)
    ncp_candi <- res$delta * sqrt(n_candi)
    if(design == "two.sample"){
      ncp_candi <- ncp_candi / 2
    }

    power_candi <- 1 - pt(t_critical_candi, df_candi, ncp = ncp_candi)
    if(alternative == "two.sided"){
      power_candi <- power_candi + pt(-t_critical_candi, df_candi, ncp = ncp_candi)
    }

    idx <- which(power_candi >= power)[1]
    if (is.na(idx)) {
      warning(paste0("Power did not reach ", power,
                     " within N <= ", max(n_candi), "; the maximal N was returned."))
      idx <- length(n_candi)
    }

    res$n_total   <- n_candi[idx]
    res$df  <- df_candi[idx]
    res$power     <- power_candi[idx]
    res$t_critical <- t_critical_candi[idx]
    res$ncp       <- ncp_candi[idx]

    return(structure(res, class = c("cal_n", "data.frame")))
  }

  ## -------- Solve alpha (given power) --------
  if (is.null(alpha)) {
    res$ncp <- res$delta * sqrt(res$n_total)
    if(design == "two.sample"){
      res$ncp <- res$ncp / 2
    }

    if(alternative == "one.sided"){
      res$t_critical <- qt(1 - res$power, res$df, ncp = res$ncp)
      res$alpha     <- 1 - pt(res$t_critical, res$df)

    } else{
      f <- function(x) {
        (pt(x, res$df, ncp = res$ncp) - pt(-x, res$df, ncp = res$ncp)) - (1 - power)
      }
      upper <- qt(1 - 1e-12, df = res$df)
      root <- uniroot(f, lower = 0, upper = upper, tol = 1e-12, maxiter = 10000)$root
      res$alpha <- 2 * (1 - pt(root, df = res$df))
      res$t_critical <- root
    }
    return(structure(res, class = c("cal_alpha", "data.frame")))
  }

  ## -------- Solve effect size (given N, alpha, power) --------
  if (is.null(delta) && is.null(cohensf) && is.null(peta2)) {
    if(alternative == "one.sided"){
      res$t_critical <- qt(1 - res$alpha, res$df)
    } else{
      res$t_critical <- qt(1 - res$alpha/2, res$df)
    }

    froot <- function(x) 1 - pt(res$t_critical, res$df, ncp = x) - res$power
    upper <- 100
    val_u <- froot(upper)
    while (val_u < 0 && upper < 1e6) {  # increase upper until achievable
      upper <- upper * 2
      val_u <- froot(upper)
    }
    res$ncp <- uniroot(froot, lower = 0, upper = upper, tol = 1e-12, maxiter = 10000)$root

    res$delta <- res$ncp / sqrt(res$n_total)
    if(design == "two.sample"){
      res$delta <- res$delta * 2
    }
    res$cohensf <- res$delta / 2
    res$peta2   <- cohensf_to_peta2(res$cohensf)
    return(structure(res, class = c("cal_es", "data.frame")))
  }
}
