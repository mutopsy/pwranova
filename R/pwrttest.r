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
#' @param n_total Integer scalar. Total sample size.
#'   If \code{NULL}, the function solves for \code{n_total}.
#' @param alpha Numeric in \eqn{(0,1)}. If \code{NULL}, it is solved for given the other inputs.
#' @param power Numeric in \eqn{(0,1)}. If \code{NULL}, it is computed; if \code{n_total} is \code{NULL},
#'   \code{n_total} is solved to attain this power.
#' @param delta Numeric (non-negative). Cohen's \eqn{d}-type effect size.
#'   If \code{NULL}, it is derived from \code{cohensf} or \code{peta2} when available.
#'   If all three effect-size arguments (\code{delta}, \code{cohensf}, \code{peta2})
#'   are \code{NULL}, then the effect size is treated as the unknown quantity and is
#'   solved for given \code{n_total}, \code{alpha}, and \code{power}.
#'   The exact definition depends on the design:
#'   \itemize{
#'     \item \emph{One-sample}: Cohen's \eqn{d = (\mu - \mu_0)/\sigma}.
#'     \item \emph{Paired}: Cohen's \eqn{d_z = \bar{d}/s_d}, i.e., the mean of the difference scores
#'           divided by their standard deviation.
#'     \item \emph{Two-sample (equal allocation)}: Cohen's \eqn{d} is defined as the mean difference
#'           divided by the pooled standard deviation; internally related to \eqn{f} via \eqn{d = 2f}.
#'   }
#'   If \code{NULL}, \code{delta} is derived from \code{cohensf} or \code{peta2} when available.
#' @param cohensf Numeric (non-negative). Cohen's \eqn{f}.
#'   If \code{NULL}, it can be derived from \code{delta};
#'   if \code{delta} is supplied, \code{cohensf} is ignored.
#'   Effect-size relations by design:
#'   \itemize{
#'     \item Two-sample (equal allocation): \eqn{d = 2f}
#'     \item Paired: \eqn{d_z = f}
#'     \item One-sample: \eqn{f} and \eqn{\eta_p^2} are not supported
#'   }
#' @param peta2 Numeric in \eqn{(0,1)}. Partial eta squared.
#'   If \code{NULL}, it can be derived from \code{cohensf};
#'   if \code{delta} is supplied, \code{peta2} is ignored.
#'   Not defined for one-sample designs.
#' @param alternative Character. Either \code{"two.sided"} or \code{"one.sided"}.
#' @param nlim Integer vector of length 2. Search range of total \code{n} when solving sample size.
#'
#' @details
#' \itemize{
#'   \item If multiple effect-size arguments are supplied (\code{delta}, \code{cohensf}, \code{peta2}),
#'         precedence is \code{delta} \eqn{>} \code{cohensf} \eqn{>} \code{peta2}; the rest are ignored with a warning.
#'   \item For the two-sample design, equal allocation is assumed; \code{n_total} must be even when provided,
#'         and the solved \code{n_total} will be an even number.
#'   \item For the paired design, the effect size is interpreted as \eqn{d_z}.
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
    paired = FALSE, onesample = FALSE,
    n_total = NULL, alpha = NULL, power = NULL, delta = NULL, cohensf = NULL, peta2 = NULL,
    alternative = c("two.sided", "one.sided"),
    nlim = c(2, 10000)
) {
  ## -------- Initial checks & conversions --------

  alternative <- match.arg(alternative, c("two.sided", "one.sided"))

  if (!is.logical(paired) || length(paired) != 1L) {
    stop("'paired' must be a single logical value.")
  }

  if (!is.null(n_total)) {
    if (!is.numeric(n_total) || length(n_total) != 1L || n_total %% 1 != 0) {
      stop("'n_total' must be a single positive integer.")
    }
    if (!paired && !onesample && n_total <= 3) {
      stop("'n_total' must be >= 4 for a two-sample t-test.")
    }
    if ((paired || onesample) && n_total <= 1) {
      stop("'n_total' must be >= 2 for a one-sample or paired t-test.")
    }
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
    warning("Both 'delta' and 'peta2' were supplied; 'peta2' was ignored.")
  }

  if(!is.null(delta)){
    delta <- abs(delta)
    if(design == "two.sample"){
      cohensf <- abs(delta) / 2
      peta2 <- cohensf_to_peta2(cohensf)
    } else if(design == "paired") {
      cohensf <- abs(delta)
      peta2 <- cohensf_to_peta2(cohensf)
    } else{
      cohensf <- NA_real_
      peta2 <- NA_real_
    }
  } else if(!is.null(cohensf)){
    if(design == "two.sample"){
      peta2 <- cohensf_to_peta2(cohensf)
      delta <- cohensf * 2
    } else if(design == "paired"){
      peta2 <- cohensf_to_peta2(cohensf)
      delta <- cohensf
    } else{
      stop("'cohensf' cannot be defined for one-sample t-tests.")
    }

  } else if(!is.null(peta2)){
    if(design == "two.sample"){
      cohensf <- peta2_to_cohensf(peta2)
      delta <- cohensf * 2
    } else if(design == "paired"){
      cohensf <- peta2_to_cohensf(peta2)
      delta <- cohensf
    } else{
      stop("'peta2' cannot be defined for one-sample t-tests.")
    }
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

  if(alternative == "two.sided"){
    divider <- 2
  } else{
    divider <- 1
  }

  if (!is.null(n_total) && design == "two.sample" && (n_total %% 2L != 0)) {
    stop("'n_total' must be even for the two-sample design (equal allocation).")
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

    if(design == "paired") names(res)[names(res)=="delta"] <- "delta_z"
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

    if(design == "paired") names(res)[names(res)=="delta"] <- "delta_z"
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

      # For stable and consistent computation, use pwranova()

      if(design == "two.sample"){
        nlevels_b <- 2
        nlevels_w <- NULL
      } else{
        nlevels_b <- NULL
        nlevels_w <- 2
      }

      if(design == "one.sample"){
        cohensf_tmp <- res$delta
      } else{
        cohensf_tmp <- res$cohensf
      }

      res_anova <- pwranova(
        nlevels_b = nlevels_b, nlevels_w = nlevels_w,
        n_total = n_total, cohensf = cohensf_tmp, power = power)

      res$alpha <- res_anova$alpha
      res$t_critical <- sqrt(res_anova$F_critical)
      # res$t_critical <- qt(1 - res$alpha/2, res$df)

      # f <- function(x) {
      #   (pt(x, res$df, ncp = res$ncp) - pt(-x, res$df, ncp = res$ncp)) - (1 - res$power)
      # }
      #
      # upper <- qt(1 - 1e-12, df = res$df)
      # root <- tryCatch(
      #   uniroot(g, lower = 0, upper = 1 - 1e-12,
      #           tol = 1e-12, maxiter = 10000)$root,
      #   error = function(e) {
      #     stop("Failed to solve alpha: inputs are too extreme (power/N/effect size combination).")
      #   }
      # )
      # res$alpha <- 2 * (1 - pt(root, df = res$df))
      # res$t_critical <- root
    }
    if(design == "paired") names(res)[names(res)=="delta"] <- "delta_z"
    return(structure(res, class = c("cal_alpha", "data.frame")))
  }

  ## -------- Solve effect size (given N, alpha, power) --------
  if (is.null(delta) && is.null(cohensf) && is.null(peta2)) {
    if(alternative == "one.sided"){
      res$t_critical <- qt(1 - res$alpha, res$df)
    } else{
      res$t_critical <- qt(1 - res$alpha/2, res$df)
    }

    if (alternative == "one.sided") {
      froot <- function(x) 1 - pt(res$t_critical, res$df, ncp = x) - res$power
    } else {
      froot <- function(x) {
        1 - (pt(res$t_critical, res$df, ncp = x) -
               pt(-res$t_critical, res$df, ncp = x)) - res$power
      }
    }

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
      res$cohensf <- res$delta / 2
      res$peta2   <- cohensf_to_peta2(res$cohensf)
    } else if(design == "paired"){
      res$cohensf <- res$delta
      res$peta2   <- cohensf_to_peta2(res$cohensf)
    } else{
      res$cohensf <- NA_real_
      res$peta2   <- NA_real_
    }

    if(design == "paired") names(res)[names(res)=="delta"] <- "delta_z"
    return(structure(res, class = c("cal_es", "data.frame")))
  }
}


# やること：fの変換を、pairedの場合はf=deltaとする。1標本の場合は定義できない。ドキュメントもなおす。deltaの定義を説明。
# unirootがたまにエラーを吐く。組み合わせによる？pwranovaもチェック。


