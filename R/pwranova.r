#' Power Analysis for Between-, Within-, or Mixed-Factor ANOVA
#'
#' Computes power, required total sample size, alpha, or minimal detectable
#' effect size for fixed-effects terms in between-/within-/mixed-factor ANOVA
#' designs.
#'
#' @param nlevels_b Integer vector. Numbers of levels for between-participant factors.
#'   Omit or set \code{NULL} if there is no between-participant factor.
#' @param nlevels_w Integer vector. Numbers of levels for within-participant factors.
#'   Omit or set \code{NULL} if there is no within-participant factor.
#' @param n_total Integer or integer vector. Total sample size(s) across all groups.
#'   If \code{NULL}, the function solves for \code{n_total}.
#' @param cohensf Numeric. Cohen's \eqn{f} (non-negative). If \code{NULL}, it is
#'   derived from \code{peta2} when available.
#' @param peta2 Numeric in \eqn{(0,1)}. Partial eta squared. If \code{NULL}, it is
#'   derived from \code{cohensf} when available.
#' @param alpha Numeric in \eqn{(0,1)}. If \code{NULL}, it is solved for.
#' @param power Numeric in \eqn{(0,1)}. If \code{NULL}, it is computed; if \code{n_total}
#'   is \code{NULL}, \code{n_total} is solved to achieve this power.
#' @param epsilon Numeric in \eqn{(0,1]}. Nonsphericity (Greenhouse–Geisser) parameter
#'   applied to within-participant terms with \eqn{\mathrm{df}_1 \ge 2}. Ignored if no
#'   within-participant factor or if all within factors have 2 levels.
#' @param target Character vector of term labels to compute (e.g., \code{"B1"}, \code{"W1"},
#'   \code{"B1:W1"}, ...). If \code{NULL}, all terms are returned.
#' @param max_nfactor Integer. Safety cap for the total number of factors.
#' @param nlim Integer length-2. Search range of total \code{n} when solving sample size.
#'
#' @details
#' - Numerator df are adjusted by \code{epsilon} for within terms with \eqn{\mathrm{df}_1 \ge 2}.
#' - Denominator df follow standard mixed ANOVA formulas, multiplied by the same
#'   \code{epsilon} for within terms.
#' - Critical values use the central F; power uses the noncentral F with
#'   \eqn{\lambda = f^2 \cdot n_\mathrm{total}}.
#'
#' @return A data frame with S3 class:
#'   \itemize{
#'     \item \code{"cal_power"} when power is calculated given \code{n_total}, \code{alpha}, and effect size;
#'     \item \code{"cal_n"} or \code{"cal_ns"} when \code{n_total} is solved;
#'     \item \code{"cal_alpha"} or \code{"cal_alphas"} when \code{alpha} is solved;
#'     \item \code{"cal_es"} when minimal detectable effect sizes are solved.
#'   }
#'   Columns include \code{term}, \code{df_num}, \code{df_denom}, \code{n_total},
#'   \code{alpha}, \code{power}, \code{cohensf}, \code{peta2}, \code{criticalF}, \code{ncp}, \code{epsilon}.
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.).
#' Hillsdale, NJ: Lawrence Erlbaum Associates.
#'
#' @examples
#' # One between factor (k=3), one within factor (m=4), compute power
#' pwranova(nlevels_b = 3, nlevels_w = 4, n_total = 60,
#'          cohensf = 0.25, alpha = 0.05, power = NULL, epsilon = 0.8)
#'
#' # Solve required total N for target power
#' pwranova(nlevels_b = 2, nlevels_w = NULL, n_total = NULL,
#'          peta2 = 0.06, alpha = 0.05, power = 0.8)
#'
#' @export
pwranova <- function(
    nlevels_b = NULL, nlevels_w = NULL,
    n_total = NULL, cohensf = NULL, peta2 = NULL, alpha = NULL, power = NULL,
    epsilon = 1.0, target = NULL, max_nfactor = 6, nlim = c(2, 10000)
) {
  ## ---------------- Initial checks & conversions ----------------

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

  if (is.null(nlevels_b) && is.null(nlevels_w)) {
    stop("Either 'nlevels_b' or 'nlevels_w' must be specified.")
  }

  if (!is.null(nlevels_b)) {
    if (!is.numeric(nlevels_b) || any(nlevels_b %% 1 != 0)) {
      stop("Elements of 'nlevels_b' must be integers.")
    }
    if (any(nlevels_b < 2)) stop("Elements of 'nlevels_b' must be 2 or more.")
  }
  if (!is.null(nlevels_w)) {
    if (!is.numeric(nlevels_w) || any(nlevels_w %% 1 != 0)) {
      stop("Elements of 'nlevels_w' must be integers.")
    }
    if (any(nlevels_w < 2)) stop("Elements of 'nlevels_w' must be 2 or more.")
  }

  if (length(nlim) != 2 || any(nlim %% 1 != 0)) stop("'nlim' must be two integers.")
  if (nlim[1] < 2) stop("'nlim[1]' must be 2 or larger.")
  if (nlim[1] >= nlim[2]) stop("'nlim[1]' must be smaller than 'nlim[2]'.")

  if (!is.null(alpha)) {
    if (length(alpha) != 1 && !is.null(target) && length(alpha) != length(target)) {
      stop("'alpha' must be length 1 or match the number of target terms.")
    }
    if (any(alpha <= 0 | alpha >= 1)) stop("'alpha' must be in (0, 1).")
  }
  if (!is.null(power)) {
    if (length(power) != 1 && !is.null(target) && length(power) != length(target)) {
      stop("'power' must be length 1 or match the number of target terms.")
    }
    if (any(power <= 0 | power >= 1)) stop("'power' must be in (0, 1).")
  }

  if (epsilon > 1 || epsilon <= 0) {
    stop("Nonsphericity parameter 'epsilon' must be in (0, 1].")
  }

  if ((is.null(n_total) + is.null(cohensf) + is.null(alpha) + is.null(power)) != 1) {
    stop("Exactly one of 'n_total', 'cohensf' (or 'peta2'), 'alpha', or 'power' must be NULL.")
  }

  ## ---------------- Common preprocessing ----------------

  nlevels <- c(if (is.null(nlevels_b)) integer(0) else nlevels_b,
               if (is.null(nlevels_w)) integer(0) else nlevels_w)

  nfactors_b <- length(if (is.null(nlevels_b)) integer(0) else nlevels_b)
  nfactors_w <- length(if (is.null(nlevels_w)) integer(0) else nlevels_w)
  nfactors   <- nfactors_b + nfactors_w

  if (nfactors > max_nfactor) {
    stop(paste0("The number of factors (", nfactors, ") exceeds 'max_nfactor'. ",
                "If you still wish to proceed, increase 'max_nfactor' manually."))
  }

  nfixef <- 2^nfactors - 1

  ngroups <- if (is.null(nlevels_b)) 1L else prod(nlevels_b)

  if (!is.null(n_total)) {
    if (any(n_total %% ngroups != 0)) {
      stop(paste0("'n_total' must be a multiple of the number of groups = ", ngroups, "."))
    }
  }

  # Factor labels (main effects)
  factor_lab_b <- if (nfactors_b == 0) character(0) else paste0("B", seq_len(nfactors_b))
  factor_lab_w <- if (nfactors_w == 0) character(0) else paste0("W", seq_len(nfactors_w))
  factor_lab   <- c(factor_lab_b, factor_lab_w)

  # Design matrix of terms (exclude intercept)
  designmat <- replicate(nfactors, c(0L, 1L), simplify = FALSE)
  designmat <- expand.grid(designmat)
  designmat <- as.matrix(designmat)
  designmat <- designmat[-1, , drop = FALSE] # drop intercept row

  # Order by interaction order (0=main, 1=2-way, ...)
  interaction_order <- rowSums(designmat) - 1L
  ord <- order(interaction_order)
  designmat <- designmat[ord, , drop = FALSE]
  interaction_order <- interaction_order[ord]

  # Terms including any within factor?
  if (nfactors_b == 0) {
    is_inc_within <- rep(1L, nfixef)
  } else {
    tmp <- designmat
    if (nfactors_b > 0) tmp[, seq_len(nfactors_b)] <- 0L
    is_inc_within <- as.integer(rowSums(tmp) > 0)
  }

  # Term labels
  term_lab <- character(nfixef)
  for (i in seq_len(nfixef)) {
    on <- which(designmat[i, ] == 1L)
    term_lab[i] <- paste0(factor_lab[on], collapse = ":")
  }

  # Numerator df (main effects) and for interactions
  df_main_num <- nlevels - 1L
  dfmat <- designmat
  for (i in seq_len(nfactors)) {
    dfmat[, i] <- designmat[, i] * df_main_num[i]
  }
  dfmat[dfmat == 0] <- 1L
  df_num <- apply(dfmat, 1L, prod)

  # Numerator df only from Within-factor
  dfmat_onlywithin <- dfmat
  if(nfactors_b == 0){
    df_num_onlywithin <- df_num
  } else if(nfactors_w >= 1){
    dfmat_within <- dfmat
    dfmat_within[,1:nfactors_b] <- 1
    df_num_onlywithin <- apply(dfmat_within, 1L, prod)
  } else{
    df_num_onlywithin <- rep(1, length(df_num))
  }

  # Result scaffold
  res <- data.frame(
    term      = term_lab,
    df_num    = df_num,
    df_denom  = NA_real_,
    n_total   = NA_real_,
    alpha     = NA_real_,
    power     = NA_real_,
    cohensf   = if (is.null(cohensf)) NA_real_ else cohensf,
    peta2     = if (is.null(peta2))   NA_real_ else peta2,
    criticalF = NA_real_,
    ncp       = NA_real_,
    epsilon   = 1
  )

  # Epsilon assignment (within terms with df1>=2)
  if (any(is_inc_within == 1L)) {
    res$epsilon[is_inc_within == 1L] <- epsilon
    res$epsilon[is_inc_within == 1L & res$df_num == 1] <- 1
  }

  # Adjust numerator df by epsilon
  res$df_num <- res$df_num * res$epsilon

  # Prepare known columns
  if (!is.null(n_total)) res$n_total <- n_total
  if (!is.null(alpha))   res$alpha   <- alpha
  if (!is.null(power))   res$power   <- power

  # Filter by target
  if (!is.null(target)) {
    if (!any(res$term %in% target)) stop("'target' contains no valid term label.")
    is_inc_within <- is_inc_within[res$term %in% target]
    df_num_onlywithin <- df_num_onlywithin[res$term %in% target]
    res <- res[res$term %in% target, , drop = FALSE]
  }
  nrow_res <- nrow(res)

  ## ---------------- Power (given N) ----------------
  if (is.null(power)) {
    df_denom_base <- n_total - ngroups
    tmp <- is_inc_within * df_num_onlywithin
    tmp[tmp == 0] <- 1
    if (length(res$n_total) == 1L) {
      res$df_denom <- rep(df_denom_base, nrow_res) * tmp * res$epsilon
    } else {
      res$df_denom <- df_denom_base * tmp * res$epsilon
    }
    res$criticalF <- qf(1 - res$alpha, res$df_num, res$df_denom)
    res$ncp <- res$cohensf^2 * res$n_total
    res$power <- 1 - pf(res$criticalF, res$df_num, res$df_denom, ncp = res$ncp)
    return(structure(res, class = c("cal_power", "data.frame")))
  }

  ## ---------------- Solve N (given target power) ----------------
  if (is.null(n_total)) {
    for (i in seq_len(nrow_res)) {
      tgt_power <- if (length(power) == 1L) power else power[i]
      nmin <- ceiling(nlim[1] / ngroups) * ngroups
      if (nmin <= ngroups) nmin <- ngroups * 2L
      n_candi <- seq.int(nmin, nlim[2], by = ngroups)

      df_denom_base <- n_candi - ngroups
      tmp <- is_inc_within[i] * df_num_onlywithin[i]
      if (tmp == 0) tmp <- 1
      df_denom_candi <- df_denom_base * tmp * res$epsilon[i]

      criticalF_candi <- qf(1 - res$alpha[i], res$df_num[i], df_denom_candi)
      ncp_candi <- res$cohensf[i]^2 * n_candi
      power_candi <- 1 - pf(criticalF_candi, res$df_num[i], df_denom_candi, ncp = ncp_candi)

      idx <- which(power_candi >= tgt_power)[1]
      if (is.na(idx)) {
        warning(paste0("Power for ", res$term[i], " did not reach ", tgt_power,
                       " within N <= ", max(n_candi), "; the maximal N was returned."))
        idx <- length(n_candi)
      }

      res$n_total[i]   <- n_candi[idx]
      res$df_denom[i]  <- df_denom_candi[idx]
      res$power[i]     <- power_candi[idx]
      res$criticalF[i] <- criticalF_candi[idx]
      res$ncp[i]       <- ncp_candi[idx]
    }

    cls <- if (nrow_res == 1L) "cal_n" else "cal_ns"
    return(structure(res, class = c(cls, "data.frame")))
  }

  ## ---------------- Solve alpha (given power) ----------------
  if (is.null(alpha)) {
    df_denom_base <- n_total - ngroups
    tmp <- is_inc_within * df_num_onlywithin
    tmp[tmp == 0] <- 1
    res$df_denom <- rep(df_denom_base, nrow_res) * tmp * res$epsilon
    res$ncp <- res$cohensf^2 * res$n_total
    # critical F such that power is achieved
    res$criticalF <- qf(1 - res$power, res$df_num, res$df_denom, ncp = res$ncp)
    res$alpha <- 1 - pf(res$criticalF, res$df_num, res$df_denom)
    cls <- if (nrow_res == 1L) "cal_alpha" else "cal_alphas"
    return(structure(res, class = c(cls, "data.frame")))
  }

  ## ---------------- Solve effect size (given N, alpha, power) ----------------
  if (is.null(cohensf) && is.null(peta2)) {
    df_denom_base <- n_total - ngroups
    tmp <- is_inc_within * df_num_onlywithin
    tmp[tmp == 0] <- 1
    res$df_denom <- rep(df_denom_base, nrow_res) * tmp * res$epsilon
    res$criticalF <- qf(1 - res$alpha, res$df_num, res$df_denom)

    # ncp root solve with adaptive bracketing
    for (i in seq_len(nrow_res)) {
      froot <- function(x) 1 - pf(res$criticalF[i], res$df_num[i], res$df_denom[i], ncp = x) - res$power[i]
      upper <- 100
      val_u <- froot(upper)
      while (val_u < 0 && upper < 1e6) {  # increase until power >= target is attainable
        upper <- upper * 2
        val_u <- froot(upper)
      }
      res$ncp[i] <- uniroot(froot, lower = 0, upper = upper)$root
    }

    res$cohensf <- sqrt(res$ncp / res$n_total)
    res$peta2   <- cohensf_to_peta2(res$cohensf)
    return(structure(res, class = c("cal_es", "data.frame")))
  }
}
