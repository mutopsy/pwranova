#' Power Analysis for Pearson Correlation
#'
#' Computes statistical \strong{power}, required total sample size, \eqn{\alpha},
#' or the minimal detectable correlation coefficient for a Pearson correlation test.
#' Two computational methods are supported: exact noncentral \emph{t} (\code{method = "t"})
#' and Fisher's \emph{z} transformation with normal approximation (\code{method = "z"}).
#'
#' @param alternative Character. Either \code{"two.sided"} or \code{"one.sided"}.
#' @param n_total Integer scalar. Total sample size (\eqn{n}).
#'   If \code{NULL}, the function solves for \code{n_total}.
#'   Must be \eqn{\ge 3} for \code{method = "t"} and \eqn{\ge 4} for \code{method = "z"}.
#' @param alpha Numeric in \eqn{(0,1)}. If \code{NULL}, it is solved for.
#' @param power Numeric in \eqn{(0,1)}. If \code{NULL}, it is computed; if \code{n_total} is \code{NULL},
#'   \code{n_total} is solved to attain this power.
#' @param rho Numeric correlation coefficient in \eqn{(-1,1)}, nonzero.
#'   If \code{NULL}, \code{rho} is solved for given the other inputs.
#' @param method Character. Either \code{"t"} (noncentral-\emph{t} distribution)
#'   or \code{"z"} (Fisher's \emph{z} transformation with normal approximation).
#' @param bias_correction Logical. Applies only to \code{method = "z"}.
#'   If \code{TRUE}, uses the bias-corrected Fisher \emph{z} transformation
#'   \eqn{z_p = \operatorname{atanh}(r) + r/(2(n-1))}.
#' @param nlim Integer vector of length 2. Search range of \code{n_total} when solving sample size.
#'
#' @details
#' - Exactly one of \code{n_total}, \code{rho}, \code{alpha}, or \code{power} must be \code{NULL}.
#' - For \code{method = "t"}, computations are based on the noncentral \emph{t}
#'   distribution with noncentrality parameter
#'   \eqn{\lambda = \tfrac{\rho}{\sqrt{1-\rho^2}} \sqrt{n}}.
#' - For \code{method = "z"}, computations use Fisher's \emph{z} transformation
#'   with variance 1 under the alternative hypothesis. The returned \code{ncp}
#'   corresponds to the mean of this normal distribution,
#'   \eqn{\mu = \sqrt{n-3}\, z_p}, where \eqn{z_p} is the (possibly bias-corrected)
#'   transformed correlation.
#' - Note: results from \code{method = "z"} will not exactly match \code{pwr::pwr.r.test},
#'   because \code{pwr} uses a hybrid approach combining the Fisher-\emph{z} approximation
#'   with a \emph{t}-based critical value.
#'
#' @return A one-row \code{data.frame} with class
#'   \code{"cal_power"}, \code{"cal_n"}, \code{"cal_alpha"}, or \code{"cal_es"},
#'   depending on the solved quantity. Columns:
#'   \itemize{
#'     \item \code{df} (only for \code{method = "t"})
#'     \item \code{n_total}, \code{alpha}, \code{power}
#'     \item \code{rho}, \code{t_critical} or \code{z_critical}
#'     \item \code{ncp} (noncentrality parameter or mean under the alternative:
#'           see Details)
#'   }
#'
#' @examples
#' # (1) Compute power for rho = 0.3, N = 50, two-sided test
#' pwrcortest(alternative = "two.sided", n_total = 50, rho = 0.3, alpha = 0.05)
#'
#' # (2) Solve required N for target power, using Fisher-z method
#' pwrcortest(method = "z", rho = 0.2, alpha = 0.05, power = 0.8)
#'
#' # (3) Solve minimal detectable correlation
#' pwrcortest(n_total = 60, alpha = 0.05, power = 0.9, rho = NULL)
#'
#' @importFrom stats qt pt qnorm pnorm uniroot
#' @export

pwrcortest <- function(
   alternative = c("two.sided", "one.sided"),
   n_total = NULL, alpha = NULL, power = NULL, rho = NULL,
   method = c("t", "z"), bias_correction = FALSE,
   nlim = c(2, 10000)
) {
  ## -------- Initial checks & conversions --------

  alternative <- match.arg(alternative, c("two.sided", "one.sided"))
  method <- match.arg(method, c("t", "z"))

  df <- NA_real_

  if (!is.null(n_total)) {
    if (!is.numeric(n_total) || length(n_total) != 1L || n_total %% 1 != 0) {
      stop("'n_total' must be a single positive integer.")
    }

    if (n_total < 3) {
      stop("'n_total' must be >= 3 for a Pearson correlation test.")
    }

    if (n_total < 4 && method == "z") {
      stop("'n_total' must be >= 4 when method = 'z'.")
    }

    df <- n_total - 2
  }

  if(bias_correction && method == "t"){
    warning("'bias_correction' applies only when method = 'z'; argument was ignored.")
    bias_correction <- FALSE
  }

  if (length(nlim) != 2L || any(nlim %% 1 != 0)) stop("'nlim' must be two integers.")
  if (nlim[1] < 2) stop("'nlim[1]' must be 2 or larger.")
  if (nlim[1] >= nlim[2]) stop("'nlim[1]' must be smaller than 'nlim[2]'.")

  if(!is.null(rho)){
    if (length(rho) != 1L) stop("'rho' must be length 1.")
    if (rho <= -1 || rho >= 1) stop("'rho' must be in (-1, 1).")
    if (rho == 0) stop("'rho' must be nonzero.")
    rho <- abs(rho)
  }

  if (!is.null(alpha)) {
    if (length(alpha) != 1L) stop("'alpha' must be length 1.")
    if (alpha <= 0 || alpha >= 1) stop("'alpha' must be in (0, 1).")
  }
  if (!is.null(power)) {
    if (length(power) != 1L) stop("'power' must be length 1.")
    if (power <= 0 || power >= 1) stop("'power' must be in (0, 1).")
  }

  if ((is.null(n_total) + is.null(rho) + is.null(alpha) + is.null(power)) != 1) {
    stop("Exactly one of 'n_total', 'rho', 'alpha', or 'power' must be NULL.")
  }

  if(alternative == "two.sided"){
    divider <- 2
  } else{
    divider <- 1
  }

  ## -------- Result scaffold --------
  res <- data.frame(
    df        = df,
    n_total   = NA_real_,
    alpha     = NA_real_,
    power     = NA_real_,
    rho       = if (is.null(rho)) NA_real_ else rho,
    t_critical = NA_real_,
    ncp       = NA_real_
  )

  if (!is.null(n_total)) res$n_total <- n_total
  if (!is.null(alpha))   res$alpha   <- alpha
  if (!is.null(power))   res$power   <- power

  if(method == "z"){
    colnames(res)[colnames(res)=="t_critical"] <- "z_critical"
    res <- subset(res, select = -df)
    }

  ## -------- Power (given N) --------
  if (is.null(power)) {

    if(method == "t"){
      res$t_critical <- qt(1 - res$alpha / divider, res$df)
      res$ncp <- (res$rho / sqrt(1 - res$rho^2)) * sqrt(res$n_total)
      res$power <- 1 - pt(res$t_critical, res$df, ncp = res$ncp)
      if(alternative == "two.sided"){
        res$power <- res$power + pt(-res$t_critical, res$df, ncp = res$ncp)
      }
    } else if(method == "z"){

      if(bias_correction){
        zp <- atanh(res$rho) + res$rho / (2 * (res$n_total - 1))
      } else{
        zp <- atanh(res$rho)
      }

      mu <- sqrt(res$n_total - 3) * zp

      res$z_critical <- qnorm(1 - res$alpha / divider, mean = 0, sd = 1)
      res$ncp <- mu
      res$power <- 1 - pnorm(res$z_critical,  mean = mu, sd = 1)
      if(alternative == "two.sided"){
        res$power <- res$power + pnorm(-res$z_critical, mean = mu, sd = 1)
      }
    }

    return(structure(res, class = c("cal_power", "data.frame")))
  }

  ## -------- Solve N (given target power) --------
  if (is.null(n_total)) {
    if(method == "t"){
      nmin <- max(ceiling(nlim[1]),3)
      n_candi <- seq.int(nmin, nlim[2], by = 1L)
      df_candi <- n_candi - 2
      t_critical_candi <- qt(1 - res$alpha / divider, df_candi)
      ncp_candi <- (res$rho / sqrt(1 - res$rho^2)) * sqrt(n_candi)
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

    } else if(method == "z"){

      nmin <- max(ceiling(nlim[1]),4)
      n_candi <- seq.int(nmin, nlim[2], by = 1L)

      if(bias_correction){
        zp_candi <- atanh(res$rho) + res$rho / (2 * (n_candi - 1))
      } else{
        zp_candi <- atanh(res$rho)
      }

      mu_candi <- sqrt(n_candi - 3) * zp_candi
      z_critical <- qnorm(1 - res$alpha / divider, mean = 0, sd = 1)

      power_candi <- 1 - pnorm(z_critical, mean = mu_candi, sd = 1)
      if(alternative == "two.sided"){
        power_candi <- power_candi + pnorm(-z_critical, mean = mu_candi, sd = 1)
      }

      idx <- which(power_candi >= power)[1]
      if (is.na(idx)) {
        warning(paste0("Power did not reach ", power,
                       " within N <= ", max(n_candi), "; the maximal N was returned."))
        idx <- length(n_candi)
      }

      res$n_total   <- n_candi[idx]
      res$power     <- power_candi[idx]
      res$z_critical <- z_critical
      res$ncp       <- mu_candi[idx]
    }

    return(structure(res, class = c("cal_n", "data.frame")))
  }

  ## -------- Solve alpha (given power) --------
  if (is.null(alpha)) {
    if(method == "t"){
      res$ncp <- (res$rho / sqrt(1 - res$rho^2)) * sqrt(res$n_total)

      if(alternative == "one.sided"){
        res$t_critical <- qt(1 - res$power, res$df, ncp = res$ncp)
        res$alpha     <- 1 - pt(res$t_critical, res$df)

      } else{

        # For stable and consistent computation, use pwranova()

        res_anova <- pwranova(
          nlevels_b = 2,
          n_total = n_total, cohensf = (res$rho / sqrt(1 - res$rho^2)), power = power)

        res$alpha <- res_anova$alpha
        res$t_critical <- sqrt(res_anova$F_critical)
      }

    } else if(method == "z"){

      if(bias_correction){
        zp <- atanh(res$rho) + res$rho / (2 * (res$n_total - 1))
      } else{
        zp <- atanh(res$rho)
      }

      mu <- sqrt(res$n_total - 3) * zp
      res$ncp <- mu

      if(alternative == "one.sided"){
        res$z_critical <- qnorm(1 - res$power, mean = mu, sd = 1)
        res$alpha     <- 1 - pnorm(res$z_critical)

      } else{

        f_z <- function(x){
          pnorm(x, mean = mu, sd = 1) - pnorm(-x, mean = mu, sd = 1) - (1 - res$power)
        }

        upper <- 100
        val_u <- f_z(upper)
        while (val_u < 0 && upper < 1e6) {  # increase upper until achievable
          upper <- upper * 2
          val_u <- f_z(upper)
        }

        root <- uniroot(f_z, lower = 0, upper = upper, tol = 1e-12, maxiter = 10000)$root

        res$z_critical <- root
        res$alpha <- (1 - pnorm(root)) * 2

      }

    }

    return(structure(res, class = c("cal_alpha", "data.frame")))
  }

  ## -------- Solve effect size (given N, alpha, power) --------
  if (is.null(rho)) {
    if(method == "t"){
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

      tmp <- res$ncp / sqrt(res$n_total)
      res$rho <- tmp / sqrt(1 + tmp^2)

    } else if(method == "z"){

      if(alternative == "one.sided"){
        res$z_critical <- qnorm(1 - res$alpha)
      } else{
        res$z_critical <- qnorm(1 - res$alpha/2)
      }

      if (alternative == "one.sided") {
        froot <- function(x) 1 - pnorm(res$z_critical, x, 1) - res$power
      } else {
        froot <- function(x) {
          1 - (pnorm(res$z_critical, x, 1) -
               pnorm(-res$z_critical, x, 1)) - res$power
        }
      }

      upper <- 100
      val_u <- froot(upper)
      while (val_u < 0 && upper < 1e6) {  # increase upper until achievable
        upper <- upper * 2
        val_u <- froot(upper)
      }
      mu <- uniroot(froot, lower = 0, upper = upper, tol = 1e-12, maxiter = 10000)$root
      res$ncp <- mu
      zp <- mu / sqrt(res$n_total - 3)

      if(bias_correction){
        f <- function(rho) atanh(rho) + rho/(2*(res$n_total-1)) - zp
        res$rho <- uniroot(f, interval = c(-0.9999, 0.9999))$root

      } else{
        res$rho <- tanh(zp)
      }
    }
    return(structure(res, class = c("cal_es", "data.frame")))
  }
}

# やること: rho_to_deltaつくる？
