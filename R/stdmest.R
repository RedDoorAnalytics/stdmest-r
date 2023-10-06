#' @title Regression Standardisation for Two-Levels Hierarchical Survival Models
#'
#' @description A short description needs to go here...
#'
#' @param t A vector or single-column matrix with time points at which standardised
#'     survival (or contrasts thereof) is to be computed.
#' @param X Fixed-effects design matrix, likely produced by [Stata.model.matrix()],
#'     including all fixed effects to standardised over.
#' @param beta Estimated fixed-effects, likely produced by [read_e()].
#' @param Sigma Estimated variance-covariance matrix, likely produced by [read_e()].
#' @param b Value of the random intercept to fix.
#' @param bse Standard error of `'b'`, required for the algorithm used when computing
#'     confidence intervals.
#' @param bref Value of the random intercept to fix as reference level. Defaults
#'     to 0. Required when the argument `contrast = TRUE`.
#' @param brefse Standard error of `'bref'`, required for the algorithm used when computing
#'     confidence intervals.
#' @param contrast Should a contrast be computed (e.g., survival difference comparing
#'     `'b'` vs `'bref`)? Defaults to `TRUE`.
#' @param distribution Baseline hazard distribution. Possible values are `"exponential"`
#'     and `"weibull"`.
#' @param conf.int Should confidence intervals for any standardised quantity be
#'     computed? Defaults to `FALSE`, as the algorithm is computationally intensive.
#' @param B How many replications should the algorithm for computing confidence
#'     intervals use? Defaults to 100.
#' @param method Method used to obtain confidence intervals. Possible values are
#'     the percentile method (`method = "percentile"`) or the normal approximation
#'     method (`method = "normal"`).
#' @param alpha Confidence level. Defaults to 0.05 for 95% confidence intervals.
#'
#' @export
stdmest <- function(t, X, beta, Sigma, b, bse, bref = 0, brefse = 0, contrast = FALSE, distribution, conf.int = FALSE, B = 100, method = "percentile", alpha = 0.05) {
  # t <- matrix(.times, ncol = 1)
  # X <- modm$fixed
  # beta <- estimation_results$eb
  # Sigma <- estimation_results$eV
  # b <- 1
  # bse <- 1
  # bref <- 0
  # brefse <- 0
  # contrast <- TRUE
  # distribution <- "weibull"
  # B <- 50
  # method <- "percentile"
  # alpha <- 0.05
  # conf.int <- TRUE

  # Match
  method <- match.arg(arg = method, choices = c("percentile", "normal"), several.ok = FALSE)
  distribution <- match.arg(arg = distribution, choices = c("exponential", "weibull"), several.ok = FALSE)

  #
  t = matrix(data = t, ncol = 1)

  #
  betaX <- beta[, colnames(beta) %in% colnames(X), drop = FALSE]
  #
  if (distribution == "exponential") {
    ln_p <- matrix(data = 0.0)
  } else if (distribution == "weibull") {
    ln_p <- matrix(data = beta[, "ln_p"])
  }
  if (!is.finite(ln_p)) stop("Something went wrong, could not extract 'ln_p' from 'beta'. Check your input values.", call. = FALSE)

  # Point estimate
  S <- predictMeanSurv1R(t = t, X = X, betaX = betaX, b = b, ln_p = ln_p)
  # Contrast, if required
  if (contrast) {
    Sref <- predictMeanSurv1R(t = t, X = X, betaX = betaX, b = bref, ln_p = ln_p)
    Sdiff <- S - Sref
  }
  # CIs
  if (conf.int) {
    # For this, resample new parameters/random effects B times
    new_beta <- mvtnorm::rmvnorm(n = B, mean = beta, sigma = Sigma)
    colnames(new_beta) <- colnames(beta)
    new_b <- matrix(data = stats::rnorm(n = B, mean = b, sd = bse), ncol = 1)
    new_betaX <- new_beta[, colnames(new_beta) %in% colnames(X), drop = FALSE]
    #
    if (distribution == "exponential") {
      new_ln_p <- matrix(data = 0.0, nrow = B)
    } else if (distribution == "weibull") {
      new_ln_p <- new_beta[, "ln_p", drop = FALSE]
    }
    if (contrast) {
      new_bref <- matrix(data = stats::rnorm(n = B, mean = bref, sd = brefse), ncol = 1)
    }
    new_S <- predictMeanSurv1(t = t, X = X, betaX = new_betaX, b = new_b, ln_p = new_ln_p)
    if (contrast) {
      new_Sref <- predictMeanSurv1(t = t, X = X, betaX = new_betaX, b = new_bref, ln_p = new_ln_p)
      new_Sdiff <- new_S - new_Sref
    }
    if (method == "percentile") {
      S_conf.low <- matrixStats::rowQuantiles(x = new_S, probs = alpha / 2)
      S_conf.high <- matrixStats::rowQuantiles(x = new_S, probs = 1 - alpha / 2)
      if (contrast) {
        Sref_conf.low <- matrixStats::rowQuantiles(x = new_Sref, probs = alpha / 2)
        Sref_conf.high <- matrixStats::rowQuantiles(x = new_Sref, probs = 1 - alpha / 2)
        Sdiff_conf.low <- matrixStats::rowQuantiles(x = new_Sdiff, probs = alpha / 2)
        Sdiff_conf.high <- matrixStats::rowQuantiles(x = new_Sdiff, probs = 1 - alpha / 2)
      }
    } else if (method == "normal") {
      z.crit <- stats::qnorm(p = 1 - alpha / 2)
      S_conf.low <- S - matrixStats::rowSds(x = new_S) * z.crit
      S_conf.high <- S + matrixStats::rowSds(x = new_S) * z.crit
      if (contrast) {
        Sref_conf.low <- Sref - matrixStats::rowSds(x = new_Sref) * z.crit
        Sref_conf.high <- Sref + matrixStats::rowSds(x = new_Sref) * z.crit
        Sdiff_conf.low <- Sdiff - matrixStats::rowSds(x = new_Sdiff) * z.crit
        Sdiff_conf.high <- Sdiff + matrixStats::rowSds(x = new_Sdiff) * z.crit
      }
    }
    #
    out <- data.frame(
      t = t,
      S = S,
      S_conf.low = S_conf.low,
      S_conf.high = S_conf.high
    )
    if (contrast) {
      out <- cbind.data.frame(
        out,
        data.frame(
          Sref = Sref,
          Sref_conf.low = Sref_conf.low,
          Sref_conf.high = Sref_conf.high,
          Sdiff = Sdiff,
          Sdiff_conf.low = Sdiff_conf.low,
          Sdiff_conf.high = Sdiff_conf.high
        )
      )
    }
  } else {
    out <- data.frame(
      t = t,
      S = S
    )
    if (contrast) {
      out <- cbind.data.frame(
        out,
        data.frame(
          Sref = Sref,
          Sdiff = Sdiff
        )
      )
    }
  }
  return(out)
}
