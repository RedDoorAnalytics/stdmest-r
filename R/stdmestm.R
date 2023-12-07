#' @title Partially Marginal Regression Standardisation for Hierarchical Survival Models
#'
#' @description This can be used to compute standardised survival probabilities
#'     (over the case-mix covariates) while fixing the value of one of the random
#'     intercepts while marginalising over the others.
#'     Note that currently only three-levels models are supported, e.g., with a
#'     fixed random intercept and a random intercept that we marginalise over.
#'     For comparison, [stdmest()] fixes all random intercepts instead.
#'
#' @param t A vector or single-column matrix with time points at which standardised
#'     survival (or contrasts thereof) is to be computed.
#' @param X Fixed-effects design matrix, likely produced by [Stata.model.matrix()],
#'     including all fixed effects to standardised over.
#' @param beta Estimated fixed-effects, likely produced by [read_e()].
#' @param Sigma Estimated variance-covariance matrix, likely produced by [read_e()].
#' @param b Value of the random intercept to fix. Must be a single scalar.
#' @param bse Standard error of `'b'`, required for the algorithm used when computing
#'     confidence intervals. Must be a single scalar.
#' @param bref Value of the random intercept to fix as reference level. Defaults
#'     to 0. Required when the argument `contrast = TRUE`. Must be a single scalar.
#' @param brefse Standard error of `'bref'`, required for the algorithm used when computing
#'     confidence intervals. Must be a single scalar.
#' @param varmargname Name/label of the coefficient in 'beta' denoting the variance of
#'     the random intercept to marginalise over (random intercept assumed assumed to
#'     be centered on the value of zero).
#' @param contrast Should a contrast be computed (e.g., survival difference comparing
#'     `'b'` vs `'bref`)? Defaults to `TRUE`.
#' @param distribution Baseline hazard distribution. Possible values are `"exponential"`
#'     and `"weibull"`.
#' @param conf.int Should confidence intervals for any standardised quantity be
#'     computed? Defaults to `FALSE`, as the algorithm is computationally intensive.
#' @param B How many replications should the algorithm for computing confidence
#'     intervals use? Defaults to 100.
#' @param cimethod Method used to obtain confidence intervals. Possible values are
#'     the percentile method (`cimethod = "percentile"`, the default) or the normal
#'     approximation method (`cimethod = "normal"`).
#' @param alpha Confidence level. Defaults to 0.05 for 95% confidence intervals.
#' @param nk Number of nodes used for the numerical integration procedure: larger
#'     values yield more accurate results, at the cost of additional computational
#'     complexity. Defaults to 7, which (anecdotally) yielded a good trade-off
#'     between accuracy and computational cost.
#'
#' @export
stdmestm <- function(t, X, beta, Sigma, b, bse, bref = 0, brefse = 0, varmargname, contrast = FALSE, distribution, conf.int = FALSE, B = 100, cimethod = "percentile", alpha = 0.05, nk = 7) {
  # dt <- read_dta(file = "data-raw/data3Lsim-pp.dta") |>
  #   zap_formats() |>
  #   zap_label() |>
  #   zap_labels() |>
  #   mutate(`0b.X3` = as.numeric(X3 == 0)) |>
  #   mutate(`1.X3` = as.numeric(X3 == 1))
  #
  # ## Load estimation results
  # estimation_results <- read_e(
  #   path_eb = "data-raw/data3Lsim-eb.xlsx",
  #   path_eV = "data-raw/data3Lsim-eV.xlsx"
  # )
  #
  # ## Times for predictions
  # t <- seq(0, max(dt$t), length.out = 5)
  # t <- matrix(t, ncol = 1)
  #
  # ## Usage in our settings:
  # modm <- Stata.model.matrix(
  #   fixed = ~ X1 + X2 + `0b.X3` + `1.X3`,
  #   random = ~ b_hospital + b_provider - 1,
  #   data = dt,
  #   eb = estimation_results$eb
  # )
  #
  # X <- modm$fixed
  # beta <- estimation_results$eb
  # Sigma <- estimation_results$eV
  # b <- 2
  # bse <- 1
  # bref <- 0
  # brefse <- 0
  # varmargname <- "var(_cons[hospital_id>provider_id])"
  # contrast <- TRUE
  # distribution <- "weibull"
  # conf.int <- TRUE
  # B <- 100
  # cimethod <- "percentile"
  # alpha <- 0.05
  # nk <- 7

  # Match
  cimethod <- match.arg(arg = cimethod, choices = c("percentile", "normal"), several.ok = FALSE)
  distribution <- match.arg(arg = distribution, choices = c("exponential", "weibull"), several.ok = FALSE)

  # Check that
  if (!is.finite(beta[, varmargname])) stop("Could not pick 'varmargname' from 'beta'. Check your input values.")
  varmarg <- beta[, varmargname, drop = FALSE]

  # Data for GH quadrature
  GHdata <- fastGHQuad::gaussHermiteData(n = nk)
  #
  t <- matrix(data = t, ncol = 1)
  b <- matrix(data = b, ncol = 1)
  bref <- matrix(data = bref, ncol = 1)
  #
  betaX <- beta[, colnames(beta) %in% colnames(X), drop = FALSE]
  #
  if (distribution == "exponential") {
    ln_p <- matrix(data = 0.0)
  } else if (distribution == "weibull") {
    ln_p <- matrix(data = beta[, "ln_p"])
  }
  if (!is.finite(ln_p)) stop("Could not pick 'ln_p' from 'beta'. Check your input values.", call. = FALSE)

  # Point estimate
  S <- predictMeanIntSurvR(t = t, X = X, betaX = betaX, b = b, ln_p = ln_p, varmarg = varmarg, GHdata = GHdata)
  # S <- predictMeanIntSurv(t = t, X = X, betaX = betaX, b = b, ln_p = ln_p, varmarg = varmarg, GHx = GHdata$x, GHw = GHdata$w)

  # Contrast, if required
  if (contrast) {
    Sref <- predictMeanIntSurvR(t = t, X = X, betaX = betaX, b = bref, ln_p = ln_p, varmarg = varmarg, GHdata = GHdata)
    # Sref <- predictMeanIntSurv(t = t, X = X, betaX = betaX, b = bref, ln_p = ln_p, varmarg = varmarg, GHx = GHdata$x, GHw = GHdata$w)
    Sdiff <- S - Sref
  }
  # CIs
  if (conf.int) {
    # For this, resample new parameters/random effects B times
    new_beta <- mvtnorm::rmvnorm(n = B, mean = beta, sigma = Sigma)
    colnames(new_beta) <- colnames(beta)
    new_b <- mvtnorm::rmvnorm(n = B, mean = b, sigma = diag(bse^2, nrow = length(bse), ncol = length(bse)))
    new_b <- matrixStats::rowSums2(x = new_b)
    new_b <- matrix(data = new_b, ncol = 1)
    new_betaX <- new_beta[, colnames(new_beta) %in% colnames(X), drop = FALSE]
    #
    if (distribution == "exponential") {
      new_ln_p <- matrix(data = 0.0, nrow = B)
    } else if (distribution == "weibull") {
      new_ln_p <- new_beta[, "ln_p", drop = FALSE]
    }
    #
    new_varmarg <- new_beta[, varmargname, drop = FALSE]
    # Make sure new_varmarg are >= 0
    new_varmarg <- pmax(new_varmarg, 0.0)
    #
    if (contrast) {
      new_bref <- mvtnorm::rmvnorm(n = B, mean = bref, sigma = diag(brefse^2, nrow = length(brefse), ncol = length(brefse)))
      new_bref <- matrixStats::rowSums2(x = new_bref)
      new_bref <- matrix(data = new_bref, ncol = 1)
    }
    new_S <- predictMeanIntSurvR(t = t, X = X, betaX = new_betaX, b = new_b, ln_p = new_ln_p, varmarg = new_varmarg, GHdata = GHdata)
    # new_S <- predictMeanIntSurv(t = t, X = X, betaX = new_betaX, b = new_b, ln_p = new_ln_p, varmarg = new_varmarg, GHx = GHdata$x, GHw = GHdata$w)

    if (contrast) {
      new_Sref <- predictMeanIntSurvR(t = t, X = X, betaX = new_betaX, b = new_bref, ln_p = new_ln_p, varmarg = new_varmarg, GHdata = GHdata)
      # new_Sref <- predictMeanIntSurv(t = t, X = X, betaX = new_betaX, b = new_bref, ln_p = new_ln_p, varmarg = new_varmarg, GHx = GHdata$x, GHw = GHdata$w)
      new_Sdiff <- new_S - new_Sref
    }
    if (cimethod == "percentile") {
      S_conf.low <- matrixStats::rowQuantiles(x = new_S, probs = alpha / 2)
      S_conf.high <- matrixStats::rowQuantiles(x = new_S, probs = 1 - alpha / 2)
      if (contrast) {
        Sref_conf.low <- matrixStats::rowQuantiles(x = new_Sref, probs = alpha / 2)
        Sref_conf.high <- matrixStats::rowQuantiles(x = new_Sref, probs = 1 - alpha / 2)
        Sdiff_conf.low <- matrixStats::rowQuantiles(x = new_Sdiff, probs = alpha / 2)
        Sdiff_conf.high <- matrixStats::rowQuantiles(x = new_Sdiff, probs = 1 - alpha / 2)
      }
    } else if (cimethod == "normal") {
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
