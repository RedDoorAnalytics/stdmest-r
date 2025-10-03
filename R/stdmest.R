#' @title Regression Standardisation for Hierarchical Survival Models
#'
#' @description This can be used to compute standardised survival probabilities
#'     (over the case-mix covariates) while fixing the value of one (or more)
#'     random intercepts.
#'
#' @param t A vector or single-column matrix with time points at which standardised
#'     survival (or contrasts thereof) is to be computed.
#' @param X Fixed-effects design matrix, likely produced by [Stata.model.matrix()],
#'     including all fixed effects to standardised over.
#' @param beta Estimated fixed-effects, likely produced by [read_e()].
#' @param Sigma Estimated variance-covariance matrix, likely produced by [read_e()].
#' @param b Value of the random intercept to fix. Can be a vector if multiple
#'     hierarchical levels are present.
#' @param bse Standard error of `'b'`, required for the algorithm used when computing
#'     confidence intervals. Can be a vector if multiple hierarchical levels are
#'     present.
#' @param bref Value of the random intercept to fix as reference level. Defaults
#'     to 0. Required when the argument `contrast = TRUE`. Can be a vector if multiple
#'     hierarchical levels are present.
#' @param brefse Standard error of `'bref'`, required for the algorithm used when computing
#'     confidence intervals. Can be a vector if multiple hierarchical levels are
#'     present.
#' @param contrast Should a contrast be computed (e.g., survival difference comparing
#'     `'b'` vs `'bref`)? Defaults to `TRUE`.
#' @param distribution Baseline hazard distribution. Possible values are `"exponential"`
#'     and `"weibull"`.
#' @param conf.int Should confidence intervals for any standardised quantity be
#'     computed? Defaults to `FALSE`, as the algorithm is computationally intensive.
#' @param B How many replications should the algorithm for computing confidence
#'     intervals use? Defaults to 1000.
#' @param cimethod Method used to obtain confidence intervals. Possible values are
#'     the percentile method (`cimethod = "percentile"`, the default) or the normal
#'     approximation method (`cimethod = "normal"`).
#' @param cicloglog Should calculations for confidence intervals based on `cimethod = "normal"`
#'     method be performed on the cloglog scale, to avoid confidence intervals outside the
#'     range of possible values? Defaults to `FALSE`. Note that confidence intervals are
#'     symmetric around the point estimates when `cicloglog = FALSE`, and they are not
#'     symmetric otherwise.
#' @param alpha Confidence level. Defaults to 0.05 for 95% confidence intervals.
#'
#' @export
stdmest <- function(t, X, beta, Sigma, b, bse, bref = 0, brefse = 0, contrast = FALSE, distribution, conf.int = FALSE, B = 1000, cimethod = "percentile", cicloglog = FALSE, alpha = 0.05) {
  # library(tidyverse)
  # library(haven)
  # devtools::load_all()
  # dt <- read_dta(file = "data-raw/data3CIA-pp.dta") |>
  #   zap_formats() |>
  #   zap_label() |>
  #   zap_labels() |>
  #   mutate(`0b.mmrc` = as.numeric(mmrc == 0)) |>
  #   mutate(`1.mmrc` = as.numeric(mmrc == 1)) |>
  #   mutate(`2.mmrc` = as.numeric(mmrc == 2)) |>
  #   mutate(`3.mmrc` = as.numeric(mmrc == 3)) |>
  #   mutate(`4.mmrc` = as.numeric(mmrc == 4))
  # estimation_results <- read_e(path = "data-raw/data3CIA-ebV.xlsx")
  # .times <- seq(0, max(dt$months), length.out = 10)
  # modm <- Stata.model.matrix(
  #   fixed = ~ age + fev1pp + `0b.mmrc` + `1.mmrc` + `2.mmrc` + `3.mmrc` + `4.mmrc`,
  #   random = ~ b - 1,
  #   data = dt,
  #   eb = estimation_results$eb
  # )
  # best <- filter(dt, b == min(b)) |>
  #   distinct(b, bse)
  # t <- matrix(.times, ncol = 1)
  # beta <- estimation_results$eb
  # X <- modm$fixed
  # Sigma <- estimation_results$eV
  # b <- best$b
  # bse <- best$bse
  # distribution <- "weibull"
  # contrast <- TRUE
  # conf.int <- TRUE
  # bref <- 0
  # brefse <- 0
  # B <- 10
  # cimethod <- "percentile"
  # cicloglog <- FALSE
  # alpha <- 0.05

  # Match
  cimethod <- match.arg(arg = cimethod, choices = c("percentile", "normal"), several.ok = FALSE)
  distribution <- match.arg(arg = distribution, choices = c("exponential", "weibull"), several.ok = FALSE)

  # Check that length of b, bse, bref, brefse is the same
  if (!(length(b) == length(bse) & length(b) == length(bref) & length(b) == length(brefse))) {
    stop("'b', 'bse', 'bref', and 'brefse' must have the same number of elements.", call. = FALSE)
  }

  #
  t <- matrix(data = t, ncol = 1)

  #
  betaX <- beta[, colnames(beta) %in% colnames(X), drop = FALSE]
  #
  if (distribution == "exponential") {
    ln_p <- matrix(data = 0.0)
  } else if (distribution == "weibull") {
    ln_p <- matrix(data = beta[, "ln_p"])
  }
  if (!is.finite(ln_p)) stop("Could not pick 'ln_p' from 'beta'. Check your input values.", call. = FALSE)

  # Create a single combined b to use
  b_sum <- sum(b)
  b_sum <- matrix(data = b_sum, ncol = 1)

  # Point estimate
  S <- predictMeanSurv1R(t = t, X = X, betaX = betaX, b = b_sum, ln_p = ln_p)
  # Contrast, if required
  if (contrast) {
    bref_sum <- sum(bref)
    bref_sum <- matrix(data = bref_sum, ncol = 1)
    Sref <- predictMeanSurv1R(t = t, X = X, betaX = betaX, b = bref_sum, ln_p = ln_p)
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
    if (contrast) {
      new_bref <- mvtnorm::rmvnorm(n = B, mean = bref, sigma = diag(brefse^2, nrow = length(brefse), ncol = length(brefse)))
      new_bref <- matrixStats::rowSums2(x = new_bref)
      new_bref <- matrix(data = new_bref, ncol = 1)
    }
    new_S <- predictMeanSurv1R(t = t, X = X, betaX = new_betaX, b = new_b, ln_p = new_ln_p)
    if (contrast) {
      new_Sref <- predictMeanSurv1R(t = t, X = X, betaX = new_betaX, b = new_bref, ln_p = new_ln_p)
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
      if (cicloglog) {
        new_S <- cloglog(new_S)
        new_S[new_S == -Inf | new_S == Inf] <- 0.0
        S_conf.low <- cloglog(S) + matrixStats::rowSds(x = new_S) * z.crit
        S_conf.high <- cloglog(S) - matrixStats::rowSds(x = new_S) * z.crit
        S_conf.low <- cexpexp(S_conf.low)
        S_conf.high <- cexpexp(S_conf.high)
      } else {
        S_conf.low <- S - matrixStats::rowSds(x = new_S) * z.crit
        S_conf.high <- S + matrixStats::rowSds(x = new_S) * z.crit
      }
      if (contrast) {
        if (cicloglog) {
          new_Sref <- cloglog(new_Sref)
          new_Sref[new_Sref == -Inf | new_Sref == Inf] <- 0.0
          Sref_conf.low <- cloglog(Sref) + matrixStats::rowSds(x = new_Sref) * z.crit
          Sref_conf.high <- cloglog(Sref) - matrixStats::rowSds(x = new_Sref) * z.crit
          Sref_conf.low <- cexpexp(Sref_conf.low)
          Sref_conf.high <- cexpexp(Sref_conf.high)
        } else {
          Sref_conf.low <- Sref - matrixStats::rowSds(x = new_Sref) * z.crit
          Sref_conf.high <- Sref + matrixStats::rowSds(x = new_Sref) * z.crit
        }
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
