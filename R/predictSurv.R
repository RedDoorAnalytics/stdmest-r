#' @keywords internal
predictSurv1R <- function(t, X, betaX, b, ln_p = 0.0) {
  S <- exp(-exp(as.numeric(tcrossprod(X, betaX)) + b) * (t^exp(ln_p)))
  # Return survival
  return(S)
}

#' @keywords internal
predictMeanSurv1R <- function(t, X, betaX, b, ln_p) {
  if (nrow(betaX) != nrow(ln_p)) stop("Dimensions of 'betaX' and 'ln_p' do not match.")
  output <- matrix(0, nrow = length(t), ncol = nrow(betaX))
  for (i in 1:nrow(betaX)) {
    output[, i] <- vapply(X = t, FUN = function(tt) mean(predictSurv1R(t = tt, X = X, betaX = betaX[i, , drop = FALSE], b = b[i, ], ln_p = ln_p[i, ])), FUN.VALUE = numeric(length = 1L))
  }
  return(output)
}

#' @keywords internal
predictIntSurvR <- function(t, X, betaX, b, ln_p, varmarg, GHX, GHW, S2) {
  S1 <- predictSurv1R(t = t, X = X, betaX = betaX, b = b + GHX, ln_p = ln_p)
  Sint <- (sqrt(varmarg) * sqrt(2)) * ((S1 * S2) %*% GHW)
  Sint <- as.numeric(Sint)
  return(Sint)
}

#' @keywords internal
predictMeanIntSurvR <- function(t, X, betaX, b, ln_p, varmarg, GHdata) {
  if (nrow(betaX) != nrow(ln_p)) stop("Dimensions of 'betaX' and 'ln_p' do not match.")
  if (nrow(betaX) != nrow(varmarg)) stop("Dimensions of 'betaX' and 'varmarg' do not match.")
  output <- matrix(0, nrow = length(t), ncol = nrow(betaX))
  nn <- nrow(X)
  for (i in 1:nrow(betaX)) {
    GHX <- GHdata$x * sqrt(varmarg[i, ]) * sqrt(2)
    GHW <- GHdata$w * exp(GHdata$x^2)
    S2 <- stats::dnorm(x = GHX, sd = sqrt(varmarg[i, ]))
    # Turn into matrices
    GHX <- matrix(data = GHX, nrow = nn, ncol = length(GHX), byrow = TRUE)
    GHW <- matrix(data = GHW, nrow = length(GHW), ncol = 1)
    S2 <- matrix(data = S2, nrow = nrow(GHX), ncol = ncol(GHX), byrow = TRUE)
    # Do the calculations
    output[, i] <- vapply(X = t, FUN = function(tt) {
      mean(predictIntSurvR(t = tt, X = X, betaX = betaX[i, , drop = FALSE], b = b[i, ], ln_p = ln_p[i, ], varmarg = varmarg[i, ], GHX = GHX, GHW = GHW, S2 = S2))
    }, FUN.VALUE = numeric(length = 1L))
  }
  return(output)
}
