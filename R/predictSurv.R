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
}
