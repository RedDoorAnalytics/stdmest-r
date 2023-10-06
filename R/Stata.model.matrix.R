#' @title Create a Model Matrix Compatible With Stata
#'
#' @param fixed Formula for the fixed effects part of the model, passed to
#'     [stats::model.matrix()]. Ideally, create indicator (one-hot encoding) variables
#'     by hand, thus avoiding calling, e.g., `factor()`.
#' @param random Formula for the random effects part of the model, passed to
#'     [stats::model.matrix()]. Remove the overall intercept by having `[...] - 1`
#'     in the formula.
#' @param data Dataset in which formulae are to be evaluated.
#' @param eb Fitted model coefficients, i.e., after using [read_e()].
#'
#' @return A model matrix.
#'
#' @export
Stata.model.matrix <- function(fixed, random, data, eb) {
  ## Fixed effects part:
  # Create R-style model matrix
  mm <- stats::model.matrix(fixed, data = data)
  # Adjust intercept name
  colnames(mm)[grepl(pattern = "Intercept", x = colnames(mm))] <- "_cons"
  # Adjust names that are not syntactically valid in R (but okay in Stata)
  idx <- which(x = grepl(pattern = "`", x = colnames(mm)))
  colnames(mm)[idx] <- stringr::str_sub(string = colnames(mm)[idx], start = 2L, end = -2L)
  # Reorder design matrix
  mm <- mm[, colnames(eb)[colnames(eb) %in% colnames(mm)]]
  ## Random effects part:
  # Create R-style model matrix
  rm <- stats::model.matrix(random, data = data)
  # Return
  out <- list(fixed = mm, random = rm)
  return(out)
}
