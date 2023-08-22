#' @title Create a Model Matrix Compatible With Stata
#'
#' @param ... Everything that is passed to [stats::model.matrix()], e.g., model
#'     formula and data. Ideally, create indicator (one-hot encoding) variables
#'     by hand, thus avoiding calling, e.g., `factor()`.
#' @param eb Fitted model coefficients, i.e., after using [read_e()].
#'
#' @return A model matrix.
#'
#' @export
Stata.model.matrix <- function(..., eb) {
  # Create R-style model matrix
  mm <- stats::model.matrix(...)
  # Adjust intercept name
  colnames(mm)[grepl(pattern = "Intercept", x = colnames(mm))] <- "_cons"
  # Adjust names that are not syntactically valid in R (but okay in Stata)
  idx <- which(x = grepl(pattern = "`", x = colnames(mm)))
  colnames(mm)[idx] <- stringr::str_sub(string = colnames(mm)[idx], start = 2L, end = -2L)
  # Reorder design matrix
  mm <- mm[, names(eb)[names(eb) %in% colnames(mm)]]
  # Return
}
