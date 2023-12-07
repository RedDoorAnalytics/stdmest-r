#' @title Read Estimation Results From Stata
#'
#' @param path Path of .xlsx file containing `mestreg` estimation results as
#'     exported by `mestreg_export` in Stata (bundled with the `stdmest`
#'     command).
#'
#' @return A list with two elements, one for `e(b)` (named `eb`) and one for `e(V)`
#'     (named `eV`).
#'
#' @export
#'
read_e <- function(path) {
  # Read e(b)
  tbl <- readxl::read_excel(path = path, skip = 1, sheet = "e(b)")
  tbl$`...1` <- NULL
  eb <- matrix(data = as.numeric(tbl[1, ]), nrow = 1)
  colnames(eb) <- names(tbl)
  # Read e(V)
  eV <- readxl::read_excel(path = path, col_names = FALSE, sheet = "e(V)")
  eV <- as.matrix(eV)
  colnames(eV) <- colnames(eb)
  rownames(eV) <- names(eb)
  # Return
  out <- list(eb = eb, eV = eV)
  return(out)
}
