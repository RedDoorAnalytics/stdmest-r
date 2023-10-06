#' @title Read Estimation Results From Stata
#'
#' @param path_eb Path of .xlsx file containing `e(b)` estimation results as
#'     exported by `putexcel` in Stata.
#' @param path_eV Path of .xlsx file containing `e(V)` estimation results as
#'     exported by `putexcel` in Stata.
#'
#' @return A list with two elements, one for `e(b)` (named `eb`) and one for `e(V)`
#'     (named `eV`).
#' @export
#'
read_e <- function(path_eb, path_eV) {
  # Read e(b)
  tbl <- readxl::read_excel(path = path_eb, skip = 1)
  tbl$`...1` <- NULL
  eb <- matrix(data = as.numeric(tbl[1, ]), nrow = 1)
  colnames(eb) <- names(tbl)
  # Read e(V)
  eV <- readxl::read_excel(path = path_eV, col_names = FALSE)
  eV <- as.matrix(eV)
  colnames(eV) <- colnames(eb)
  rownames(eV) <- names(eb)
  # Return
  out <- list(eb = eb, eV = eV)
  return(out)
}
