#' @title Read Estimation Results From Stata
#'
#' @param path Path of .xlsx file containing `mestreg` or `stmixed` estimation
#'     results as exported by `modexpt` in Stata (bundled with the `stdmest`
#'     command).
#'
#' @return A list with different components (depending on whether we are reading
#'     in results from `mestreg` or `stmixed`. Common elements will be `eb` (for
#'     `e(b)`) and `eV` (for `e(V)`), for the estimated coefficients and their
#'     variance-covariance matrix, among others. Check the input .xlsx file and
#'     `ereturn list` in Stata (after `mestreg`, `stmixed`) for more details.
#'
#' @export
#'
read_e <- function(path) {
  # path <- "data-raw/data3CIA-ebV-weibull.xlsx"
  # path <- "data-raw/data3CIA-ebV-rp3-orthog.xlsx"
  # path <- "data-raw/data3CIA-ebV-rp3-noorthog.xlsx"

  # Identify model
  cmd <- readxl::read_excel(path = path, col_names = FALSE, sheet = "cmd")
  cmd <- unlist(cmd, use.names = FALSE)
  # Read e(b)
  r1 <- readxl::read_excel(path = path, sheet = "e(b)", col_names = FALSE, n_max = 1)
  r1 <- unlist(r1, use.names = FALSE)
  r2 <- readxl::read_excel(path = path, sheet = "e(b)", col_names = FALSE, skip = 1, n_max = 1)
  r2 <- unlist(r2, use.names = FALSE)
  r3 <- readxl::read_excel(path = path, sheet = "e(b)", col_names = FALSE, skip = 2, n_max = 1)
  r3$`...1` <- NULL
  eb <- matrix(data = as.numeric(r3[1, ]), nrow = 1)
  # Fix names if -stmixed- models
  if (cmd == "stmixed") {
    namesfix <- readxl::read_excel(path = path, col_names = FALSE, sheet = "e(cmplabels1)")
    namesfix <- unlist(namesfix)
    namesfix <- stringr::str_split_1(string = namesfix, pattern = " ")
    r2[1:length(namesfix)] <- namesfix
  }
  # Assign names
  colnames(eb) <- paste0(r1, ":", r2)
  # Read e(V)
  eV <- readxl::read_excel(path = path, col_names = FALSE, sheet = "e(V)")
  eV <- as.matrix(eV)
  colnames(eV) <- colnames(eb)
  rownames(eV) <- colnames(eb)
  # Read family
  family <- readxl::read_excel(path = path, col_names = FALSE, sheet = "family")
  family <- unlist(family, use.names = FALSE)
  # For -mestreg- models, read parametrisation
  if (cmd == "mestreg") {
    frm <- readxl::read_excel(path = path, col_names = FALSE, sheet = "frm")
    frm <- unlist(frm, use.names = FALSE)
  }
  # Read cmdline
  cmdline <- readxl::read_excel(path = path, col_names = FALSE, sheet = "cmdline")
  cmdline <- unlist(cmdline, use.names = FALSE)
  # For RP models, extract other stuff:
  if (family == "rp") {
    # Knots:
    knots <- readxl::read_excel(path = path, col_names = FALSE, sheet = "e(knots1)")
    knots <- unlist(knots, use.names = FALSE)
    knots <- stringr::str_split_1(string = knots, pattern = " ")
    knots <- as.numeric(knots)
    # Orthog:
    orthog <- readxl::read_excel(path = path, col_names = FALSE, sheet = "e(orthog1)")
    orthog <- unlist(orthog, use.names = FALSE)
    if (is.na(orthog)) orthog <- ""
    # RCS rmat:
    if (orthog != "") {
      rcsrmat <- readxl::read_excel(path = path, col_names = FALSE, sheet = "e(rcsrmat_1)")
      rcsrmat <- as.matrix(rcsrmat)
    }
  }
  # Return
  out <- list(eb = eb, eV = eV, cmd = cmd, cmdline = cmdline, family = family)
  if (cmd == "mestreg") {
    out$frm <- frm
  }
  if (family == "rp") {
    out$knots <- knots
    out$orthog <- orthog
    if (orthog != "") {
      out$rcsrmat <- rcsrmat
    }
  }
  return(out)
}
