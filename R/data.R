#' @title Individual Participant Data Meta-Analysis of COPD Patients
#'
#' @description A random subset of the 3CIA database of COPD patients.
#'
#' @format A data frame the following columns:
#' * `months`, length of follow-up in months;
#' * `status`, indicator of the status at the end of follow-up (where
#'   0: alive, 1: death);
#' * `age`, age in years;
#' * `fev1pp`, FEV1 measurement;
#' * `mmrc`, dyspnea score (mMRC);
#' * `cohort`, study indicator.
#'
#' @references
#' * J. B. Soriano, B. Lamprecht, A. S. Ramirez et al. Mortality prediction in
#'   chronic obstructive pulmonary disease comparing the GOLD 2007 and 2011 staging
#'   systems: A pooled analysis of individual patient data. The Lancet Respiratory
#'   Medicine, 3(6):443--450, 2015. <doi:10.1016/S2213-2600(15)00157-5>
#'
#' @keywords datasets
#'
#' @examples
#' data("data3CIA", package = "stmepred")
"data3CIA"
