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
#' data("data3CIA", package = "stdmest")
"data3CIA"

#' @title Simulated Three-Level Data
#'
#' @description A simulated dataset with a three-levels hierarchical structure and
#'     a time to event outcome.
#'
#' @format A data frame the following columns:
#' * `hospital_id`, ID of the 3^rd^ level of hierarchy (highest);
#' * `provider_id`, ID of the 2^nd^ level of hierarchy;
#' * `patient_id`, ID of the 1^st^ level of hierarchy (lowest);
#' * `X1`, continuous covariate recorded at a patient-level;
#' * `X2`, a second continuous covariate recorded at a patient-level;
#' * `X3`, a binary covariate recorded at a patient-level;
#' * `t`, time to event data;
#' * `d`, event indicator variable, where 1 = event and 0 = censoring.
#'
#' @keywords datasets
#'
#' @examples
#' data("data3Lsim", package = "stdmest")
"data3Lsim"
