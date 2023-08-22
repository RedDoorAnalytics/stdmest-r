## The dataset is a !!!_random subset_!!! of the 3CIA database
## of COPD patients
##
## Variable definitions
##  follow_up_months: length of follow-up in months
##  alive1_death2: indicator of the status at the end of follow-up
##                 (1: alive, 2: death)
##  age: age in years
##  fev1pp: FEV1 measurement
##  mmrc: dyspnea score (mMRC)
##  ages, fev1pps, mmrcs: rescaled versions of the three covariates
##                        (using their ranges)
##  cohort: study indicator
##  su: Surv object treating survival times as right-censored
##  sui: Surv object taking interval censoring into account

## Reference for the data:
## J. B. Soriano, B. Lamprecht, A. S. Ramirez, P. Martinez-Camblor, B. Kaiser,
## I. Alfageme, P. Almagro, C. Casanova, C. Esteban, J. J. Soler-Cataluna,
## J. P. de Torres, M. Miravitlles, B. R. Celli, J. M. Marin, M. A. Puhan,
## P. Sobradillo, P. Lange, A. L. Sternberg, J. Garcia-Aymerich,
## A. M. Turner, M. K. Han, A. Langhammer, L. Leivseth, P. Bakke,
## A. Johannessen, N. Roche, and D. D. Sin. -- Mortality prediction in chronic
## obstructive pulmonary disease comparing the GOLD 2007 and 2011 staging
## systems: A pooled analysis of individual patient data.
## The Lancet Respiratory Medicine, 3(6):443--450, 2015.
## <doi:10.1016/S2213-2600(15)00157-5>

## Packages
library(tidyverse)
library(tramME)

## Load data
load(system.file(file.path("demo-data", "ipd.rda"), package = "tramME"))

## Remove extra columns
data3CIA <- data3CIA |>
  select(-sui, -su, -mmrcs, -fev1pps, -ages, -mmrc0_4)

## Rename
data3CIA <- data3CIA |>
  rename(months = follow_up_months, status = alive1_death2)

## Recode status
data3CIA <- data3CIA |>
  mutate(status = status - 1)

## Reformat cohort indicator
data3CIA <- data3CIA |>
  mutate(cohort = as.numeric(cohort))

## Export in Stata format
usethis::use_data(data3CIA, overwrite = TRUE)
