library(haven)
library(tidyverse)
library(readxl)
library(ragg)
library(mvtnorm)
library(matrixStats)
library(furrr)

## Common number of "parametric bootstrap" replications
.B <- 200

## Load data
dt <- read_dta(file = "~/Documents/proj-bsw-1cvd/data/91-data3CIA-pp.dta") |>
  zap_formats() |>
  zap_label() |>
  zap_labels() |>
  mutate(`0b.mmrc0_4` = as.numeric(mmrc0_4 == 0)) |>
  mutate(`1.mmrc0_4` = as.numeric(mmrc0_4 == 1)) |>
  mutate(`2.mmrc0_4` = as.numeric(mmrc0_4 == 2)) |>
  mutate(`3.mmrc0_4` = as.numeric(mmrc0_4 == 3)) |>
  mutate(`4.mmrc0_4` = as.numeric(mmrc0_4 == 4))

## Load estimation results
estimation_results = read_e(
  path_eb = "~/Documents/proj-bsw-1cvd/data/91-eb.xlsx",
  path_eV = "~/Documents/proj-bsw-1cvd/data/91-eV.xlsx"
  )

## Times for predictions
.times <- seq(0, max(dt$follow_up_months), length.out = 100)

# Usage in our settings:
modm <- Stata.model.matrix(~ age + fev1pp + `0b.mmrc0_4` + `1.mmrc0_4` + `2.mmrc0_4` + `3.mmrc0_4` + `4.mmrc0_4`, data = dt, eb = estimation_results$eb)
