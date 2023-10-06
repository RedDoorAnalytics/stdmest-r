library(haven)
library(tidyverse)
library(readxl)
library(ragg)
library(mvtnorm)
library(matrixStats)
library(furrr)
devtools::load_all()

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
estimation_results <- read_e(
  path_eb = "~/Documents/proj-bsw-1cvd/data/91-eb.xlsx",
  path_eV = "~/Documents/proj-bsw-1cvd/data/91-eV.xlsx"
)

## Times for predictions
.times <- seq(0, max(dt$follow_up_months), length.out = 200)

## Usage in our settings:
modm <- Stata.model.matrix(
  fixed = ~ age + fev1pp + `0b.mmrc0_4` + `1.mmrc0_4` + `2.mmrc0_4` + `3.mmrc0_4` + `4.mmrc0_4`,
  random = ~ b - 1,
  data = dt,
  eb = estimation_results$eb
)

best <- filter(dt, b == min(b)) |>
  distinct(b, bse)

a <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = best$b, bse = best$bse, distribution = "weibull", contrast = TRUE, conf.int = TRUE)
a |>
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = S_conf.low, ymax = S_conf.high, fill = "S"), alpha = 0.1) +
  geom_ribbon(aes(ymin = Sref_conf.low, ymax = Sref_conf.high, fill = "Sref"), alpha = 0.1) +
  geom_line(aes(y = S, color = "S")) +
  geom_line(aes(y = Sref, color = "Sref"))
a |>
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = Sdiff_conf.low, ymax = Sdiff_conf.high), alpha = 0.1) +
  geom_line(aes(y = Sdiff))
