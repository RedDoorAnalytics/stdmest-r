library(haven)
library(tidyverse)
library(readxl)
library(ragg)
library(mvtnorm)
library(matrixStats)
library(furrr)
devtools::load_all()

###
### Two-levels example
###

## Load data
dt <- read_dta(file = "data-raw/data3CIA-pp.dta") |>
  zap_formats() |>
  zap_label() |>
  zap_labels() |>
  mutate(`0b.mmrc` = as.numeric(mmrc == 0)) |>
  mutate(`1.mmrc` = as.numeric(mmrc == 1)) |>
  mutate(`2.mmrc` = as.numeric(mmrc == 2)) |>
  mutate(`3.mmrc` = as.numeric(mmrc == 3)) |>
  mutate(`4.mmrc` = as.numeric(mmrc == 4))

## Load estimation results
estimation_results <- read_e(
  path_eb = "data-raw/data3CIA-eb.xlsx",
  path_eV = "data-raw/data3CIA-eV.xlsx"
)

## Times for predictions
.times <- seq(0, max(dt$months), length.out = 100)

## Usage in our settings:
modm <- Stata.model.matrix(
  fixed = ~ age + fev1pp + `0b.mmrc` + `1.mmrc` + `2.mmrc` + `3.mmrc` + `4.mmrc`,
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

## Comparison with Stata

## Times for predictions
.times <- seq(0, 200, length.out = 5)

## Usage in our settings:
modm <- Stata.model.matrix(
  fixed = ~ age + fev1pp + `0b.mmrc` + `1.mmrc` + `2.mmrc` + `3.mmrc` + `4.mmrc`,
  random = ~ b - 1,
  data = dt,
  eb = estimation_results$eb
)

# Comparison with Stata
Smin <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = -1.006262, bse = .2222539, distribution = "weibull", contrast = TRUE, conf.int = TRUE, B = 2000, cimethod = "normal")
Smin2 <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed[dt$cohort == 18, ], Sigma = estimation_results$eV, b = -1.006262, bse = .2222539, distribution = "weibull", contrast = TRUE, conf.int = TRUE, B = 2000, cimethod = "normal")
data.frame(tt = .times, Smin = Smin$S, Smin_lower = Smin$S_conf.low, Smin_upper = Smin$S_conf.high, Smin2 = Smin2$S, Smin2_lower = Smin2$S_conf.low, Smin2_upper = Smin2$S_conf.high)
#    tt      Smin Smin_lower Smin_upper     Smin2 Smin2_lower Smin2_upper
# 1   0 1.0000000  1.0000000  1.0000000 1.0000000   1.0000000   1.0000000
# 2  50 0.9318754  0.8997094  0.9640413 0.9245237   0.8874684   0.9615790
# 3 100 0.8372385  0.7688946  0.9055824 0.8198216   0.7416359   0.8980073
# 4 150 0.7427501  0.6470008  0.8384994 0.7159733   0.6077461   0.8242006
# 5 200 0.6560073  0.5421759  0.7698387 0.6216718   0.4948546   0.7484891

###
### Three-levels example
###

## Load data
dt <- read_dta(file = "data-raw/data3Lsim-pp.dta") |>
  zap_formats() |>
  zap_label() |>
  zap_labels() |>
  mutate(`0b.X3` = as.numeric(X3 == 0)) |>
  mutate(`1.X3` = as.numeric(X3 == 1))

## Load estimation results
estimation_results <- read_e(
  path_eb = "data-raw/data3Lsim-eb.xlsx",
  path_eV = "data-raw/data3Lsim-eV.xlsx"
)

## Times for predictions
.times <- seq(0, max(dt$t), length.out = 200)

## Usage in our settings:
modm <- Stata.model.matrix(
  fixed = ~ X1 + X2 + `0b.X3` + `1.X3`,
  random = ~ b_hospital + b_provider - 1,
  data = dt,
  eb = estimation_results$eb
)

reffs <- distinct(dt, hospital_id, b_hospital, bse_hospital, provider_id, b_provider, bse_provider) |>
  left_join(distinct(dt, hospital_id, b_hospital) |> mutate(rank_hospital = rank(b_hospital))) |>
  filter(rank_hospital <= 5 | rank_hospital >= max(rank_hospital) - 5)

plan(multisession)
a <- future_map(.x = seq(nrow(reffs)), .f = function(i) {
  out <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = c(reffs$b_hospital[i], reffs$b_provider[i]), bse = c(reffs$bse_hospital[i], reffs$bse_provider[i]), bref = c(0, 0), brefse = c(0, 0), distribution = "weibull", contrast = TRUE, conf.int = TRUE)
  out$hospital_id <- reffs$hospital_id[i]
  out$b_hospital <- reffs$b_hospital[i]
  out$provider_id <- reffs$provider_id[i]
  out$b_provider <- reffs$b_provider[i]
  return(out)
}, .progress = TRUE, .options = furrr_options(seed = 9943756))
plan(sequential)
a <- bind_rows(a)
a |>
  ggplot(aes(x = t, group = provider_id)) +
  geom_ribbon(aes(ymin = S_conf.low, ymax = S_conf.high, fill = "S"), alpha = 0.05) +
  geom_ribbon(aes(ymin = Sref_conf.low, ymax = Sref_conf.high, fill = "Sref"), alpha = 0.05) +
  geom_line(aes(y = S, color = "S")) +
  geom_line(aes(y = Sref, color = "Sref")) +
  facet_wrap(~hospital_id, labeller = label_both)
a |>
  ggplot(aes(x = t, group = provider_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_ribbon(aes(ymin = Sdiff_conf.low, ymax = Sdiff_conf.high), alpha = 0.1) +
  geom_line(aes(y = Sdiff)) +
  facet_wrap(~hospital_id, labeller = label_both)
# animated version
library(gganimate)
a |>
  ggplot(aes(x = t, group = provider_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_ribbon(aes(ymin = Sdiff_conf.low, ymax = Sdiff_conf.high), alpha = 0.1) +
  geom_line(aes(y = Sdiff)) +
  transition_states(hospital_id) +
  enter_fade() +
  exit_fade() +
  ease_aes("linear") +
  labs(title = "Hospital: {closest_state}")
