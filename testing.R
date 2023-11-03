library(haven)
library(tidyverse)
library(readxl)
library(ragg)
library(mvtnorm)
library(matrixStats)
library(furrr)
library(fastGHQuad)
library(ggtext)
devtools::load_all()
options(scipen = 100)

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
dt_to_use <- subset(dt, dt$months > 0)
modm <- Stata.model.matrix(
  fixed = ~ age + fev1pp + `0b.mmrc` + `1.mmrc` + `2.mmrc` + `3.mmrc` + `4.mmrc`,
  random = ~ b - 1,
  data = dt_to_use,
  eb = estimation_results$eb
)

# Comparison with Stata
set.seed(20231018)
Smin <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = -1.006262, bse = .2222539, distribution = "weibull", conf.int = TRUE, B = 2000)
Smin2 <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed[dt_to_use$cohort == 18, ], Sigma = estimation_results$eV, b = -1.006262, bse = .2222539, distribution = "weibull", conf.int = TRUE, B = 2000)
data.frame(tt = .times, Smin = Smin$S, Smin_lower = Smin$S_conf.low, Smin_upper = Smin$S_conf.high, Smin2 = Smin2$S, Smin2_lower = Smin2$S_conf.low, Smin2_upper = Smin2$S_conf.high)
#    tt      Smin Smin_lower Smin_upper     Smin2 Smin2_lower Smin2_upper
# 1   0 1.0000000  1.0000000  1.0000000 1.0000000   1.0000000   1.0000000
# 2  50 0.9319127  0.8929531  0.9571315 0.9245237   0.8830277   0.9533136
# 3 100 0.8373119  0.7583575  0.8940225 0.8198216   0.7327097   0.8853843
# 4 150 0.7428443  0.6347484  0.8260286 0.7159733   0.5979398   0.8129497
# 5 200 0.6561102  0.5315389  0.7571369 0.6216718   0.4896364   0.7419416

###
### Drawing multivariate normal samples using SVD
###

n.values <- 1e6
eb <- estimation_results$eb
eV <- estimation_results$eV
std_norm <- matrix(data = rnorm(ncol(eV) * n.values), nrow = n.values, ncol = ncol(eV))
svdd <- svd(eV)
rand_data_svd <- t(svdd$u %*% (diag(sqrt(svdd$d))) %*% t(std_norm))
for (i in 1:ncol(eV)) {
  rand_data_svd[, i] <- rand_data_svd[, i] + eb[i]
}
cov(rand_data_svd) - eV

###
### Three-levels example
###

## Load data
dt <- read_dta(file = "data-raw/data3Lsim-pp.dta") |>
  zap_formats() |>
  zap_label() |>
  zap_labels() |>
  mutate(`0b.X3` = as.numeric(X3 == 0)) |>
  mutate(`1.X3` = as.numeric(X3 == 1)) |>
  filter(`_st` == 1)

## Load estimation results
estimation_results <- read_e(
  path_eb = "data-raw/data3Lsim-eb.xlsx",
  path_eV = "data-raw/data3Lsim-eV.xlsx"
)

## Times for predictions
.times <- seq(0, max(dt$t), length.out = 5)

## Usage in our settings:
modm <- Stata.model.matrix(
  fixed = ~ X1 + X2 + `0b.X3` + `1.X3`,
  random = ~ b_hospital + b_provider - 1,
  data = dt,
  eb = estimation_results$eb
)

set.seed(349856)
reffs <- distinct(dt, hospital_id, b_hospital, bse_hospital, provider_id, b_provider, bse_provider) |>
  filter(b_hospital %in% sample(unique(dt$b_hospital), 2)) |>
  arrange(b_hospital, b_provider) |>
  group_by(hospital_id) |>
  slice_sample(n = 1) |>
  ungroup()
reffs
# # A tibble: 2 × 6
#   hospital_id b_hospital bse_hospital provider_id b_provider bse_provider
#         <dbl>      <dbl>        <dbl>       <dbl>      <dbl>        <dbl>
# 1          88      0.878        0.311        4420      0.122        0.570
# 2         397      2.65         0.407        3078     -0.298        0.522
a1 <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = c(0.878, 0.122), bse = c(0.311, 0.570), bref = c(0, 0), brefse = c(0, 0), distribution = "weibull", conf.int = TRUE, B = 2000, cimethod = "normal")
a2 <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = c(2.65, -0.298), bse = c(0.407, 0.522), bref = c(0, 0), brefse = c(0, 0), distribution = "weibull", conf.int = TRUE, B = 2000, cimethod = "normal")
data.frame(tt = .times, Sa1 = a1$S, Sa1_lower = a1$S_conf.low, Sa1_upper = a1$S_conf.high, Sa2 = a2$S, Sa2_lower = a2$S_conf.low, Sa2_upper = a2$S_conf.high)
#     tt         Sa1   Sa1_lower  Sa1_upper           Sa2    Sa2_lower   Sa2_upper
# 1  0.0 1.000000000  1.00000000 1.00000000 1.00000000000  1.000000000 1.000000000
# 2  2.5 0.206141020 -0.05329716 0.46557920 0.02954803823 -0.079408427 0.138504503
# 3  5.0 0.027802364 -0.06700799 0.12261272 0.00114208034 -0.015675085 0.017959246
# 4  7.5 0.004786736 -0.02944316 0.03901663 0.00009295082 -0.003372755 0.003558657
# 5 10.0 0.001041533 -0.01251203 0.01459510 0.00001306697 -0.000879363 0.000905497

###
### Three-levels example, partially marginal
###

# First, calculate some "fully conditional" predictions
set.seed(934867)
.times <- seq(0, max(dt$t), length.out = 100)
reffs <- distinct(dt, hospital_id, b_hospital, bse_hospital, provider_id, b_provider, bse_provider) |>
  filter(b_hospital %in% sample(unique(dt$b_hospital), 9))
plan(multisession)
a <- future_map(.x = seq(nrow(reffs)), .f = function(i) {
  preds <- stdmest(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = c(reffs$b_hospital[i], reffs$b_provider[i]), bse = c(reffs$bse_hospital[i], reffs$bse_provider[i]), bref = c(0, 0), brefse = c(0, 0), distribution = "weibull", contrast = TRUE, conf.int = TRUE)
  preds$hospital_id <- reffs$hospital_id[i]
  preds$provider_id <- reffs$provider_id[i]
  return(preds)
}, .progress = TRUE, .options = furrr_options(seed = 32475))
plan(sequential)
a <- bind_rows(a)

#
reffsm <- distinct(reffs, hospital_id, b_hospital, bse_hospital)

plan(multisession)
am <- future_map(.x = seq(nrow(reffsm)), .f = function(i) {
  out <- stdmestm(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = reffsm$b_hospital[i], bse = reffsm$bse_hospital[i], bref = 0, brefse = 0, varmargname = "var(_cons[hospital_id>provider_id])", distribution = "weibull", contrast = TRUE, conf.int = TRUE)
  out$hospital_id <- reffsm$hospital_id[i]
  out$b_hospital <- reffsm$b_hospital[i]
  return(out)
}, .progress = TRUE, .options = furrr_options(seed = 45986))
plan(sequential)
am <- bind_rows(am)

a |>
  ggplot(aes(x = t, group = provider_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = Sdiff_conf.low, ymax = Sdiff_conf.high), alpha = 0.05) +
  geom_line(aes(y = Sdiff)) +
  geom_ribbon(data = am, aes(ymin = Sdiff_conf.low, ymax = Sdiff_conf.high, group = hospital_id), alpha = 0.1, fill = "red") +
  geom_line(data = am, aes(y = Sdiff, group = hospital_id), color = "red") +
  facet_wrap(~hospital_id, labeller = label_both) +
  theme_bw(base_size = 12, base_family = "Atkinson Hyperlegible") +
  theme(plot.title = element_textbox_simple()) +
  labs(x = "Time", y = "Standardised Survival Difference", title = "Black lines denote provider-specific standardised survival differences (with the theoretical overall average as the reference). Red lines denote the hospital-specific standardised survival differences, marginally over providers.")
ggsave(filename = "testing.png", device = ragg::agg_png, width = 8, height = 7, dpi = 300)

###
### Code for comparing with Stata's implementation
###

.times <- seq(0, 10, length.out = 5)
out <- stdmestm(t = matrix(.times, ncol = 1), beta = estimation_results$eb, X = modm$fixed, Sigma = estimation_results$eV, b = 0.878, bse = 0.311, bref = 0, brefse = 0, varmargname = "var(_cons[hospital_id>provider_id])", distribution = "weibull", contrast = TRUE, conf.int = TRUE, B = 2000, cimethod = "normal")
out
#      t           S    S_conf.low S_conf.high       Sref Sref_conf.low Sref_conf.high       Sdiff Sdiff_conf.low Sdiff_conf.high
# 1  0.0 1.000000000  1.0000000000  1.00000000 1.00000000    1.00000000     1.00000000  0.00000000     0.00000000     0.000000000
# 2  2.5 0.252036599  0.1164155373  0.38765766 0.45942711    0.40991674     0.50893748 -0.20739051    -0.33674710    -0.078033924
# 3  5.0 0.052759040  0.0002934933  0.10522459 0.15662493    0.12366405     0.18958581 -0.10386589    -0.15691719    -0.050814596
# 4  7.5 0.013903120 -0.0054939075  0.03330015 0.05815525    0.04131102     0.07499948 -0.04425213    -0.06597111    -0.022533158
# 5 10.0 0.004391529 -0.0033833830  0.01216644 0.02382806    0.01537816     0.03227796 -0.01943654    -0.02917959    -0.009693476
