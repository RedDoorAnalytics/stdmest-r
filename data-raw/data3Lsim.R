###
### Simulated data for a three-levels model
###

## Packages
library(tidyverse)
library(glue)
library(matrixStats)
library(usethis)
library(haven)

## Seed
set.seed(945738)

## Setup
n <- 25000 # subjects
m <- 1000 # hospitals
s <- 10 # surgeons per hospital

## Base data
base <- tibble(id = seq(n))
base <- base |>
  mutate(
    X1 = rnorm(n = n),
    X2 = rnorm(n = n, sd = 5),
    X3 = rbinom(n = n, size = 1, prob = 0.4)
  )

## For multinomial logistic cluster assignment, see:
## --> https://library.virginia.edu/data/articles/simulating-multinomial-logistic-regression-data

## Scenarios based on Chen et al. (2020), published in SMMR
gammas <- rnorm(n = m, sd = sqrt(0.5))
phi1s <- rnorm(n = m, sd = sqrt(0.1))
phi2s <- rnorm(n = m, sd = sqrt(0.2))
phi3s <- rnorm(n = m, sd = sqrt(0.3))
lps <- map(.x = 2:m, .f = function(j) gammas[j] + phi1s[j] * base$X1 + phi2s[j] * base$X2 + phi3s[j] * base$X3)
lps <- do.call(cbind, lps)
den <- 1 + rowSums2(exp(lps))
p1 <- 1 / den
ps <- map(.x = 1:ncol(lps), .f = function(j) exp(lps[, j]) / den)
ps <- do.call(cbind, ps)
ps <- cbind(p1, ps)
base$hospital <- apply(ps, MARGIN = 1, function(x) sample(x = seq(ncol(ps)), size = 1, prob = x))

## Then do the same, but for surgeons
gammas <- rnorm(n = s, sd = sqrt(0.5))
phi1s <- rnorm(n = s, sd = sqrt(1.5))
phi2s <- rnorm(n = s, sd = sqrt(1.5))
phi3s <- rnorm(n = s, sd = sqrt(2))
lps <- map(.x = 2:s, .f = function(j) gammas[j] + phi1s[j] * base$X1 + phi2s[j] * base$X2 + phi3s[j] * base$X3)
lps <- do.call(cbind, lps)
den <- 1 + rowSums2(exp(lps))
p1 <- 1 / den
ps <- map(.x = 1:ncol(lps), .f = function(j) exp(lps[, j]) / den)
ps <- do.call(cbind, ps)
ps <- cbind(p1, ps)
base$surgeon <- apply(ps, MARGIN = 1, function(x) sample(x = seq(ncol(ps)), size = 1, prob = x))

## Random effects
bh <- distinct(base, hospital)
bh <- mutate(bh, bh = rnorm(n = nrow(bh), sd = sqrt(.5)))
bs <- distinct(base, surgeon)
bs <- mutate(bs, bs = rnorm(n = nrow(bs), sd = sqrt(.1)))

## Survival parameters
lambda <- 0.1
p <- 2

## Simulate survival outcomes
dt <- left_join(base, bs, by = "surgeon") |>
  left_join(bh, by = "hospital") |>
  mutate(u = runif(n = nrow(base))) |>
  mutate(nu = X1 * 0.1 + X2 * 0.2 + X3 * (-0.5) + bs + bh) |>
  mutate(true = (-log(u) / (lambda * exp(nu)))^(1 / p)) |>
  mutate(t = pmin(true, 10)) |>
  mutate(d = as.numeric(t < 10)) |>
  arrange(hospital, surgeon, id)

## Make surgeon ID that is unique overall (not only within cluster)
dt <- dt |>
  mutate(surgeon = as.numeric(factor(paste0(hospital, "_", surgeon))))

## Rename IDs
dt <- dt |>
  rename(hospital_id = hospital, provider_id = surgeon, patient_id = id)

## Export
data3Lsim <- dt |>
  select(hospital_id, provider_id, patient_id, X1, X2, X3, t, d) |>
  arrange(hospital_id, provider_id, patient_id)
use_data(data3Lsim, overwrite = TRUE)
write_dta(data = data3Lsim, path = "data-raw/data3Lsim.dta")
