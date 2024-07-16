// Working directory
cd "~/R-dev/stdmest"
adopath ++ "~/Stata-dev/stdmest/modexpt"

// Two levels: IPD meta-analysis
use "data-raw/data3CIA", clear

// stset
stset months, failure(status == 1)

// Weibull PH, from -mestreg-
* Adjusting for age, FEV1, and dyspnea score
mestreg c.age c.fev1pp ib0.mmrc || cohort: , dist(weibull) nohr
* export e(b) and e(V)
modexpt, filename("data-raw/data3CIA-ebV-weibull.xlsx") replace
* predict random effect, with their SEs, and export that
predict b_wei, reffects reses(bse_wei)

// Log-logistic AFT, from -mestreg-
* Adjusting for age, FEV1, and dyspnea score
mestreg c.age c.fev1pp ib0.mmrc || cohort: , dist(loglogistic) time
* export e(b) and e(V)
modexpt, filename("data-raw/data3CIA-ebV-logl.xlsx") replace
* predict random effect, with their SEs, and export that
predict b_logl, reffects reses(bse_logl)

// Export
compress
save "data-raw/data3CIA-pp.dta", replace
