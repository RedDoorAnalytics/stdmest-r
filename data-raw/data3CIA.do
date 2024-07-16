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

// RP(3), from -stmixed-, with orthogonalised splines
* Adjusting for age, FEV1, and dyspnea score
tabulate mmrc, generate(immrc)
stmixed age fev1pp immrc2 immrc3 immrc4 immrc5 || cohort: , dist(rp) df(3)
* export e(b) and e(V)
modexpt, filename("data-raw/data3CIA-ebV-rp3-orthog.xlsx") replace
* predict random effect, with their SEs, and export that
predict b_rp3orthog, reffects
* SEs not currently supported by -stmixed-: will give SEs of zero (for now)
generate bse_rp3orthog = 0.0
drop _rcs*

// RP(3), from -stmixed-, with non-orthogonalised splines
* Adjusting for age, FEV1, and dyspnea score
stmixed age fev1pp immrc2 immrc3 immrc4 immrc5 || cohort: , dist(rp) df(3) noorthog
* export e(b) and e(V)
modexpt, filename("data-raw/data3CIA-ebV-rp3-noorthog.xlsx") replace
* predict random effect, with their SEs, and export that
predict b_rp3noorthog, reffects
* SEs not currently supported by -stmixed-: will give SEs of zero (for now)
generate bse_rp3noorthog = 0.0
drop _rcs*

// Export
compress
save "data-raw/data3CIA-pp.dta", replace
