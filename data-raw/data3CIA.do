// Working directory
cd "~/R-dev/stdmest"

// Two levels: IPD meta-analysis
use "data-raw/data3CIA", clear

// stset
stset months, failure(status == 1)

// mestreg Weibull model adjusting for age, FEV1, and dyspnea score
// we loop over PH distributions to get 'best fitting' model
mestreg c.age c.fev1pp ib0.mmrc || cohort: , dist(weibull) nohr

// export e(b) and e(V)
*
matrix eb = e(b)
putexcel set "data-raw/data3CIA-eb.xlsx", replace
putexcel A1 = matrix(eb), names
*
matrix eV = e(V)
putexcel set "data-raw/data3CIA-eV.xlsx", replace
putexcel A1 = matrix(eV)

// predict random effect, with their SEs, and export that
predict b, reffects reses(bse)
compress
save "data-raw/data3CIA-pp.dta", replace
