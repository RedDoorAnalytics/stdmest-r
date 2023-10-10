// Working directory
cd "~/R-dev/stdmest"

// Three levels: simulated data
use "data-raw/data3Lsim", clear

// stset
stset t, failure(d)

// mestreg Weibull model
mestreg c.X1 c.X2 i.X3 || hospital_id: || provider_id:, dist(weibull)

// export e(b) and e(V)
*
matrix eb = e(b)
putexcel set "data-raw/data3Lsim-eb.xlsx", replace
putexcel A1 = matrix(eb), names
*
matrix eV = e(V)
putexcel set "data-raw/data3Lsim-eV.xlsx", replace
putexcel A1 = matrix(eV)

// predict random effect, with their SEs, and export that
predict b_h b_h_s, reffects reses(bse_h bse_h_s)
compress
save "data-raw/data3Lsim-pp.dta", replace
