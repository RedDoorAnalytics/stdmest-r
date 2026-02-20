// Require Stata's stdmest
// net install stdmest, from("https://raw.githubusercontent.com/RedDoorAnalytics/stdmest/main/")

// Three levels: simulated data
use "data-raw/data3Lsim", clear

// stset
stset t, failure(d)

// mestreg Weibull model
mestreg c.X1 c.X2 i.X3 || hospital_id: || provider_id:, dist(weibull)

// export e(b) and e(V)
modexpt, filename("data-raw/data3Lsim-ebV.xlsx") replace

// predict random effect, with their SEs, and export that
predict b_hospital b_provider, reffects reses(bse_hospital bse_provider)
compress
save "data-raw/data3Lsim-pp.dta", replace
