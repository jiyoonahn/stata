capture log close _all
log using ipw_bi.log, replace

clear

use cohort
// load Stata dataset
// if you want to load a SAS file, would use--import sas filename.sas7bdat

global covar "rec_age65 rec_female ib1.rec_bmic"
// Covariates that you would want to adjust for is assigned to 'covar'.
// ib2.rec_eth: rec_eth=2 would be reference


/*-----------------------------
Estimate stabilized inverse probability weights; binary treatments
- treatment(exposure): atg (1: anti-thymocyte globulin, 0: basiliximab)
- covariates (L)
  rec_age65 (1: 65+, 0: 18-64)
  rec_female (1: female, 0: male)
  ib2.rec_eth (1: white, 2: black [reference], 3:other, 4:hispanic)
  rec_bmic: categorized BMI (1: <20, 2: 20-29.9, 3: 30+)
  rec_bmic_miss: missing indicator for rec_bmic
-----------------------------*/
***** Estimate probability of getting treated with atg conditional on covariates
***** P[atg=1 | L]
logit atg $covar
predict p_atg, pr

***** Estimate probability of getting treated with atg
***** P[atg=1]
logit atg
predict _p, pr
// if interested in EMM by age, P[atg=1 | rec_age65] need to be used as below
*logit atg rec_age65
*predict _p, pr

***** Calculate stabilized IPW
gen sipw = _p / p_atg if atg == 1
replace sipw = (1 - _p) / (1 - p_atg) if atg == 0
summarize sipw
// check if the mean of SIPW is close to 1!

***** Check positivity assumption: propensity score distribution
kdensity p_atg if atg == 1, lc("gs4") lpattern(solid) lwidth(medthick) plot(kdensity p_atg if atg == 0, lc("gs9") lpattern(dash) lwidth(medthick)) legend(order(1 "rATG" 2 "Basiliximab") col(1) pos(10) ring(0)) graphregion(color(white)) title("") xtitle("Propensity score") xscale(r(0(0.2)1)) xlabel(.0 0.2 0.4 0.6 0.8 1.0)
graph export psdist_temp.pdf, replace

***** Check balance after weighted
bys atg: summarize $covar [iweight = sipw]



/*-----------------------------
Estimate hazard ratio using Cox model
-----------------------------*/
***** Setting time-to-event data
capture drop event
capture drop exit
gen exit = min(t1, clm_from, tdate + 365.25 * 3)
replace exit = exit + 0.1 if exit == tdate
gen event = inrange(clm_from, tdate, exit)

***** Unadjusted
stset exit, origin(tdate) f(event) scale(365.25)
stcox atg
// Overall
stcox atg##rec_age65
// EMM by age

***** Weighted
stset exit [pweight = sipw], origin(tdate) f(event) scale(365.25)
stcox atg
stcox atg##rec_age65



log close _all
