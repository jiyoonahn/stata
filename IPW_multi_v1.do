capture log close _all
log using ipw_multi.log, replace

clear

use cohort
// load Stata dataset
// if you want to load a SAS file, would use--import sas filename.sas7bdat

global covar "rec_age65 rec_female ib1.rec_bmic"
// Covariates that you would want to adjust for is assigned to 'covar'.
// ib2.rec_eth: rec_eth=2 would be reference


/*-----------------------------
Estimate stabilized inverse probability weights; multiple treatments
- treatment(exposure): atgsw 
					   (1: basiliximab; continued steroid maintenance (CSM)
					    2: atg; CSM
						3: basiliximab; early steroid withdrawal (ESW)
						4: atg; ESW) 
- covariates (L)
  rec_age65 (1: 65+, 0: 18-64)
  rec_female (1: female, 0: male)
  ib2.rec_eth (1: white, 2: black [reference], 3:other, 4:hispanic)
  rec_bmic: categorized BMI (1: <20, 2: 20-29.9, 3: 30+)
  rec_bmic_miss: missing indicator for rec_bmic
-----------------------------*/
gen atgsw = .
replace atgsw = 1 if atg == 0 & sw == 0
replace atgsw = 2 if atg == 1 & sw == 0
replace atgsw = 3 if atg == 0 & sw == 1
replace atgsw = 4 if atg == 1 & sw == 1

***** Estimate probability of receiving each treatment conditional on covariates
***** P[atgsw=1 | L] P[atgsw=2 | L] P[atgsw=3 | L] P[atgsw=4 | L]
mlogit atgsw $covar, nolog
predict p1 p2 p3 p4, pr
list atgsw p1-p4 in 1/10

***** Estimate probability of receiving each treatment
***** P[atgsw=1] P[atgsw=2] P[atgsw=3] P[atgsw=4]
mlogit atgsw, nolog
predict _p1 _p2 _p3 _p4, pr
// if interested in EMM by age, P[atgsw | rec_age65] need to be used as below
*mlogit atgsw rec_age65, nolog
*predict _p1 _p2 _p3 _p4, pr

gen sipw = .
forvalues y = 1 / 4 {
	replace sipw = _p`y' / p`y' if atgsw == `y'
}
summarize sipw



/*-------------------------------
Check positivity assumption
: Distribution of propensity score of each treatment 
-------------------------------*/
graph box p1, over(atgsw)
graph export psdist_atgsw_1.pdf, replace
graph box p2, over(atgsw)
graph export psdist_atgsw_2.pdf, replace
graph box p3, over(atgsw)
graph export psdist_atgsw_3.pdf, replace
graph box p4, over(atgsw)
graph export psdist_atgsw_4.pdf, replace

***** Check balance after weighted
bys atgsw: summarize $covar [iweight = sipw]




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
stcox ib1.atgsw
// Overall
stcox ib1.atg##rec_age65
// EMM by age

***** Weighted
stset exit [pweight = sipw], origin(tdate) f(event) scale(365.25)
stcox ib1.atgsw
stcox ib1.atgsw##rec_age65



log close _all
