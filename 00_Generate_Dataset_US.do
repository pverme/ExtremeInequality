global pc "/Users/lidia/Dropbox/Collaborazioni/Ceriani_Verme_BottomIncomes"
global data "${pc}/Analysis/Data"
global figures "${pc}/Paper/Figures"
global data_new "${pc}/Analysis/Data_new"
*===============================================================================
* Ceriani - Verme - Bottom Incomes 
*===============================================================================

use "${data}/CPS2016.dta", clear
rename h_year year
gen y = htotval/h_numper
decode gestfips, gen(state)
levelsof state, local(state)
levelsof year, local(year)



*-------------------------------------------------------------------------------
* Summary Stats
*-------------------------------------------------------------------------------

foreach y of local year{
matrix M_`y' = . , . , . , . , . , . , . , . , . , . , . , . 
		foreach s of local state{
			preserve
			keep if state=="`s'" & year==`y'
			count if y<0
			
			matrix Neg_`s' = r(N)
			drop if y<0
			
			sum y [aw=hsup_wgt], det
			matrix Max_`s' = (r(max)/r(mean))			
			matrix Min_`s' = (r(min)/r(mean))			
			matrix Obs_`s'  = r(N)
			matrix Pop_`s'  = r(sum_w)
			local  mean = r(mean)
			local  sd = r(sd)
			local  min = r(min)
			local  max = r(max)
			
			conindex y [fw=round(hsup_wgt)], truezero	
			
			matrix M_`y'_`s' = `y', Obs_`s', Pop_`s', Neg_`s',`mean', `sd', `min', `max', Min_`s' , Max_`s' ,  r(CI), r(CIse)
			
			matrix rownames M_`y'_`s'= `s'
			matrix M_`y' = M_`y' \ M_`y'_`s'
			restore
		}
}

matrix colnames M_2016 = year obs wgt neg mean sd  min max gamma lambda gini se_gini

putexcel set "${pc}/Table.xlsx", sheet("M_US_Stats", replace)  modify
putexcel A1 = matrix(M_2016), names



*-------------------------------------------------------------------------------
* Import Matrix with relevant information
*-------------------------------------------------------------------------------
import excel "${pc}/Table.xlsx", sheet("M_US_Stats") cellrange(A1:M53) firstrow  clear
drop if year==.
rename A state

* H1
gen beta_H1 = ((1-gamma)-gini*(gamma+1))/(gini*gamma)
gen gini_H1 = (1-gamma)/(1+gamma)

* H2
gen tau_H2 = ((lambda-1)-gini*(lambda+1))/(gini*lambda)
gen gini_H2 = (lambda-1)/(lambda+1)

* H3
gen tau_ex = 0.01
gen A = gini+tau_ex*(lambda-1)
gen B = (1+tau_ex)*(1+tau_ex*lambda)
gen beta_H3 = (A-gini*B)/(B-A)
gen tau_H3  = ((lambda - 1) - gini*(lambda+1))/gini*lambda
gen gini_H3 = (lambda-1)/(lambda+1)

forvalues i=1(1)100{
gen beta_ex_`i' = `i'/100

* H4
gen gini_H4_`i' = (tau_ex*(lambda-1)+beta_ex_`i'*(1-gamma))/((1+beta_ex_`i'+tau_ex)*(1+beta_ex_`i'*gamma+tau_ex*lambda)-1)
gen gini_BT_`i' = (gini+tau_ex*(lambda-1)+beta_ex_`i'*(1-gamma))/((1+beta_ex_`i'+tau_ex)*(1+beta_ex_`i'*gamma+tau_ex*lambda))
}


forvalues i=1(1)100{
assert gini_BT_`i'>gini
}

forvalues j=1(1)20{
gen tau_ex_`j' = `j'/100
gen gini_T_`j' = (gini+tau_ex_`j'*(lambda-1))/((1+tau_ex_`j')*(1+tau_ex_`j'*lambda))
	forvalues i=1(1)20{
	gen gini_BT_`i'_`j' = (gini+tau_ex_`j'*(lambda-1)+beta_ex_`i'*(1-gamma))/((1+beta_ex_`i'+tau_ex_`j')*(1+beta_ex_`i'*gamma+tau_ex_`j'*lambda))	
}
}

forvalues j=1(1)20{
	forvalues i=1(1)20{
		assert gini_BT_`i'_`j'< gini_T_`j'
		}
	}	

	order state year gini mean gamma lambda *_H1 *_H2 A B *_H3 *_H4_* gini_BT* beta_*

export excel using "${pc}/Table.xlsx", sheet("Hypotheses", replace) firstrow(var) 

*-------------------------------------------------------------------------------
* Generate new distribution adding to the bottom and to the top
*-------------------------------------------------------------------------------

levelsof state, local(state)
	foreach s of local state{
		preserve
		keep if state=="`s'"
		drop if y<0
		

		gen y_B = y
		gen y_T = y

		gen wgt_B = hsup_wgt
		gen wgt_T = hsup_wgt

		count 
		local obs = r(N) + 1
			
		sum y [aw=hsup_wgt], det

		local min = r(min)
		local max = r(max)
		local pop = r(sum_w)
					
		local b  = 0.03*`pop'			
		local t  = 0.01*`pop'

		set obs `obs'	
	
		* add to the bottom			
		replace y_B = `min' if y == .
		replace wgt_B = `b' if wgt_B==.
		
		* add to the top
		replace y_T = `max' if y == .
		replace wgt_T = `t' if wgt_T==.
		
		* add state name and year
		replace state = "`s'" if state==""
		replace year = 2016 if year==.
	
*-------------------------------------------------------------------------------
* Generate new distribution adding to the bottom and to the top			*-------------------------------------------------------------------------------
		gen y_BT = y_B
		gen wgt_BT = wgt_B

		count 
		local obs = r(N) + 1
					
		* add to the bottom
		set obs  `obs'
					
		replace y_BT = `max' if y_BT == .
		replace wgt_BT = `t' if wgt_BT==.
		
		* add state name
		replace state = "`s'" if state==""
		replace year = 2016 if year==.
		tempfile data_`s'
		save `data_`s''
restore
}

*-------------------------------------------------------------------------------
* Append together datasets
*-------------------------------------------------------------------------------*
levelsof state, local(state)
use `data_AK', clear
gen temp = 1		
	foreach s of local state{
	append using `data_`s''	
}
drop if temp==1
drop temp

save "${data_new}/Dataset.dta", replace
