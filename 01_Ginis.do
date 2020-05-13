global pc "/Users/lidia/Dropbox/Collaborazioni/Ceriani_Verme_BottomIncomes"
global data "${pc}/Analysis/Data"
global figures "${pc}/Paper/Figures"
global data_new "${pc}/Analysis/Data_new"
*===============================================================================
* Ceriani - Verme - Bottom Incomes 
*===============================================================================

use "${data_new}/Dataset.dta", clear

* condindex (jacknife?)
* ineqdeco con bootstrap

rename hsup_wgt wgt_I
label var wgt_I "Initial distribution weights - hsup_wgt"
rename y y_I
label var y_I "Initial distribution incomes - y"

*-------------------------------------------------------------------------------
* Gini Stats
*-------------------------------------------------------------------------------

* Conindex
*-------------------------------------------------------------------------------

levelsof state, local(state)

mat drop _all

matrix M = . , . , . , . , . , . , . , .
	foreach s of local state{
		preserve
		keep if state=="`s'"
		drop if y_I<0
			
		local dist "I B T BT"
			foreach d of local dist{
				conindex y_`d' [fw=round(wgt_`d')], truezero	
				local gini_conindex_`d' = r(CI)
				local se_conindex_`d'   = r(CIse)
				matrix M_`s'_`d' = `gini_conindex_`d'' , `se_conindex_`d''
				matrix rownames M_`s'_`d'= `s'
			}			
			matrix M = M \ M_`s'_I ,  M_`s'_B ,  M_`s'_T ,  M_`s'_BT
			restore
		}

matrix colnames M = gini_I_con  se_I_con  gini_B_con  se_B_con  gini_T_con  se_T_con  gini_BT_con  se_BT_con
putexcel set "${pc}/Table.xlsx", sheet("M_US_conindex", replace)  modify
putexcel A1 = matrix(M), names


/* FastGini(1)
*-------------------------------------------------------------------------------

levelsof state, local(state)

mat drop _all

matrix M = . , . , . , . , . , . , . , .
	foreach s of local state{
		preserve
		keep if state=="`s'"
		drop if y_I<0
			
		local dist "I B T BT"
			foreach d of local dist{
				fastgini1 y_`d' [pw=wgt_`d'], jk
				local gini_fast_`d' = r(gini)
				local se_fast_`d'   = r(se)
				matrix M_`s'_`d' = `gini_fast_`d'', `se_fast_`d''
				matrix names M_`s'_`d'= `s'
			}			
			matrix M = M \ M_`s'_I ,  M_`s'_B ,  M_`s'_T ,  M_`s'_BT
			restore
		}

matrix colnames M = gini_I_fast se_I_fast  gini_B_fast  se_B_fast gini_T_fast  se_T_fast  gini_BT_fast  se_BT_fast
putexcel set "${pc}/Table.xlsx", sheet("M_US_fast", replace)  modify
putexcel A1 = matrix(M), names
*/


/* Ineqdeco with bootstrap errors
*-------------------------------------------------------------------------------

* Initial
cap program drop ineqdec0_se
	program ineqdec0_se, rclass
				qui ineqdec0 y_I [fw=round(wgt_I)]
				return scalar gini = r(gini)
end


mat I_SE_I = . , . 
levelsof state, local(state)
	foreach s of local state{
	preserve
		keep if state=="`s'"
		drop if y_I<0
		bootstrap r(gini), reps(500) nodots : ineqdec0_se
		mat I_SE_I_`s' = r(table)[1,1] , r(table)[2,1]
		mat names I_SE_I_`s'= `s'
		mat I_SE_I = I_SE_I \ I_SE_I_`s'
	restore
}

* Bottom

cap program drop ineqdec0_se
	program ineqdec0_se, rclass
				qui ineqdec0 y_B [fw=round(wgt_B)]
				return scalar gini = r(gini)
end

mat I_SE_B = . , . 
levelsof state, local(state)
	foreach s of local state{
	preserve
		keep if state=="`s'"
		drop if y_I<0
		bootstrap r(gini), reps(500) nodots : ineqdec0_se
		mat I_SE_B_`s' = r(table)[1,1] , r(table)[2,1]
		mat I_SE_B = I_SE_B \ I_SE_B_`s'
		mat names I_SE_B_`s'= `s'
	restore
}

* Top
cap program drop ineqdec0_se
	program ineqdec0_se, rclass
				qui ineqdec0 y_T [fw=round(wgt_T)]
				return scalar gini = r(gini)
end

mat I_SE_T = . , . 
levelsof state, local(state)
	foreach s of local state{
	preserve
		keep if state=="`s'"
		drop if y_I<0
		bootstrap r(gini), reps(500) nodots : ineqdec0_se
		mat I_SE_T_`s' = r(table)[1,1] , r(table)[2,1]
		mat I_SE_T = I_SE_T \ I_SE_T_`s'
		mat names I_SE_T_`s'= `s'
	restore
}

* Bottom and Top
cap program drop ineqdec0_se
	program ineqdec0_se, rclass
				qui ineqdec0 y_BT [fw=round(wgt_BT)]
				return scalar gini = r(gini)
end

mat I_SE_BT = . , . 
levelsof state, local(state)
	foreach s of local state{
	preserve
		keep if state=="`s'"
		drop if y_I<0
		bootstrap r(gini), reps(500) nodots : ineqdec0_se
		mat I_SE_BT_`s' = r(table)[1,1] , r(table)[2,1]
		mat I_SE_BT = I_SE_BT \ I_SE_BT_`s'
		mat names I_SE_BT_`s'= `s'
	restore
}

mat I_SE = I_SE_I , I_SE_B , I_SE_T , I_SE_BT

matrix colnames I_SE = gini_I_ineq  se_I_ineq  gini_B_ineq  se_B_ineq  gini_T_ineq  se_T_ineq  gini_BT_ineq  se_BT_ineq
putexcel set "${pc}/Table.xlsx", sheet("M_US_ineqdec0", replace)  modify
putexcel A1 = matrix(I_SE), names
*/
*-------------------------------------------------------------------------------
* Import in Excel variables
*-------------------------------------------------------------------------------

import excel "${pc}/Table.xlsx", sheet("M_US_conindex") cellrange(A1:I5) firstrow clear

*-------------------------------------------------------------------------------
* T-Test
* https://nces.ed.gov/nationsreportcard/tdw/analysis/2004_2005/infer_compare2_overlap.aspx
*-------------------------------------------------------------------------------
rename A state
drop if gini_I==.

gen t_IB = abs(gini_B_con-gini_I_con) / (se_B_con^2 + (1-2*0.97)*se_I_con^2)^0.5
gen t_IT = abs(gini_T_con-gini_I_con) / (se_T_con^2 + (1-2*0.99)*se_I_con^2)^0.5
gen t_IBT = abs(gini_BT_con-gini_I_con) / (se_BT_con^2 + (1-2*0.96)*se_I_con^2)^0.5

gen delta_BI = gini_B_con-gini_I
gen delta_TI = gini_T_con-gini_I
gen delta_BTI = gini_BT_con-gini_I

order state gini_I_con se_I_con gini_B_con se_B_con delta_BI t_IB gini_T_con se_T_con delta_TI t_IT gini_BT_con se_BT_con delta_BTI t_IBT

export excel using "${pc}/Table.xlsx", sheet("Ginis_TTest", replace) firstrow(var) 
