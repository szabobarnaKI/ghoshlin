// program to perform competing risk adjustment on mi set data, 
// according to Geskus method using inverse probability of censoring weights

program define ghoshlin
	version 16.1
	syntax [anything] [if] [in], frame(string asis) failure(string) competing(string) [maxevents(integer 0) keep(varlist fv) matchonly(varlist fv)]
	local id="`_dta[st_id]'"
	local tstop="`_dta[st_bt]'"
	// if maxevents is not specified (<=0), change it to more than max observed events (number of rows is always more than max potential events)
	if `maxevents'<=0 local maxevents=_N
	
	capture u_mi_assert_set
	if _rc==0 {
		local mi_set="mi"
	}
	else local mi_set=""
	
	foreach i of local keep {
		assert "`i'"!="tstart"
		assert "`i'"!="tstop"
		assert "`i'"!="event"
		assert "`i'"!="weight_c"
	}

	if "`id'"=="" {
		di "You have to stset the data with an id() variable before using ghoshlin."
		error 198
	}
	
	de_factorize `keep', local(keep)
	de_factorize `matchonly', local(matchonly)

	gettoken framename replace : frame, parse(", ")
	local replace=subinstr("`replace'"," ","",.)
	local replace=subinstr("`replace'",",","",.)

	if "`replace'"=="replace" {
	    capture frame drop `framename'
	}
	quietly frame
	local currentframe="`r(currentframe)'"
	quietly frame copy `r(currentframe)' `framename'
	quietly frame `framename': {
		if "`in/'"!="" {
			keep `in'
		}
		if "`if/'"!="" {
			keep `if'
		}
		if "`mi_set'"=="mi" {
			mi extract 0, clear
		}
		drop if _st!=1
		
		// only observations that are stset and that should be used are kept
		
		generate double _old_t=_t
		generate double _old_t0=_t0
		generate _old_fail=(`failure')
		tempfile keepfile
		quietly save `keepfile'
		
		// store list of time indexes for each event, terminal or not, for later use
		tempname eventtimes
		frame copy `framename' `eventtimes'

		// create a single row dataset for the last time point, with minimal data kept
		tempname terminal_event time_max failure_event reachmaxevent eventcount maxfup
		
		stgen `eventcount'=count0(`failure') // this contains the number of failures until the actual time interval´s end
		drop if (`eventcount'>`maxevents') | (`eventcount'==`maxevents' & !(`failure'))
		// now only observations up to maxevents number of failures are kept
		stgen `reachmaxevent'=ever(`eventcount'==`maxevents')
		stgen `terminal_event'=ever(`competing')
		stgen double `time_max'=max(_t)
		summarize `time_max', meanonly
		scalar `maxfup'=r(max)
		generate `failure_event'=.
		replace `failure_event'=(`failure') if `time_max'==_t

		keep if `time_max'==_t
		keep `id' `time_max' `reachmaxevent' `terminal_event' `failure_event'

		// process the event time list
		frame `eventtimes' : {
			rename _t `time_max'
			keep `time_max' //_time_max
			// this frame now contains the list of event times (any terminal event or failure event)
			set obs `=_N+2'
			replace `time_max'=`maxfup' if _n==_N
			replace `time_max'=0 if _n==_N-1
			duplicates drop
			sort `time_max'
			generate cens_at_risk0=. // number of subjects at risk for being censored at this time point
			generate censored=. // number of subjects censored until next time index, among those at risk here
			generate double inv_cens_prob=. // variable for the calculated censoring probability for the time interval ending at time_max
			tempname time_index time_index_1
			forvalues i=1/`=_N' {
				scalar `time_index'=`time_max'[`i']
				scalar `time_index_1'=`time_max'[`=`i'+1']
				frame `framename': count if `time_max'>=`time_index'
				replace cens_at_risk0=r(N) if _n==`i'
				frame `framename': count if `time_max'>=`time_index' & `time_max'<`time_index_1' & `terminal_event'==0 & `failure_event'==0
				replace censored=r(N) if _n==`i'
			}
			drop if censored==0 & `time_max'!=0
			//**********************************************
			replace inv_cens_prob=1-censored/cens_at_risk0
			//**********************************************
			generate double ln_ipc=log(inv_cens_prob)
			
			generate index=_n
			tempname time_index_count
			scalar `time_index_count'=_N
		}
		
		// reduce the dataset to subjects who need to remain in the risk set after a terminal event
		keep if `reachmaxevent'==0 & `terminal_event'==1
		generate expand_num=0
		tempname subj_tstop
		forvalues i=1/`=_N' { // loop over subjects
			scalar `subj_tstop'=`time_max'[`i']
			// calculate the number of added time intervals after a terminal event, for each subject
			frame `eventtimes': count if `time_max'>`subj_tstop'
			replace expand_num=r(N) if _n==`i' // expand_num now contains the number of time intervals until the end of the longest follow-up
		}
		
		// expand and calculate weights
		expand =expand_num
		sort `id'

		generate next_t=.
		generate next_index=.
		replace next_index=`time_index_count'-expand_num+1 if `id'!=`id'[_n-1]
		replace next_index=next_index[_n-1]+1 if next_index==.

		frlink m:1 next_index, frame(`eventtimes' index)
		frget ln_ipc, from(`eventtimes')
		frget tstop=`time_max', from(`eventtimes')

		generate cum_log_ipc=.
		replace cum_log_ipc=ln_ipc if `id'!=`id'[_n-1]
		replace cum_log_ipc=cum_log_ipc[_n-1]+ln_ipc if cum_log_ipc==.

		generate weight_c=exp(cum_log_ipc)
		drop if weight_c==.
		
		generate double tstart=.
		replace tstart=`time_max' if `id'!=`id'[_n-1]
		replace tstart=tstop[_n-1] if tstart==.
		
		// cleanup
		keep `id' tstart tstop weight_c
		order `id' weight_c tstart tstop 
		generate event=0
	}	
	
	// add all sub-observations from the original dataset, until the terminal events
	quietly frame `framename': {
		drop if tstart==0
		append using `keepfile', keep(`id' _old_t0 _old_t _old_fail)
		
		replace tstart=_old_t0 if tstart==.
		replace tstop=_old_t if tstop==.
		replace event=1 if _old_fail==1
		
		drop if tstart==. | tstop==.
		drop _old_t0 _old_t _old_fail
		replace event=0 if event==.
		replace weight_c=1 if weight_c==.
	}

	// re-add multiple imputations if needed
	if "`mi_set'"=="mi" {
		// add M imputations
		local imputations=`_dta[_mi_M]'
		quietly frame `framename': {
			mi set flong
			mi set M=`imputations'
		}
	}
	
	// link the new dataset to the original using frlink and add "keep" variables
	preserve // original frame
	quietly drop if `tstop'==.
	quietly frame `framename': {
		if "`keep'"!="" {
			local mi_matching=cond("`mi_set'"=="mi","_mi_m ","")
//			frlink m:1 `id' `mi_matching' tstop, frame(`currentframe' `id' `mi_matching' `tstop') generate(lnk_`currentframe')
			frlink m:1 `id' `mi_matching' tstop, frame(`currentframe' `id' `mi_matching' _t) generate(lnk_`currentframe')
			sort `mi_matching' `id' tstart
			// add all covariates in option "keep", matching each imputation
			frget `keep', from(lnk_`currentframe')
			drop lnk_`currentframe'
		}
	}
	restore

	// fill in gaps in the frlink-ed data, unless the variable is listed in matchonly()
	quietly frame `framename': {
		foreach i of local keep {
		    local tofill=1
		    foreach j of local matchonly {
			    if "`i'"=="`j'" {
				    local tofill=0
					continue, break
				}
			}
			if `tofill'==1 {
				capture confirm numeric variable `i'
				if _rc==0 {
					replace `i'=`i'[_n-1] if `id'==`id'[_n-1] & `i'==.
				}
				else {
					replace `i'=`i'[_n-1] if `id'==`id'[_n-1] & `i'==""
				}
			}
		}
	}
	
	frame `framename': noisily `mi_set' stset tstop [pw=weight_c], failure(event) enter(tstart)
	di as result "Dataset successfully converted to a Fine & Gray or Ghosh-Lin model into frame `framename'."
	di as result "Database is {it:`mi_set' stset} for the outcome defined by the {it:failure()} option."
	di as result _n "To stset the data for the failure event, use the following code:"
	di as result "{it:`mi_set' stset tstop [pw=weight_c], failure(event==1) enter(tstart)}"
	di as result `"Do not forget to use {it:"vce(cluster `id')"} or {it:"vce(robust)"} in your regression models!"'
end

program define de_factorize
	syntax [anything(name=factorvars)], local(string asis)
	local factorvars=subinstr("`factorvars'","#"," ",.) // split interactions to separate variables
	local x:word count `factorvars'
	local output=""
	forvalues i=1/`x' {
		local s : word `i' of `factorvars'
		if strpos("`s'",".")!=0 { // factor variable notation found
			local s=substr("`s'",strpos("`s'",".")+1,strlen("`s'")-strpos("`s'","."))
			local output="`output' `s'"
		}
		else local output="`output' `s'"
	}
	c_local `local'="`output'"
end
