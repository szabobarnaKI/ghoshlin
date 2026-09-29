// program to perform competing risk adjustment on mi set data, 
// according to Ghosh and Lin's IPCW method (2002), which coincides with
// Geskus method (2011) for Fine-Gray models using time to first events only

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
		
		// create a single row dataset for the last time point, with minimal data kept
		tempname terminal_event time_max time_min failure_event reachmaxevent eventcount

		stgen `eventcount'=count0(`failure') // this contains the number of failures until the actual time interval´s end
		drop if (`eventcount'>`maxevents') | (`eventcount'==`maxevents' & !(`failure'))
		// now only observations up to maxevents number of failures are kept
		stgen `reachmaxevent'=ever(`eventcount'==`maxevents')
		stgen `terminal_event'=ever(`competing')
		stgen double `time_max'=max(_t)
		stgen double `time_min'=min(_t0)
		generate `failure_event'=.
		replace `failure_event'=(`failure') if `time_max'==_t

		keep if `time_max'==_t
		keep `id' `time_max' `time_min' `reachmaxevent' `terminal_event' `failure_event'

		// build the censoring-probability time grid and the inverse-probability-
		// of-censoring weights entirely in Mata (ghoshlin_grid_weights, backed
		// by ghoshlin_calc_ipcw at the bottom of this file): a Kaplan-Meier-type
		// estimator of the censoring distribution that groups subjects by
		// exact, double-precision equality of `time_max' (avoiding the
		// float-precision copy that was the original source of incorrect
		// weights when ties occurred), accounts for left truncation via
		// `time_min', and follows the "event < censor < entry" tie order used
		// by the reference Geskus-redistribution implementation. The grid
		// (columns: time_max, at-risk count, censored count,
		// 1-censored/at-risk, log-IPC increment) is kept entirely inside
		// Mata's own workspace (an external Mata matrix, not a Stata
		// matrix and not a temporary frame): with a large number of
		// distinct censoring times, a genuine Stata matrix would be capped
		// by matsize (798 in Stata/IC, 11,000 in SE/MP, regardless of how
		// matsize is set), while a Mata matrix is bounded only by memory.
		ghoshlin_grid_weights `time_max' `time_min' `terminal_event' `failure_event'

		// reduce the dataset to subjects who need to remain in the risk set after a terminal event
		keep if `reachmaxevent'==0 & `terminal_event'==1

		// number of synthetic post-terminal-event rows each subject needs,
		// computed directly against the Mata-resident grid in a single Mata
		// pass (ghoshlin_expand_count, backed by ghoshlin_calc_expand_num)
		// instead of one frame switch and risk-set scan per subject.
		tempvar expand_num
		ghoshlin_expand_count `time_max', gen(`expand_num')

		// expand and fill in the synthetic rows
		expand =`expand_num'
		sort `id'

		tempvar newgrp
		generate byte `newgrp'=(`id'!=`id'[_n-1])

		// fill tstart/tstop/weight_c for every expanded row in a single Mata
		// pass (ghoshlin_fill_expanded, backed by ghoshlin_store_fill),
		// replacing the previous per-row frlink/frget lookup against a
		// temporary frame and the row-by-row cumulative-sum replace logic.
		// id itself never enters Mata, so it may be numeric or string.
		ghoshlin_fill_expanded `time_max' `expand_num' `newgrp', tstart(tstart) tstop(tstop) weight(weight_c)
		capture mata drop ghoshlin__grid // free the Mata-resident grid now that it has been fully consumed
		drop if weight_c==.

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

// ---------------------------------------------------------------------------
// ghoshlin_grid_weights
//
// Computes the censoring-probability time grid used to extend a subject's
// risk set after a terminal/competing event: a Kaplan-Meier-type estimator
// of the censoring distribution (see ghoshlin_calc_ipcw at the bottom of
// this file for the tie-handling and left-truncation details).
//
// Must be called with the per-subject frame (one row per id, containing
// each subject's overall follow-up window and event status) as the current
// frame. It requires four existing double-precision variables in that
// frame: the subject's last follow-up time ("time_max"), the subject's
// stset entry/left-truncation time ("time_min"), an indicator for ever
// having a terminal/competing event, and the failure-event indicator at
// that last time.
//
// Stores the grid (columns: time_max, at-risk count, censored count,
// 1-censored/at-risk, log-IPC increment) into the external Mata matrix
// ghoshlin__grid (via ghoshlin_store_grid at the bottom of this file) -
// deliberately never as a Stata matrix (st_matrix() results are capped by
// matsize - 798 in Stata/IC, 11,000 in SE/MP - which a large number of
// distinct censoring times could exceed) and never as a temporary frame.
// ghoshlin_expand_count and ghoshlin_fill_expanded read it back directly
// from Mata; ghoshlin__grid is dropped again once ghoshlin_fill_expanded
// has consumed it (see the main ghoshlin program).
// ---------------------------------------------------------------------------
program define ghoshlin_grid_weights
	syntax varlist(min=4 max=4)
	tokenize `varlist'
	local timevar  `1'
	local entryvar `2'
	local termvar  `3'
	local failvar  `4'

	mata: ghoshlin_store_grid(st_data(.,"`timevar'"), st_data(.,"`entryvar'"), st_data(.,"`termvar'"), st_data(.,"`failvar'"))
end

// ---------------------------------------------------------------------------
// ghoshlin_expand_count
//
// For each subject in the current frame, counts how many grid points in
// the Mata-resident grid (ghoshlin__grid, set by ghoshlin_grid_weights)
// lie strictly after that subject's own exit time - i.e. how many synthetic
// risk-set-extension rows expand() should create for them. Done as a single
// Mata pass over the grid (ghoshlin_calc_expand_num at the bottom of this
// file), rather than one frame switch and risk-set scan per subject.
// ---------------------------------------------------------------------------
program define ghoshlin_expand_count
	syntax varlist(min=1 max=1), gen(string)
	local timevar `varlist'

	quietly generate double `gen'=.
	mata: st_store(., "`gen'", ghoshlin_calc_expand_num_ext(st_data(.,"`timevar'")))
end

// ---------------------------------------------------------------------------
// ghoshlin_fill_expanded
//
// Fills tstart/tstop/weight_c for every synthetic post-terminal-event row
// created by expand(). Must be called after expand() and sort(id), with
// three existing variables: each row's own subject-level exit time
// (duplicated across the expansion by expand()), that subject's
// ghoshlin_expand_count() result (also duplicated), and a 0/1 flag marking
// the first row of each subject's block (e.g. id!=id[_n-1]) - this avoids
// ever comparing id itself inside Mata, so id may be numeric or string.
//
// Computes, in a single Mata pass (ghoshlin_store_fill at the bottom of
// this file) against the Mata-resident grid (ghoshlin__grid), each row's
// position in the grid from its rank within its subject's block, replacing
// the previous per-row frlink/frget lookup against a temporary frame and
// the row-by-row cumulative-sum replace logic.
// ---------------------------------------------------------------------------
program define ghoshlin_fill_expanded
	syntax varlist(min=3 max=3), tstart(string) tstop(string) weight(string)
	tokenize `varlist'
	local timevar   `1'
	local expandvar `2'
	local newgrpvar `3'

	quietly generate double `tstart'=.
	quietly generate double `tstop'=.
	quietly generate double `weight'=.
	mata: ghoshlin_store_fill_ext(st_data(.,"`timevar'"), st_data(.,"`expandvar'"), st_data(.,"`newgrpvar'"), "`tstart'", "`tstop'", "`weight'")
end

mata:
// Kaplan-Meier estimator of the censoring distribution, evaluated at every
// distinct subject-level exit time ("time_max"), with left truncation: a
// subject only counts toward the risk set from its own stset entry time
// ("time_min") onward, matching the follow-up window actually declared to
// stset (enter()/id()), rather than assuming everyone is at risk from t=0.
// Returns a k x 5 matrix with columns (time_max, number at risk, number
// censored, 1 - censored/at risk, log of that quantity), one row per
// distinct time at which someone was purely censored, plus boundary rows at
// t=0 and t=maxfup.
//
// Ties are handled by construction, following the same convention as the
// reference "Geskus redistribution" implementation for this class of model
// (R's survival::finegray, which cites Geskus 2011 - the method this file's
// own header already names): at a shared timestamp, a subject's own
// failure/terminal event is treated as occurring just before that instant,
// so it is removed from the censoring risk set there, while censoring at
// that same instant is treated as occurring after any such event, so
// censored subjects still count in both the risk set and the numerator at
// that time, and a subject entering exactly at that instant is not yet
// counted in the risk set there ("event < censor < entry", per finegray's
// own documentation of this tie order). Subjects sharing an identical
// time_max value are grouped and their at-risk/censored counts computed
// together in one comparison, rather than through hand-built interval
// boundaries that depended on prior deduplication of a value copied at a
// different numeric precision (the source of the original tie bug).
real matrix ghoshlin_calc_ipcw(real colvector t, real colvector entry,
                                real colvector term, real colvector fail,
                                real scalar maxfup)
{
	real colvector cens, ut, natrisk, ncens, invp, lnipc, keepidx, entered
	real scalar    k, i

	cens = (term :== 0) :& (fail :== 0) // subjects with neither a terminal nor a failure event at time_max are purely censored

	ut = uniqrows(t) // sorted, unique subject-level exit times
	if (rows(ut)==0 | ut[1] != 0)     ut = 0 \ ut
	if (ut[rows(ut)] != maxfup)       ut = ut \ maxfup
	ut = uniqrows(ut) // re-dedup in case 0 or maxfup already occurred

	k = rows(ut)
	natrisk = J(k,1,.)
	ncens   = J(k,1,.)
	invp    = J(k,1,.)
	lnipc   = J(k,1,.)

	for (i=1; i<=k; i++) {
		// at risk of censoring at ut[i]: entered strictly before ut[i] (an
		// entry exactly at ut[i] has not happened yet, per "entry" last in
		// the tie order), and either still under observation past ut[i], or
		// exiting exactly at ut[i] via pure censoring (a same-time
		// failure/terminal event has already left the risk set, per
		// "event" first in the tie order)
		entered    = (entry :< ut[i])
		natrisk[i] = sum(entered :& ((t :> ut[i]) :| ((t :== ut[i]) :& cens)))
		ncens[i]   = sum(entered :& (t :== ut[i]) :& cens) // number purely censored exactly at ut[i]
		if (natrisk[i] > 0) {
			invp[i]  = 1 - ncens[i]/natrisk[i]
			lnipc[i] = ln(invp[i])
		}
		else {
			// nobody has (strictly) entered the risk set yet at this
			// boundary time - always true at the forced t=0 row (entry < 0
			// is never true), and possibly elsewhere under left truncation;
			// there is no one to censor, so the weight does not change
			invp[i]  = 1
			lnipc[i] = 0
		}
	}

	// keep only rows where the weight actually changes, plus the two
	// boundary rows (t=0 anchors cumulation, t=maxfup bounds the risk-set
	// expansion after a terminal event)
	keepidx = selectindex((ncens :> 0) :| (ut :== 0) :| (ut :== maxfup))

	return ((ut[keepidx], natrisk[keepidx], ncens[keepidx], invp[keepidx], lnipc[keepidx]))
}

// Computes the grid via ghoshlin_calc_ipcw (maxfup taken as max(t), so it
// never needs to be passed through a syntax-parsed Stata option) and stores
// it into the external Mata matrix ghoshlin__grid, where it stays entirely
// inside Mata's own workspace rather than becoming a Stata matrix (whose
// size ghoshlin_calc_ipcw's own row count - one per distinct censoring
// time - could exceed matsize for a dataset with many distinct event
// times) or a temporary frame. "external" makes ghoshlin__grid persist
// across separate mata: blocks within the same Stata session, so
// ghoshlin_calc_expand_num_ext and ghoshlin_store_fill_ext (below) can read
// back whatever this call last stored.
void ghoshlin_store_grid(real colvector t, real colvector entry,
                          real colvector term, real colvector fail)
{
	external real matrix ghoshlin__grid

	ghoshlin__grid = ghoshlin_calc_ipcw(t, entry, term, fail, max(t))
}

// Counts, for each subject's own exit time in `t', how many grid points in
// `ut' (sorted ascending) lie strictly after it - the number of synthetic
// risk-set-extension rows that subject needs. A plain Mata loop over
// subjects, each iteration a vectorized comparison against the (typically
// much shorter) grid vector - no frame switching and no per-subject Stata
// command dispatch, unlike the loop this replaces.
real colvector ghoshlin_calc_expand_num(real colvector t, real colvector ut)
{
	real scalar k, i, m
	real colvector en

	k = rows(ut)
	m = rows(t)
	en = J(m,1,.)
	for (i=1; i<=m; i++) {
		en[i] = k - sum(ut :<= t[i])
	}
	return(en)
}

// Reads back the grid time column from ghoshlin__grid (set by
// ghoshlin_store_grid) and delegates to ghoshlin_calc_expand_num, so the
// Stata-side ghoshlin_expand_count subroutine never has to touch the grid
// directly.
real colvector ghoshlin_calc_expand_num_ext(real colvector t)
{
	external real matrix ghoshlin__grid

	return(ghoshlin_calc_expand_num(t, ghoshlin__grid[.,1]))
}

// Fills tstart/tstop/weight_c for every row produced by expand() +
// sort(id). `t' is each row's own subject-level exit time (constant within
// a subject's block), `en' that subject's expand count (also constant
// within the block), `newgrp' a 0/1 flag marking the first row of each
// block. Each subject's block of rows corresponds, in order, to the grid
// points strictly after their own exit time; this walks that
// correspondence with a running rank-within-block and a running cumulative
// log-IPC sum, replacing the previous per-row frlink/frget lookup and
// cum_log_ipc replace logic.
void ghoshlin_store_fill(real colvector t, real colvector en, real colvector newgrp,
                          real colvector ut, real colvector lnipc,
                          string scalar tstartname, string scalar tstopname,
                          string scalar weightname)
{
	real scalar n, k, i, rank, gidx, csum
	real colvector tstart_out, tstop_out, weight_out

	n = rows(t)
	k = rows(ut)
	tstart_out = J(n,1,.)
	tstop_out  = J(n,1,.)
	weight_out = J(n,1,.)

	rank = 0
	csum = 0
	for (i=1; i<=n; i++) {
		if (newgrp[i]) rank = 1
		else           rank = rank + 1

		gidx = k - en[i] + rank // this row's position in the grid

		if (rank==1) {
			tstart_out[i] = t[i]
			csum = lnipc[gidx]
		}
		else {
			tstart_out[i] = ut[gidx-1]
			csum = csum + lnipc[gidx]
		}
		tstop_out[i]  = ut[gidx]
		weight_out[i] = exp(csum)
	}

	st_store(., tstartname, tstart_out)
	st_store(., tstopname,  tstop_out)
	st_store(., weightname, weight_out)
}

// Reads back the grid's time and log-IPC columns from ghoshlin__grid (set
// by ghoshlin_store_grid) and delegates to ghoshlin_store_fill, so the
// Stata-side ghoshlin_fill_expanded subroutine never has to touch the grid
// directly.
void ghoshlin_store_fill_ext(real colvector t, real colvector en, real colvector newgrp,
                              string scalar tstartname, string scalar tstopname,
                              string scalar weightname)
{
	external real matrix ghoshlin__grid

	ghoshlin_store_fill(t, en, newgrp, ghoshlin__grid[.,1], ghoshlin__grid[.,5], tstartname, tstopname, weightname)
}
end
