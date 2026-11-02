{smcl}
{* *! version 1.3  01nov2026}{...}
{viewerjumpto "Syntax" "ghoshlin##syntax"}{...}
{viewerjumpto "Description" "ghoshlin##description"}{...}
{viewerjumpto "Model considerations" "ghoshlin##considerations"}{...}
{viewerjumpto "Weights" "ghoshlin##weights"}{...}
{viewerjumpto "Options" "ghoshlin##options"}{...}
{viewerjumpto "Remarks" "ghoshlin##remarks"}{...}
{viewerjumpto "Examples" "ghoshlin##examples"}{...}
{title:Title}

{phang}
{cmd:ghoshlin} {hline 2} Competing risk adjustment for time to event data, using inverse probability of censoring weights (as described by Ghosh-Lin 2002).


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:ghoshlin}
{ifin}
{cmd:,} frame({it:newframename}[,replace]) failure({it:string}) competing({it:string}) [maxevents({it:integer}) keep({it:varlist}) matchonly({it:varlist}) {opt noexp:andfw}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt frame}}The name of a new frame for the transformed dataset. {p_end}
{synopt:{opt failure}}Expression identifying failures. It needs to be enclosed in quotes or compound double quotes.{p_end}
{synopt:{opt competing}}Expression that defines competing risk (or terminal) events. The terminal events may overlap with failure events.
It needs to be enclosed in quotes or compound double quotes.{p_end}
{synopt:{opt maxevents}}Maximum number of events per subject to analyze. If not defined, all failure events without upper limit will be analyzed.{p_end}
{synopt:{opt keep}}Variables kept during dataset transformation.{p_end}
{synopt:{opt matchonly}}Variables kept during dataset transformation, but without filling the gaps.{p_end}
{syntab:Weights}
{synopt:{opt noexp:andfw}}Do not expand frequency weights into independent subjects (see {help ghoshlin##weights:Weights}).{p_end}
{synoptline}
{p 4 6 2}Weights are specified in {help stset}, not in {cmd:ghoshlin}; see {help ghoshlin##weights:Weights}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ghoshlin} performs competing risk adjustment for time to event data, somewhat similarly to the user-written {cmd:stcrprep} 
command. The transformed dataset is stored in a new frame defined by the {opt frame()} option.
{p_end}

{pstd}
{cmd:ghoshlin} supports right censored data and left truncation, handled through stset. Gaps are not checked for and are 
ignored. {cmd: ghoshlin} will still run if time gaps exist in the dataset, making the indirect assumption that no terminal/competing event can occur during gap 
times. This is true if terminal events occur only at the latest follow-up time (e.g. death), but may be a false assumption if follow-up is possible after a competing 
event (e.g. preventive breast removal prevents the occurrence of breast cancer, but is not lethal, thus follow-up after a potentially uncaptured breast removal is 
possible). It is the user's responsibility to consider whether this assumption can safely be accepted or not.
{p_end}

{pstd}
The data has to be stset before calling {cmd: ghoshlin}, with both target and competing events as 
failures. This command recognizes whether the data is multiply imputed, and acts accordingly. Be careful though, as only imputed covariates will be kept, but not imputed events or follow-up times.
{p_end}

{pstd}
For successful conversion, you need to stset your data with the id() option before calling {cmd:ghoshlin}. Note that in the current version of Stata, you can not stset with time varying weights, thus, after {cmd:ghoshlin} 
the data is stset without the use of the {opt id} option. Therefore, it is important to use cluster variance estimation ({opt vce(cluster patient_id)}) or robust variance ({opt vce(robust)}) when fitting stcox-based 
models to the transformed dataset. The use of {opt vce(cluster id)} is preferred as it provides the most similar results compared to {help stcrreg} (no difference in standard errors up to 4 digits).{p_end}

{pstd}
The command automatically calls {help stset} or {help mi stset} after dataset transformation, there is no absolute need to do it manually. If you want to change st settings, you need to re-call {help stset} or {help mi stset} 
after {cmd:ghoshlin}.
{p_end}


{marker considerations}{...}
{title:Model considerations}

{pstd}
The command {cmd:ghoshlin} transforms the dataset so that for each subject that experienced a terminal (competing) event, it will include weighted observations after the competing event to modify the risk set, as if the subject would 
still be at risk despite the terminal 
event. If used on time to first event data (or by limiting the dataset to {opt maxevents(1)}), using {cmd: ghoshlin} and subsequently fitting {cmd: stcox} results in a Fine & Gray model (as described by 
{browse "https://pubmed.ncbi.nlm.nih.gov/20377575/" :Geskus}). When used with repeated events data in a setup according to the Andersen-Gill model 
(or Lin-Wei-Yang-Ying model which is the same but analyzed with the {opt vce(robust)} or {opt vce(cluster id)} option), the transformed dataset can be used to fit a 
{browse "https://www3.stat.sinica.edu.tw/statistica/j12n3/j12n31/j12n31.htm" :Ghosh-Lin model} using the IPCW method.
{p_end}

{pstd}
Compared to {cmd: stcrreg}, {cmd: ghoshlin} provides substantial speed improvement for Fine & Gray analysis. Directly after launching Stata, using the Stata example dataset for hypoxia (as described in {help stcrreg}), 
the Stata built-in model takes 2.17sec on the author´s test computer, while the time for dataset conversion using {cmd: ghoshlin} and subsequent analysis using {cmd: stcox} takes only 
0.07sec. The time improvement is even larger for datasets with multiple imputation. The command results has been validated against Stata's built-in command {cmd:stcrreg} for time to first event, 
and against the R implementation of the Ghosh-Lin model using the HF-ACTION dataset published in the R package "mets".
{p_end}

{pstd}
{ul:Time-varying covariates:} The command extends the at-risk time of subjects beyond the terminal events, until the estimated time of non-informative 
censoring. The command itself cannot handle time-varying data beyond the terminal event on its own. Covariates that are independent from the subject (so-called exogenous variables such as season of the year, information on the president of the 
country at a specific time period) can be measured and included appropriately, and need to be added manually after running {cmd:ghoshlin} (see the {opt matchonly} 
option). For endogenous variables (those measured on the subject, such as blood pressure or the number of children), values observed at the terminal event may have had changed subsequently if the subject did not experience the terminal 
event. Therefore, the last-value-carry-forward method used as default, and that works for time-non-varying covariates, gives potentially false results, and estimates are difficult to 
interpret. Thus, the use of endogenous time-varying covariates is discouraged.
{p_end}


{marker weights}{...}
{title:Weights}

{pstd}
{cmd:ghoshlin} uses the weights declared to {help stset} ({cmd:fweight}s and {cmd:pweight}s). Weights are not specified in the {cmd:ghoshlin} command itself. The weight 
of a subject must be constant over the whole follow-up of that subject; otherwise {cmd:ghoshlin} stops with an error.
{p_end}

{pstd}
The stset weights enter the Kaplan-Meier estimator of the censoring distribution as weighted counts (number at risk and number censored). For frequency weights this is 
identical to replicating each subject as many times as its weight; for probability weights it gives the weighted Kaplan-Meier estimator. In the transformed dataset, 
{it:weight_c} is the product of the subject's stset weight and the censoring redistribution weight (on the observed records, {it:weight_c} equals the stset weight). The 
transformed dataset is always stset with {cmd:[pw=weight_c]}, because the product is generally not an integer.
{p_end}

{pstd}
{ul:Frequency weights:} by default, each record of a frequency weighted subject is replicated as many times as the frequency weight, so that every copy becomes an 
independent subject. A new identifier variable {it:repid} is created for the copies, and {it:weight_c} then contains only the censoring redistribution weight. Use 
{opt vce(cluster repid)} in subsequent {cmd:stcox} models in order to obtain correct point estimates and standard errors. Without 
expansion ({opt noexpandfw}), clustering on the original {it:id} treats all copies of a subject as one cluster; point estimates are the same (exactly with 
{opt breslow} ties), but the standard errors are not those of a frequency weighted analysis. Frequency weights must be integers for the expansion. The expansion is 
not done for {help mi} set data (a note is displayed), and it is ignored for data that are not stset with frequency weights.
{p_end}

{pstd}
{ul:Probability weights:} no expansion is done, {it:weight_c} carries the product of the probability weight and the redistribution weight, and {opt vce(cluster id)} 
should be used as for unweighted data.
{p_end}

{pstd}
Note that the standard errors from {cmd:stcox} on the transformed dataset (weighted or not) treat the censoring weights as known, i.e. they do not account for the 
estimation of the censoring distribution.
{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt frame} defines the name of a new frame for the transformed dataset. Use the suboption {it:replace} to optionally replace an existing frame.

{phang}
{opt failure} defines failures. E.g {it:failure(death_type==1)} where {it:death_type} is a variable taking the value of 1 if that subject experienced the event at interest (e.g CV death), death_type==2 for the competing absorbing event 
(e.g non-CV death) and death_type==0 for no event (censoring).

{phang}
{opt competing} defines competing risk events, using the same syntax as the {opt failure} option.

{phang}
{opt maxevents} defines the maximum number of repeated events per subject. Subjects experiencing the maximum analyzed failure count are not kept in the risk set after the last event, regardless of competing/terminal or not.

{phang}
{opt keep} Variables kept during dataset transformation. By default, gaps will be filled in by last value carry forward method.

{phang}
{opt matchonly} Variables kept during dataset transformation, that should not be filled in. If matchonly() contains variables not listed in keep(), those will not be kept.

{dlgtab:Weights}

{phang}
{opt noexpandfw} switches off the default expansion of frequency weights declared to {help stset} into independent subjects (see {help ghoshlin##weights:Weights}). 
The frequency weight is then kept as a multiplier in {it:weight_c}. The default (expansion) can also be requested explicitly with the option {opt expandfw}; 
{opt expandfw} and {opt noexpandfw} may not be combined. Both options are ignored when the data are not stset with frequency weights.


{marker remarks}{...}
{title:Remarks on factor variable notation}

{pstd}
In the {opt keep()} and {opt matchonly()} options, both regular variable names and factor variable notations are 
allowed.{cmd: ghoshlin} interprets them by removing the factor variable notation, thus you may use both forms without any concern.
{p_end}


{title:AI disclosure}

{pstd}
This command was originally developed by a human. Claude AI with human supervision was used for speed-optimization, to correct the handling of left truncation and 
ties, and to add support for stset weights. The final version of {cmd: ghoshlin} has then been extensively tested and validated by a human.


{marker examples}{...}
{title:Examples}
{pstd}

	{title:Example 1:} {it:Use ghoshlin to fit a Ghosh-Lin model.}
	{it:The {bf:events} variable in the example dataset codes outcomes and takes values of .: censored, 1: CV death, 2: Non-CV death, 3: HF hospitalization.}
	{it:The primary endpoint is total HF hospitalizations and CV death. CV and Non-CV death are terminal events.}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days, id(ID) failure(event == 1, 2, 3) exit(time .)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("inlist(event,1,3)") competing("inlist(event,1,2)") keep(age sex diabetes) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

	{title:Example 2:} {it:Use ghoshlin to fit a Fine-Gray model.}
	{it:The primary endpoint is CV death. Non-CV death is the competing event. HF hospitalizations are omitted.}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days, id(ID) failure(event == 1, 2) exit(time .)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("event==1") competing("event==2") keep(age sex diabetes) maxevents(1) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

	{title:Example 3:} {it:Use ghoshlin to fit a Fine-Gray model. Equivalent to example 2.}
	{it:The primary endpoint is CV death. Non-CV death is the competing event. HF hospitalizations are omitted.}
	{it:Note that if the maximum number of events is 1, it does not make any difference if the failure event is also included as terminal (competing) event.}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days, id(ID) failure(event == 1, 2) exit(time .)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("inlist(event,1)") competing("inlist(event,1,2)") keep(age sex diabetes) maxevents(1) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

	{title:Example 4:} {it:Re-stset after ghoshlin in case stset has been changed.}
{pstd}

{pmore2}{cmd:. frame transformed: stset tstop [pw=weight_c], failure(event) enter(tstart)}{p_end}

	{title:Example 5:} {it:Ghosh-Lin model with frequency weights.}
	{it:The variable {bf:fw} holds an integer frequency weight that is constant within each subject (here a hypothetical variable).}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days [fw=fw], id(ID) failure(event == 1, 2, 3) exit(time .)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("inlist(event,1,3)") competing("inlist(event,1,2)") keep(age sex diabetes) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster repid) nolog}{p_end}

	{title:Example 6:} {it:Same as example 5, without expanding the frequency weights.}
{pstd}

{pmore2}{cmd:. ghoshlin, failure("inlist(event,1,3)") competing("inlist(event,1,2)") keep(age sex diabetes) frame(transformed,replace) noexpandfw}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

{hline}

{marker author}{...}
{title:Author}

{pstd}
Barna Szabó-Söderberg{break}
{it:Karolinska Institutet}{break}
{ul:barna.szabo@ki.se}
{p_end}

