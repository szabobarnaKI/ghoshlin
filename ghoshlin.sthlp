{smcl}
{* *! version 1.0  18aug2025}{...}
{viewerjumpto "Syntax" "ghoshlin##syntax"}{...}
{viewerjumpto "Description" "ghoshlin##description"}{...}
{viewerjumpto "Model considerations" "ghoshlin##considerations"}{...}
{viewerjumpto "Options" "ghoshlin##options"}{...}
{viewerjumpto "Remarks" "ghoshlin##remarks"}{...}
{viewerjumpto "Examples" "ghoshlin##examples"}{...}
{title:Title}

{phang}
{cmd:ghoshlin} {hline 2} Competing risk adjustment for time to event data, using inverse probability of censoring weights (as described by Geskus).


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:ghoshlin}
{ifin}
{cmd:,} frame({it:newframename}[,replace]) failure({it:string}) competing({it:string}) [maxevents({it:integer}) keep({it:varlist}) matchonly({it:varlist})]

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

{marker description}{...}
{title:Description}

{pstd}
{cmd:ghoshlin} performs competing risk adjustment for time to event data, somewhat similarly to the user-written {cmd:stcrprep} command.
The transformed dataset is stored in a new frame defined by the {opt frame()} option.
{p_end}

{pstd}
{cmd:ghoshlin} supports only right censored data without gaps, and without left truncation.
The data has to be stset before calling {cmd: ghoshlin}, with both target and competing events as failures.
This command recognizes whether the data is multiply imputed, and acts accordingly. Be careful though, as only imputed covariates will be kept, but not imputed events.
{p_end}

{pstd}
For successful conversion, you need to stset your data with the id() option before calling {cmd:ghoshlin}. Note that in the current version of Stata, you can not stset with time varying weights, thus, after {cmd:ghoshlin} 
the data is stset without the use of the {opt id} option.
Therefore, it is important to use cluster variance estimation ({opt vce(cluster patient_id)}) or robust variance ({opt vce(robust)}) when fitting stcox-based models to the transformed dataset. 
The use of {opt vce(cluster id)} is preferred as it provides the most similar results compared to {help stcrreg} (no difference in standard errors up to 4 digits).{p_end}

{pstd}
The command automatically calls {help stset} or {help mi stset} after dataset transformation, there is no absolute need to do it manually. If you want to change st settings, you need to re-call {help stset} or {help mi stset} 
after {cmd:ghoshlin}.
{p_end}


{marker considerations}{...}
{title:Model considerations}
{pstd}
The command {cmd:ghoshlin} transforms the dataset so that for each subject that experienced a terminal (competing) event, it will include weighted observations after the competing event to modify the risk set, as if the subject would 
still be at risk despite the terminal event.
If used on time to first event data (or by limiting the dataset to {opt maxevents(1)}), a Fine & Gray model is fit using {cmd: ghoshlin} and thereafter {cmd: stcox}. When used with repeated events data in a setup according to the Andersen-Gill model 
(or Lin-Wei-Yang-Ying model which is the same but analyzed with the {opt vce(robust)} or {opt vce(cluster id)} option), the transformed dataset can be used to fit a Ghosh-Lin model.
{p_end}

{pstd}
Compared to {cmd: stcrreg}, {cmd: ghoshlin} provides substantial speed improvement for Fine & Gray analysis. Using the Stata example dataset for hypoxia (as described in {cmd: stcrreg}), 
the Stata built-in model takes 11.68sec on the author´s test computer, while the time for dataset conversion using {cmd: ghoshlin} and subsequent analysis using {cmd: stcox} takes only 0.51sec.
The time improvement is even larger for multiple imputation data. The command results has been validated against Stata's built-in command {cmd:stcrreg} for time to first event, 
and against the R implementation of the Ghosh-Lin model using the HF-ACTION dataset published in the R package "mets".
{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt frame} defines the name of a new frame for the transformed dataset. Use the suboption {it:replace} to optionally replace an existing frame.

{phang}
{opt failure} defines failures. E.g {it:failure(death_type==1)} where {it:death_type} is a variable taking the value of 1 if that subject experienced the event at interest (e.g CV death), death_type==2 for the competing absorbing event (e.g non-CV death) and death_type==0 for no event (censoring).

{phang}
{opt competing} defines competing risk events, using the same syntax as the {opt failure} option.

{phang}
{opt maxevents} defines the maximum number of repeated events per subject. Subjects experiencing the maximum analyzed failure count are not kept in the risk set after the last event, regardless of competing/terminal or not.

{phang}
{opt keep} Variables kept during dataset transformation. By default, gaps will be filled in by last value carry forward method.

{phang}
{opt matchonly} Variables kept during dataset transformation, that should not be filled in. If matchonly() contains variables not listed in keep(), those will not be kept.

{marker remarks}{...}
{title:Remarks on factor variable notation}
{pstd}
In the {opt keep()} and {opt matchonly()} options, both regular variable names and factor variable notations are allowed. {cmd: ghoshlin} interprets them by removing the factor variable notation, thus you may use both forms without any concern.
{p_end}



{marker examples}{...}
{title:Examples}
{pstd}

	{title:Example 1:} {it:Use ghoshlin to fit a Ghosh-Lin model.}
	{it:The {bf:events} variable in the example dataset codes outcomes and takes values of .: censored, 1: CV death, 2: Non-CV death, 3: HF hospitalization.}
	{it:The primary endpoint is total HF hospitalizations and CV death. CV and Non-CV death are terminal events.}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days, id(ID) failure(event == 1, 2, 3)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("inlist(event,1,3)") competing("inlist(event,1,2)") keep(age sex diabetes) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

	{title:Example 2:} {it:Use ghoshlin to fit a Fine-Gray model.}
	{it:The primary endpoint is CV death. Non-CV death is the competing event. HF hospitalizations are omitted.}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days, id(ID) failure(event == 1, 2)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("event==1") competing("event==2") keep(age sex diabetes) maxevents(1) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

	{title:Example 3:} {it:Use ghoshlin to fit a Fine-Gray model. Equivalent to example 2.}
	{it:The primary endpoint is CV death. Non-CV death is the competing event. HF hospitalizations are omitted.}
	{it:Note that if the maximum number of events is 1, it does not make any difference if the failure event is also included as terminal (competing) event.}
{pstd}

{pmore2}{cmd:. sysuse ghoshlin_example}{p_end}
{pmore2}{cmd:. stset fup_days, id(ID) failure(event == 1, 2)}{p_end}
{pmore2}{cmd:. ghoshlin, failure("inlist(event,1)") competing("inlist(event,1,2)") keep(age sex diabetes) maxevents(1) frame(transformed,replace)}{p_end}
{pmore2}{cmd:. frame transformed: stcox c.age i.sex i.diabetes, vce(cluster ID) nolog}{p_end}

	{title:Example 4:} {it:Re-stset after ghoshlin in case stset has been changed.}
{pstd}

{pmore2}{cmd:. frame transformed: stset tstop [pw=weight_c], failure(event) enter(tstart)}{p_end}

{hline}

{marker author}{...}
{title:Author}

{pstd}
Barna Szabó-Söderberg{break}
{it:Karolinska Institutet}{break}
{ul:barna.szabo@ki.se}
{p_end}

