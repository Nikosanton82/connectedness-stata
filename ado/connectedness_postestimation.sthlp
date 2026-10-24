{smcl}
{* *! version 0.8.2  24oct2026}{...}
{viewerjumpto "Description" "connectedness_postestimation##description"}{...}
{viewerjumpto "estat table" "connectedness_postestimation##table"}{...}
{viewerjumpto "estat group" "connectedness_postestimation##group"}{...}
{viewerjumpto "estat correlation" "connectedness_postestimation##correlation"}{...}
{viewerjumpto "Examples" "connectedness_postestimation##examples"}{...}
{viewerjumpto "Stored results" "connectedness_postestimation##results"}{...}
{viewerjumpto "References" "connectedness_postestimation##references"}{...}
{title:Title}

{phang}
{bf:connectedness postestimation} {hline 2} Postestimation tools for connectedness


{marker description}{...}
{title:Description}

{pstd}
The following postestimation commands are available after {helpb connectedness}:

{synoptset 20}{...}
{p2coldent:Command}Description{p_end}
{synoptline}
{synopt:{cmd:estat table}}redisplay the connectedness table{p_end}
{synopt:{cmd:estat group}}group measures: aggregated, conditional, inclusive, exclusive,
internal and external connectedness{p_end}
{synopt:{cmd:estat correlation}}conditional or partial correlations of each period{p_end}
{synopt:{helpb cngraph}}graph the time-varying measures{p_end}
{synopt:{helpb cnnetwork}}network plot of pairwise connectedness{p_end}
{synopt:{helpb estimates}}cataloging estimation results{p_end}
{synoptline}


{marker table}{...}
{title:Syntax for estat table}

{p 8 16 2}
{cmd:estat} {cmdab:tab:le} [{cmd:,} {opt dig:its(#)} {opt par:t(all|positive|negative)} {opt ban:d(#)} {opt notab:le}]

{pstd}
Typing {cmd:connectedness} without arguments has the same effect. {opt part()} shows
the tables of positive or negative changes after {cmd:connectedness, asymmetric};
{opt band(#)} shows the table of frequency band {it:#} after
{cmd:connectedness, approach(frequency)}, and {opt component()} the contemporaneous
or lagged table after {cmd:connectedness, approach(r2)}.


{marker group}{...}
{title:Syntax for estat group}

{p 8 16 2}
{cmd:estat} {cmdab:gro:up}{cmd:,} {opt t:ype(type)} {opt gr:oups(groups)}
[{opt fr:ame(name)} {opt replace} {opt dig:its(#)} {opt notab:le}]

{synoptset 16}{...}
{p2coldent:{it:type}}Description{p_end}
{synoptline}
{synopt:{cmd:aggregated}}connectedness between groups, each group treated as one unit
(Stenfors, Chatziantoniou and Gabauer 2022){p_end}
{synopt:{cmd:conditional}}connectedness within a subset of variables, renormalized as if
the subset were the whole system (Stenfors, Chatziantoniou and Gabauer 2022){p_end}
{synopt:{cmd:inclusive}}spillovers that involve the group, in either direction
(Chatziantoniou et al. 2023){p_end}
{synopt:{cmd:exclusive}}spillovers that do not involve the group
(Chatziantoniou et al. 2023){p_end}
{synopt:{cmd:internal}}spillovers within each group (Gabauer and Gupta 2018){p_end}
{synopt:{cmd:external}}spillovers between groups (Gabauer and Gupta 2018){p_end}
{synoptline}

{pstd}
{it:groups} lists the groups, each in double quotes and optionally preceded by a
name and a colon, for example
{cmd:groups("US: US_FPU US_MPU" "JP: JP_FPU JP_MPU")}. {cmd:aggregated},
{cmd:internal} and {cmd:external} require at least two groups;
{cmd:conditional}, {cmd:inclusive} and {cmd:exclusive} take one group, which need not be
quoted. A variable may belong to at most one group; variables in no group are left
as they are.

{pstd}
{cmd:estat group} works on the decomposition stored for every period in the results
frame, so {cmd:connectedness} must not have been run with {cmd:nofevd}, and it
requires the time-domain or model-free measures ({cmd:approach(time)} or
{cmd:approach(modelfree)}), as in the R package, with the default row
normalization ({cmd:normalize(row)}). For each
period it transforms the decomposition and recomputes TCI, corrected TCI, TO, FROM,
NET, NPT and NPDC; it displays the table of the averaged transformed decomposition
and stores the series in {opt frame()} (default {cmd:cn_group}), which
{helpb cngraph} and {helpb cnnetwork} accept through their {opt frame()} option.

{pstd}
The frame also holds, for {cmd:aggregated}, {cmd:tciext} and {cmd:ctciext}, the
quantity that ConnectednessApproach returns as TCI_ext, sum_g TO_g (k - k_g)/k, where
TO_g comes from the aggregated table, and its version with k - 1 in the denominator;
it is a sum over groups and is not bounded by 100; and, for {cmd:internal} and {cmd:external},
{cmd:gtci_}{it:g} and {cmd:cgtci_}{it:g}, the group TCI, sum of TO over the members of
group {it:g} divided by k or k - 1 (Gabauer and Gupta 2018).

{pstd}
{opt notable} suppresses the table; the header, the TCIs and the group TCIs are still
shown. {opt digits(#)} sets the number of decimals, default {cmd:digits(2)}.

{pstd}
Pairwise connectedness and influence indices are not reported for group measures. The R
package returns non-finite values for them in the time domain, so there is no reference
to validate against.


{marker correlation}{...}
{title:Syntax for estat correlation}

{p 8 16 2}
{cmd:estat} {cmdab:corr:elation} [{cmd:,} {opt par:tial} {opt fr:ame(name)} {opt replace}
{opt dig:its(#)} {opt notab:le}]

{pstd}
{cmd:estat correlation} turns the covariance matrices that {cmd:connectedness} stores
in the results frame, {cmd:cov_}{it:i}{cmd:_}{it:j}, into correlations, as R's
{cmd:ConditionalCorrelation()}, or with {opt partial} into partial correlations, as R's
{cmd:PartialCorrelations()}: with {bf:P} the inverse of the covariance matrix
{bf:H}, the partial correlation of {it:i} and {it:j} is
-p_ij/(p_ii p_jj)^(1/2), the correlation of the two variables given all the others.
The covariances are the residual covariance of the VAR, of each rolling window or of
each TVP-VAR period, or the DCC-GARCH conditional covariance, so the correlations are
the conditional correlations of the model of each period.

{pstd}
The series {cmd:cor_}{it:i}{cmd:_}{it:j} (or {cmd:pcor_}{it:i}{cmd:_}{it:j}), {it:i} < {it:j},
are saved in {opt frame()} (default {cmd:cn_corr}) with the time variable and the date
variable of the results frame, so {helpb tsline} graphs them directly. The command
displays the average correlation matrix over the periods, with ones on the diagonal
(R's {cmd:PartialCorrelations()} puts -1 there). It requires the covariance
matrices, so {cmd:connectedness} must not have been run with {opt nocovariance},
nor with {cmd:approach(r2)} or {cmd:approach(modelfree)} on the data, which estimate
no covariance matrix; after {cmd:model(dccgarch) approach(r2)} the conditional
covariances are available. {opt digits(#)} sets the number of decimals, default
{cmd:digits(3)}, and {opt notable} suppresses the table.


{marker examples}{...}
{title:Examples}

{pstd}Gabauer and Gupta (2018): internal and external connectedness{p_end}
{phang2}{cmd:. use cn_gg2018, clear}{p_end}
{phang2}{cmd:. tsset month}{p_end}
{phang2}{cmd:. connectedness US_FPU-JP_CPU, model(tvpvar) window(200) datevar(date) replace}{p_end}
{phang2}{cmd:. estat group, type(internal) groups("US: US_FPU US_MPU US_TPU US_CPU" "JP: JP_FPU JP_MPU JP_TPU JP_CPU") replace}{p_end}
{phang2}{cmd:. cngraph gtci, frame(cn_group)}{p_end}
{phang2}{cmd:. estat group, type(external) groups("US: US_FPU US_MPU US_TPU US_CPU" "JP: JP_FPU JP_MPU JP_TPU JP_CPU") replace}{p_end}

{pstd}Stenfors, Chatziantoniou and Gabauer (2022): conditional and aggregated connectedness{p_end}
{phang2}{cmd:. use cn_cgs2022, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness X2YCHF-X10YSEK, model(tvpvar) window(250) corrected datevar(date) replace}{p_end}
{phang2}{cmd:. estat group, type(conditional) groups(X2YCHF X2YDKK X2YEUR X2YGBP X2YNOK X2YSEK) replace}{p_end}
{phang2}{cmd:. estat group, type(aggregated) replace groups("CHF: X2YCHF X5YCHF X10YCHF" "DKK: X2YDKK X5YDKK X10YDKK" "EUR: X2YEUR X5YEUR X10YEUR" "GBP: X2YGBP X5YGBP X10YGBP" "NOK: X2YNOK X5YNOK X10YNOK" "SEK: X2YSEK X5YSEK X10YSEK")}{p_end}
{phang2}{cmd:. cngraph net, frame(cn_group)}{p_end}

{pstd}Gabauer (2020): DCC-GARCH conditional and partial correlations{p_end}
{phang2}{cmd:. use cn_g2020, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness EUR GBP CHF JPY, model(dccgarch) datevar(date) replace}{p_end}
{phang2}{cmd:. estat correlation}{p_end}
{phang2}{cmd:. frame cn_corr: tsline cor_1_2, title("EUR and GBP, DCC-GARCH conditional correlation")}{p_end}
{phang2}{cmd:. estat correlation, partial frame(cn_pcorr) replace}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:estat group} stores the following in {cmd:r()}:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(TCI)}}TCI of the averaged transformed decomposition{p_end}
{synopt:{cmd:r(cTCI)}}corrected TCI{p_end}
{synopt:{cmd:r(n)}}number of units in the table{p_end}

{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(type)}}group measure{p_end}
{synopt:{cmd:r(names)}}names of the units in the table{p_end}
{synopt:{cmd:r(groups)}}group names{p_end}
{synopt:{cmd:r(frame)}}frame holding the series{p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}connectedness table, (n+4) x (n+1){p_end}
{synopt:{cmd:r(gtci)}}averages of tciext and ctciext ({cmd:aggregated}) or of
gtci_g and cgtci_g ({cmd:internal}, {cmd:external}){p_end}

{pstd}
{cmd:estat correlation} stores the following in {cmd:r()}:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(T)}}number of periods{p_end}
{synopt:{cmd:r(k)}}number of variables{p_end}

{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(type)}}{cmd:conditional} or {cmd:partial}{p_end}
{synopt:{cmd:r(frame)}}frame holding the series{p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(mean)}}average correlation matrix, k x k{p_end}


{marker references}{...}
{title:References}

{phang}
Chatziantoniou, I., A. H. Elsayed, D. Gabauer, and G. Gozgor. 2023. Oil price shocks and
exchange rate dynamics: Evidence from decomposed and partial connectedness measures for
oil importing and exporting economies. {it:Energy Economics} 120: 106627.

{phang}
Gabauer, D. 2020. Volatility impulse response analysis for DCC-GARCH models: The role
of volatility transmission mechanisms. {it:Journal of Forecasting} 39: 788-796.

{phang}
Gabauer, D., and R. Gupta. 2018. On the transmission mechanism of country-specific and
international economic uncertainty spillovers: Evidence from a TVP-VAR connectedness
decomposition approach. {it:Economics Letters} 171: 63-71.

{phang}
Stenfors, A., I. Chatziantoniou, and D. Gabauer. 2022. Independent policy, dependent
outcomes: A game of cross-country dominoes across European yield curves.
{it:Journal of International Financial Markets, Institutions and Money} 81: 101658.


{marker author}{...}
{title:Author}

{pstd}
Nikolaos Antonakakis{break}
Department of Accounting, Economics and Finance{break}
University of Nicosia, UNIC Athens{break}
Athens, Greece{break}
{browse "mailto:antonakakis.n@unic.ac.cy":antonakakis.n@unic.ac.cy}

{pstd}
The numerical routines port the {it:ConnectednessApproach} R package by David
Gabauer, distributed under the GNU General Public License version 3; this package
is distributed under the same license.


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume ??, number ??: {browse "https://doi.org/":st0001}{p_end}

{p 7 14 2}
Help:  {helpb connectedness}, {helpb cngraph}, {helpb cnnetwork}, {helpb cnportfolio}{p_end}
