{smcl}
{* *! version 0.9.4  25oct2026}{...}
{viewerjumpto "Description" "connectedness_postestimation##description"}{...}
{viewerjumpto "estat table" "connectedness_postestimation##table"}{...}
{viewerjumpto "estat group" "connectedness_postestimation##group"}{...}
{viewerjumpto "estat correlation" "connectedness_postestimation##correlation"}{...}
{viewerjumpto "estat bootstrap" "connectedness_postestimation##bootstrap"}{...}
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
{synopt:{cmd:estat bootstrap}}bootstrap standard errors and bands, event probabilities and
detection of significant rises (VAR models){p_end}
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


{marker bootstrap}{...}
{title:Syntax for estat bootstrap}

{p 8 16 2}
{cmd:estat} {cmdab:boot:strap} [{cmd:,} {it:options}]

{synoptset 22}{...}
{p2coldent:{it:options}}Description{p_end}
{synoptline}
{synopt:{opt r:eps(#)}}draws per round; default {cmd:reps(1000)}{p_end}
{synopt:{opt nobias}}one round of draws, without bias correction{p_end}
{synopt:{opt l:evel(#)}}level of the percentile intervals; default {cmd:level(95)}{p_end}
{synopt:{opt seed(#)}}random-number seed{p_end}
{synopt:{opt even:t(dates)}}event dates, for example {cmd:event(17oct1997 11sep2001)}, or values of the time variable{p_end}
{synopt:{opt lead:s(numlist)}}leads after the event window; default {cmd:leads(0 1 5 22)}{p_end}
{synopt:{opt ref:erence(#)}}compare with the mean of the {it:#} windows before the event; default 1{p_end}
{synopt:{opt thr:eshold(#)}}probability of a change of at least {it:#} percent{p_end}
{synopt:{opt dec:rease}}probability of a fall instead of a rise{p_end}
{synopt:{opt det:ect(#)}}probability, for every window, of a rise of at least {it:#} percent from the previous window{p_end}
{synopt:{opt st:at(measure)}}measure tested by {opt event()} and {opt detect()}: {cmd:tci} (default), {cmd:to} {it:name}, {cmd:from} {it:name}, {cmd:net} {it:name}, {cmd:pairwise} {it:name1 name2}, {cmd:npdc} {it:name1 name2}, {cmd:sam} or {cmd:sam} {cmd:to}|{cmd:from}|{cmd:net} {it:name}{p_end}
{synopt:{opt cut:off(#)}}probability marked as significant; default {cmd:cutoff(90)}{p_end}
{synopt:{opt every(#)}}bootstrap every {it:#}th window for the bands; default 1{p_end}
{synopt:{opt all}}bands for TO, FROM and NET as well as the TCI{p_end}
{synopt:{opt pair:wise}}bands for the pairwise connectedness table{p_end}
{synopt:{opt npdc}}bands for net pairwise directional connectedness{p_end}
{synopt:{opt sam}}bands for the spillover asymmetry measure (asymmetric results){p_end}
{synopt:{opt pwfr:ame(name)}}frame for the pairwise bands; default {it:frame}{cmd:_pw}{p_end}
{synopt:{opt notab:le}}do not display the pairwise and NPDC tables (static results){p_end}
{synopt:{opt fr:ame(name)}}frame for the bands; default {cmd:cn_boot}{p_end}
{synopt:{opt replace}}replace an existing frame{p_end}
{synopt:{opt nodots}}suppress the progress dots{p_end}
{synopt:{opt dig:its(#)}}decimals in the static table; default {cmd:digits(2)}{p_end}
{synoptline}

{pstd}
{cmd:estat bootstrap} gives bootstrap inference for the connectedness measures of a
VAR, as Greenwood-Nimmo, Kočenda and Nguyen (2024) propose. It requires results of
{cmd:connectedness} with {cmd:model(var)} and {cmd:approach(time)}, static or
rolling, with the generalized or orthogonalized decomposition, {opt normalize()} or
{opt clusters()}, and the estimation data in memory. Quantile, LAD, penalized,
TVP-VAR and DCC-GARCH models are not covered.

{pstd}
For each window, the VAR is estimated by OLS with an intercept and new samples are
generated recursively from it: the first {it:p} observations of the window are kept
and the residuals, recentred, are drawn with replacement. Each sample is
re-estimated and its connectedness computed; draws whose VAR has a companion root of
modulus one or more are discarded and redrawn. A window whose own VAR is not stable
uses the estimates of the most recent window with a stable VAR, as the authors do,
keeping its own initial values; the estimate then reported is that window's. A first round of {opt reps()} draws
estimates the bias of each measure, the mean draw minus the estimate, which
Greenwood-Nimmo, Kočenda and Nguyen find to be positive in every window; a second
round of {opt reps()} draws minus the bias gives the corrected draws, which centre on
the estimate. The command reports the estimate, the bias, the bias-corrected
estimate (estimate minus bias), the standard error of the corrected draws and their
percentile interval recentred on the bias-corrected estimate. Static results give a
table for TCI, TO, FROM and NET; rolling results give a frame with, for every
bootstrapped window, {cmd:tci}, {cmd:tci_bias}, {cmd:tci_bc}, {cmd:tci_se},
{cmd:tci_lo} and {cmd:tci_hi} (and the same for {cmd:to_}{it:i}, {cmd:from_}{it:i}
and {cmd:net_}{it:i} with {opt all}), {cmd:ndiscard}, the number of unstable
draws, and {cmd:stable}, 0 for windows that use the estimates of an earlier window. A full rolling bootstrap refits the VAR 2 x {opt reps()} times per window;
{opt every()} bootstraps fewer windows.

{pstd}
{opt pairwise} and {opt npdc} add the pairwise connectedness table, whose element
(i, j) is the spillover from variable j to variable i, and net pairwise directional
connectedness, whose element (i, j) is the net spillover from j to i (the spillover
from j to i minus that from i to j, as {cmd:NPDC} in ConnectednessApproach). They come
from the same draws at almost no extra cost, and the bias-corrected estimates keep
the identities of the measures: the rows of the pairwise table, own shares excluded,
sum to FROM, its columns to TO, and the columns of the NPDC matrix to NET. A
bias-corrected pairwise share can be negative when the bias exceeds the estimate.
Static results store k x k matrices (see {help connectedness_postestimation##results:stored results})
and display the bias-corrected pairwise table and the NPDC matrix, with a star where
the interval of a net spillover excludes zero; {opt notable} suppresses the display.
Rolling results store the bands in a second frame, {it:frame}{cmd:_pw}
({opt pwframe()}), with one observation per window and ordered pair: the time and
date variables, {cmd:from} and {cmd:to}, the transmitting and receiving variables,
and {cmd:pw}, {cmd:pw_bias}, {cmd:pw_bc}, {cmd:pw_se}, {cmd:pw_lo}, {cmd:pw_hi}
(spillover from {cmd:from} to {cmd:to}) and the same for {cmd:npdc} (net spillover
from {cmd:from} to {cmd:to}). With 19 variables a window has 342 ordered pairs.

{pstd}
{opt event()} computes, for each event, the probability (percent) that the TCI of
the window ending on the event date, and of the windows {opt leads()} later, exceeds
its level just before the event: the share of corrected draws above the mean
corrected draw of the window before the event (Greenwood-Nimmo, Kočenda and Nguyen
2024, eq. 5), or of the {opt reference()} windows before it. With
{opt threshold(#)}, the probability is that of a rise of at least {it:#} percent
(eq. 9); {opt decrease} gives the probability of a fall. Dates are read with the date
variable of the results (or the time variable if it holds daily dates); a date
without a window, such as a weekend, is moved to the next window (Greenwood-Nimmo,
Kočenda and Nguyen date a weekend event by the previous trading day, which can be
given instead). Only the windows needed are
bootstrapped. Probabilities of at least {opt cutoff()} percent, 90 in the paper, are
marked.

{pstd}
{opt detect(#)} bootstraps every window and stores in the frame {cmd:prob}, the
probability of a rise of at least {it:#} percent from the previous window (of any
rise with {cmd:detect(0)}), and {cmd:flag}, equal to 1 when it reaches
{opt cutoff()}, to locate changes at unknown dates.

{pstd}
{opt stat()} applies {opt event()} and {opt detect()} to another measure than the
TCI: {cmd:stat(to US)}, {cmd:stat(from US)} or {cmd:stat(net US)} for one variable,
{cmd:stat(pairwise US JPN)} for the spillover from US to JPN, and
{cmd:stat(npdc US JPN)} for the net spillover from US to JPN. The probabilities are
computed as for the TCI from the corrected draws of that measure. Since NET and NPDC
can be zero or negative, rises in percent, {opt threshold()} and {opt detect(#)}
with {it:#} > 0, are not allowed for them; {cmd:detect(0)} gives the probability of
any rise.

{pstd}
{opt sam} gives bootstrap inference for the spillover asymmetry measure of Baruník,
Kočenda and Vácha (2016), SAM = TCI+ - TCI-, the TCI of the positive system minus
that of the negative system, after {cmd:connectedness} with {opt asymmetric} or
{opt positive()} and {opt negative()}. The total, positive and negative systems are
bootstrapped jointly: every draw resamples the same rows of the residuals of the
three VARs, so that the draws keep the dependence between the systems, which share
their shocks, and each system is rebuilt from its own initial values. A draw is
discarded if any of the three VARs is unstable, and a window whose VAR is unstable
takes, system by system, the estimates of the most recent window with a stable VAR.
The bias of the SAM is the bias of the difference, the mean draw of TCI+ - TCI- minus
the estimate, and the bias correction, standard error and interval are those of the
other measures. {cmd:prob} is the probability (percent) that the measure is above
zero, the share of the corrected draws, recentred as the interval, above zero: the
interval at level {it:l} excludes zero when {cmd:prob} is above 50 + {it:l}/2 or below
50 - {it:l}/2. Static results display and store ({cmd:r(sam)}) the table of SAM and its
TO, FROM and NET (the differences of TO, FROM and NET between the positive and
negative systems), with {cmd:prob}, and the TCI of each system ({cmd:r(tci_parts)}).
Rolling results add to the frame {cmd:sam}, {cmd:sam_bias}, {cmd:sam_bc},
{cmd:sam_se}, {cmd:sam_lo}, {cmd:sam_hi} and {cmd:sam_prob} (the same for
{cmd:sam_to_}{it:i}, {cmd:sam_from_}{it:i} and {cmd:sam_net_}{it:i} with {opt all}),
and {cmd:tci_pos} and {cmd:tci_neg} with their bands. The measures of the total system
come from the same draws. {cmd:stat(sam)}, or {cmd:stat(sam net} {it:name}{cmd:)} and
so on, applies {opt event()} and {opt detect()} to the SAM; as for NET, only
{cmd:detect(0)} and no {opt threshold()} are allowed.

{pstd}
The joint bootstrap extends the procedure of Greenwood-Nimmo, Kočenda and Nguyen (2024)
to the SAM. Baruník, Kočenda and Vácha (2016) judge the SAM against bands simulated
from a model without asymmetry; {opt sam} instead resamples the residuals of the data
at hand, window by window. A Monte Carlo study under a symmetric data generating
process, where the SAM is zero, is in {cmd:docs/sam-coverage-v0.9.4.txt} of the
package repository: with 500 observations the intervals of the SAM cover zero about as
often as the intervals of the TCI of a correctly specified VAR cover the true TCI,
slightly less than their nominal level (0.86 and 0.90 for 90% and 95% intervals), and
the coverage improves with 1,000 observations (0.88 and 0.92).


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

{pstd}Greenwood-Nimmo, Kočenda and Nguyen (2024): bootstrap event probabilities{p_end}
{phang2}{cmd:. use cn_gkn2024, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness US-TUR, lags(2) horizon(10) window(202) datevar(date) replace}{p_end}
{phang2}{cmd:. estat bootstrap, event(17oct1997 11sep2001 09aug2007) leads(0 1 5 22) seed(1357)}{p_end}
{phang2}{cmd:. estat bootstrap, reps(200) every(20) seed(1357)}{p_end}
{phang2}{cmd:. frame cn_boot: twoway rarea tci_lo tci_hi date || line tci_bc date}{p_end}
{phang2}{cmd:. estat bootstrap, event(17oct1997 11sep2001 09aug2007) seed(1357) stat(npdc US JPN)}{p_end}
{phang2}{cmd:. estat bootstrap, reps(200) every(20) seed(1357) pairwise npdc replace}{p_end}
{phang2}{cmd:. frame cn_boot_pw: twoway (rarea npdc_lo npdc_hi date) (line npdc_bc date) if from == "US" & to == "JPN"}{p_end}

{pstd}Bootstrap bands for the spillover asymmetry measure of Baruník, Kočenda and Vácha (2016){p_end}
{phang2}{cmd:. use cn_g2020, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness EUR GBP CHF JPY in 1/1000, asymmetric lags(1) horizon(10) replace}{p_end}
{phang2}{cmd:. estat bootstrap, sam reps(500) seed(2094)}{p_end}
{phang2}{cmd:. connectedness EUR GBP CHF JPY, asymmetric lags(1) horizon(10) window(500) datevar(date) replace}{p_end}
{phang2}{cmd:. estat bootstrap, sam reps(100) every(100) seed(2094)}{p_end}
{phang2}{cmd:. frame cn_boot: twoway rarea sam_lo sam_hi date || line sam_bc date, yline(0)}{p_end}
{phang2}{cmd:. estat bootstrap, event(15sep2008) stat(sam) reps(500) seed(2094)}{p_end}


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

{pstd}
{cmd:estat bootstrap} stores the following in {cmd:r()}:

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(reps)}}draws per round{p_end}
{synopt:{cmd:r(level)}}level of the intervals{p_end}
{synopt:{cmd:r(nwindows)}}number of windows bootstrapped{p_end}
{synopt:{cmd:r(ndiscard)}}number of unstable draws discarded{p_end}
{synopt:{cmd:r(ncarried)}}windows with an unstable VAR that use the estimates of an earlier window{p_end}
{synopt:{cmd:r(reference)}}{opt reference()} ({opt event()}){p_end}
{synopt:{cmd:r(threshold)}}{opt threshold()} ({opt event()}){p_end}
{synopt:{cmd:r(nflag)}}windows flagged ({opt detect()}){p_end}
{synopt:{cmd:r(sam_prob)}}probability (percent) that the SAM is above zero ({opt sam}, static results){p_end}
{synopt:{cmd:r(nsam)}}windows whose interval of the SAM excludes zero ({opt sam}, rolling results){p_end}

{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(bias)}}{cmd:corrected} or {cmd:none}{p_end}
{synopt:{cmd:r(events)}}event dates ({opt event()}){p_end}
{synopt:{cmd:r(leads)}}leads ({opt event()}){p_end}
{synopt:{cmd:r(frame)}}frame holding the bands{p_end}
{synopt:{cmd:r(pwframe)}}frame holding the pairwise bands ({opt pairwise}, {opt npdc}){p_end}
{synopt:{cmd:r(stat)}}measure tested ({opt event()}, {opt detect()}){p_end}

{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}estimate, bias, bias-corrected estimate, standard error and interval of TCI, TO, FROM and NET (static results){p_end}
{synopt:{cmd:r(pw_b)}}pairwise table, estimate, k x k ({opt pairwise}, static results); also {cmd:r(pw_bias)}, {cmd:r(pw_bc)}, {cmd:r(pw_se)}, {cmd:r(pw_lo)} and {cmd:r(pw_hi)}{p_end}
{synopt:{cmd:r(npdc_b)}}NPDC matrix, estimate, k x k ({opt npdc}, static results); also {cmd:r(npdc_bias)}, {cmd:r(npdc_bc)}, {cmd:r(npdc_se)}, {cmd:r(npdc_lo)} and {cmd:r(npdc_hi)}{p_end}
{synopt:{cmd:r(sam)}}estimate, bias, bias-corrected estimate, standard error, interval and probability above zero of SAM and its TO, FROM and NET ({opt sam}, static results){p_end}
{synopt:{cmd:r(tci_parts)}}estimate, bias, bias-corrected estimate, standard error and interval of the TCI of the positive and negative systems ({opt sam}, static results){p_end}
{synopt:{cmd:r(prob)}}event probabilities (percent), events x leads ({opt event()}){p_end}


{marker references}{...}
{title:References}

{phang}
Baruník, J., E. Kočenda, and L. Vácha. 2016. Asymmetric connectedness on the U.S.
stock market: Bad and good volatility spillovers. {it:Journal of Financial Markets}
27: 55-78.

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
Greenwood-Nimmo, M., E. Kočenda, and V. H. Nguyen. 2024. Detecting statistically
significant changes in connectedness: A bootstrap-based technique.
{it:Economic Modelling} 140: 106843.

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
