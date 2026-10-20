{smcl}
{* *! version 0.3.0  20oct2026}{...}
{viewerjumpto "Syntax" "connectedness##syntax"}{...}
{viewerjumpto "Description" "connectedness##description"}{...}
{viewerjumpto "Options" "connectedness##options"}{...}
{viewerjumpto "Remarks" "connectedness##remarks"}{...}
{viewerjumpto "Examples" "connectedness##examples"}{...}
{viewerjumpto "Stored results" "connectedness##results"}{...}
{viewerjumpto "Methods and formulas" "connectedness##methods"}{...}
{viewerjumpto "References" "connectedness##references"}{...}
{viewerjumpto "Author" "connectedness##author"}{...}
{title:Title}

{phang}
{bf:connectedness} {hline 2} VAR and TVP-VAR connectedness measures in the time and frequency domains


{marker syntax}{...}
{title:Syntax}

{p 8 22 2}
{cmd:connectedness} {varlist} {ifin} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt m:odel(model)}}{cmd:var} (default) or {cmd:tvpvar}{p_end}
{synopt:{opt l:ags(#)}}VAR lag order; default is {cmd:lags(1)}{p_end}
{synopt:{opt w:indow(#)}}rolling-window size ({cmd:var}) or prior training sample ({cmd:tvpvar}); default 0{p_end}

{syntab:TVP-VAR}
{synopt:{opt kappa1(#)}}forgetting factor for the coefficients; default {cmd:0.99}{p_end}
{synopt:{opt kappa2(#)}}decay factor for the error covariance; default {cmd:0.99}{p_end}
{synopt:{opt pr:ior(prior)}}{cmd:bayes} (default), {cmd:minnesota} or {cmd:uninformative}{p_end}
{synopt:{opt gam:ma(#)}}Minnesota prior tightness; default {cmd:0.01}{p_end}

{syntab:Connectedness}
{synopt:{opt app:roach(approach)}}{cmd:time} (default), {cmd:frequency}, {cmd:joint} or {cmd:extjoint}{p_end}
{synopt:{opt h:orizon(#)}}forecast horizon {it:H}; default {cmd:horizon(10)}{p_end}
{synopt:{opt gen:eralized}}generalized FEVD; the default{p_end}
{synopt:{opt orth:ogonalized}}Cholesky (orthogonalized) FEVD{p_end}
{synopt:{opt corr:ected}}report the corrected total connectedness index{p_end}
{synopt:{opt asym:metric}}also estimate the systems of positive and negative changes{p_end}

{syntab:Frequency}
{synopt:{opt band:s(numlist)}}band limits in periods, for example {cmd:bands(5 20)}{p_end}
{synopt:{opt part:ition(numlist)}}band limits in radians, in decreasing order{p_end}
{synopt:{opt scen:ario(scenario)}}{cmd:abs} (default) or {cmd:wth}{p_end}

{syntab:Results}
{synopt:{opt fr:ame(name)}}frame for the time-varying measures; default {cmd:cn_results}{p_end}
{synopt:{opt replace}}replace an existing results frame{p_end}
{synopt:{opt nofevd}}do not store the full FEVD for every period{p_end}
{synopt:{opt date:var(varname)}}copy a date variable into the results frame for graphs{p_end}

{syntab:Reporting}
{synopt:{opt notab:le}}suppress the connectedness table but show the header and TCI{p_end}
{synopt:{opt dig:its(#)}}decimals in the table; default {cmd:digits(2)}{p_end}
{synoptline}
{p 4 6 2}
You must {helpb tsset} your data before using {cmd:connectedness}; the estimation
sample must be contiguous in time. {it:varlist} may contain time-series operators.{p_end}
{p 4 6 2}
See {helpb connectedness postestimation} for {cmd:estat table} and {cmd:estat group}
(aggregated, conditional, inclusive, exclusive, internal and external connectedness),
and {helpb cngraph} and {helpb cnnetwork} for graphs.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:connectedness} estimates the connectedness measures of Diebold and Yilmaz
(2009, 2012) from a vector autoregression (VAR), either over the full sample or over
rolling windows, and the refined measures of Antonakakis, Chatziantoniou and Gabauer
(2020) from a time-varying parameter VAR (TVP-VAR) estimated with a Kalman filter
and forgetting factors. It reports the connectedness table and stores the
time-varying total (TCI), corrected total, directional (TO, FROM, NET), net pairwise
directional (NPDC), pairwise connectedness (PCI) and influence measures in a
separate frame.

{pstd}
The command is a port of the {it:ConnectednessApproach} R package (Gabauer 2025)
and reproduces its results; see {it:Remarks}.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt model(model)} selects the estimator. {cmd:var} fits a VAR({it:p}) by
ordinary least squares with an intercept, over the whole sample when
{cmd:window(0)} and over rolling windows otherwise. {cmd:tvpvar} fits the TVP-VAR of
Antonakakis, Chatziantoniou and Gabauer (2020), producing measures for every period
of the sample.

{phang}
{opt lags(#)} sets the VAR lag order {it:p}.

{phang}
{opt window(#)} has two meanings, following the R package. With {cmd:model(var)},
each rolling regression uses {it:#} - {it:p} observations, so a window of 200
with one lag uses 199 rows of data. With {cmd:model(tvpvar)} and {cmd:prior(bayes)},
the prior is an OLS VAR on the first {it:#} - {it:p} observations;
{cmd:window(0)} uses the whole sample.

{dlgtab:TVP-VAR}

{phang}
{opt kappa1(#)} and {opt kappa2(#)} are the forgetting and decay factors of Koop and
Korobilis (2013). Both must lie strictly between 0 and 1. Antonakakis,
Chatziantoniou and Gabauer (2020) use 0.99 and 0.96.

{phang}
{opt prior(prior)} sets the prior for the initial coefficients. {cmd:bayes} uses
an OLS VAR on a training sample (see {opt window()}); {cmd:minnesota} uses zero means
with variance {it:gamma}/{it:l}^2 for lag {it:l}; {cmd:uninformative} uses zero means
with variance 4.

{phang}
{opt gamma(#)} sets the Minnesota prior tightness.

{dlgtab:Connectedness}

{phang}
{opt horizon(#)} sets the forecast horizon {it:H}; the decomposition sums horizons
0 to {it:H}-1.

{phang}
{opt generalized} and {opt orthogonalized} choose the generalized forecast error
variance decomposition (Koop, Pesaran and Potter 1996; Pesaran and Shin 1998),
which is the default, or the Cholesky decomposition used by Diebold and Yilmaz (2009),
which depends on the order of {it:varlist}.

{phang}
{opt approach(approach)} selects the connectedness framework applied to the
estimated model. {cmd:time} (the default) gives the time-domain measures of
Diebold and Yilmaz (2012) and their refinements. {cmd:frequency} decomposes them
into frequency bands as in Baruník and Křehlík (2018); with
{cmd:model(tvpvar)} this is the TVP-VAR frequency approach of Chatziantoniou,
Gabauer and Gupta (2023). {cmd:joint} gives the joint spillover index of Lastrapes
and Wiesen (2021), and {cmd:extjoint} the extended joint connectedness of Balcilar,
Gabauer and Umar (2021). {cmd:joint} and {cmd:extjoint} use the generalized
decomposition and do not allow {opt orthogonalized} or {opt corrected}.

{phang}
{opt corrected} makes the corrected TCI of Chatziantoniou and Gabauer (2021) the
headline index in the table and in {helpb cngraph}. Both versions are always stored.

{dlgtab:Frequency}

{phang}
{opt bands(numlist)} sets the frequency bands by their limits in periods (cycle
lengths), each greater than 1. {cmd:bands(5 20)} gives three bands, of 1 to 5,
5 to 20 and more than 20 periods. The band limits in radians are pi/{it:#}, and the
upper limit of the first band is pi + 0.00001 so that the highest frequency is
included, as in the replications of the R package.

{phang}
{opt partition(numlist)} sets the band limits directly in radians, in decreasing
order, as the {cmd:partition} argument of the R package. Each band covers the
frequencies in [{it:lower}, {it:upper}). Without {opt bands()} or {opt partition()},
the R default {cmd:partition(}pi pi/2 0{cmd:)} is used; note that it leaves out the
frequency pi itself.

{phang}
{opt scenario(scenario)} selects how the band tables are scaled. {cmd:abs} (the
default) keeps the absolute contributions, so that the band measures add up to the
total; {cmd:wth} rescales each band's table to sum to k, giving within-band
measures.

{pstd}
The frequency decomposition uses the discrete Fourier transform of the
moving-average matrices for horizons 0 to {it:H}, so {opt horizon()} sets the
frequency resolution; replications of the R package use {cmd:horizon(100)}.

{dlgtab:Results}

{phang}
{opt frame(name)} names the frame that receives the time-varying measures; see
{it:Stored results}. The current frame is left unchanged.

{phang}
{opt replace} allows {opt frame()} to overwrite an existing frame.

{phang}
{opt asymmetric} estimates the model three times, as in Adekoya et al. (2022): on the
series in {it:varlist}, on their positive parts (negative values set to zero) and on
their negative parts (positive values set to zero). The variables are typically
returns or percentage changes. The results go to {opt frame()},
{opt frame()}{cmd:_pos} and {opt frame()}{cmd:_neg}; {cmd:connectedness, part(positive)}
or {cmd:part(negative)} redisplays the other tables, and {helpb cngraph} overlays the
three sets of measures automatically.

{phang}
{opt nofevd} stops {cmd:connectedness} from storing the full decomposition,
{cmd:fevd_}{it:i}{cmd:_}{it:j}, for every period in the results frame. The
decomposition is needed by {cmd:estat group}; omit it only when k is large and the
number of variables would exceed {help maxvar}. ({cmd:savefevd}, from version 0.1,
is still accepted and has no effect.)

{phang}
{opt datevar(varname)} copies a date variable into the results frame; {helpb cngraph}
then uses it as the horizontal axis. This is useful when the data are
{cmd:tsset} by an observation index, as with daily trading data.


{marker remarks}{...}
{title:Remarks}

{pstd}
Typing {cmd:connectedness} without arguments redisplays the results. Options
{opt digits(#)}, {opt notable}, {opt part(positive|negative)} (after
{opt asymmetric}) and {opt band(#)} (after {cmd:approach(frequency)}, to show the
table of band {it:#}) are allowed.

{pstd}
The results reproduce those of {it:ConnectednessApproach} 1.0.4 to numerical
precision, including its conventions: the forecast error variance decomposition sums
horizons 0 to {it:H}-1; rolling windows contain {it:window} - {it:p} observations;
the Bayes prior is estimated on {it:window} - 2{it:p} regression observations; the
TVP-VAR uses the updated coefficients only when their companion matrix is stable;
otherwise it keeps the previous period's coefficients for the decomposition and
shrinks the state to 0.99 times its previous value. The certification script
in the package's {cmd:test} folder checks these results against the R package.

{pstd}
Data from 18 published applications are installed with the package as
{cmd:cn_}{it:name}{cmd:.dta}, for example {cmd:cn_dy2012.dta} and
{cmd:cn_acg2020.dta}. They are ancillary files: type {cmd:net get connectedness} to copy
them to the current directory, then load them with {cmd:use}. Each is {cmd:tsset}-ready through the observation index
{cmd:t} (monthly data also carry {cmd:month}) and has a {cmd:date} variable.


{marker examples}{...}
{title:Examples}

{pstd}Diebold and Yilmaz (2012): generalized connectedness, static and rolling{p_end}
{phang2}{cmd:. use cn_dy2012, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10)}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) window(200) datevar(date) replace}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}

{pstd}Antonakakis, Chatziantoniou and Gabauer (2020): TVP-VAR connectedness{p_end}
{phang2}{cmd:. use cn_acg2020, clear}{p_end}
{phang2}{cmd:. tsset month}{p_end}
{phang2}{cmd:. connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) replace}{p_end}
{phang2}{cmd:. cngraph net}{p_end}
{phang2}{cmd:. cngraph npdc}{p_end}

{pstd}Chatziantoniou and Gabauer (2021): corrected TCI{p_end}
{phang2}{cmd:. use cn_cg2021, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness Austria-Spain, model(tvpvar) window(200) corrected datevar(date) replace}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}
{phang2}{cmd:. cnnetwork}{p_end}

{pstd}Chatziantoniou, Gabauer and Gupta (2023): TVP-VAR frequency connectedness{p_end}
{phang2}{cmd:. use cn_cgg2022, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness Europe-OPEC, model(tvpvar) window(200) horizon(100) approach(frequency) bands(5) datevar(date) replace}{p_end}
{phang2}{cmd:. connectedness, band(1)}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}
{phang2}{cmd:. cngraph net, band(2)}{p_end}

{pstd}Lastrapes and Wiesen (2021): joint spillover index{p_end}
{phang2}{cmd:. use cn_lw2021, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness energy-utilities, lags(2) horizon(30) approach(joint)}{p_end}

{pstd}Balcilar, Gabauer and Umar (2021): extended joint connectedness{p_end}
{phang2}{cmd:. use cn_bgu2021, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness CrudeOil-Wheat, model(tvpvar) horizon(20) prior(minnesota) gamma(0.1) approach(extjoint) datevar(date) replace}{p_end}
{phang2}{cmd:. cngraph net}{p_end}

{pstd}Adekoya et al. (2022): asymmetric connectedness of positive and negative changes{p_end}
{phang2}{cmd:. use cn_aaacgo2022, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. foreach v of varlist Brent-UTI {c -(}}{p_end}
{phang2}{cmd:.     generate double r_`v' = 100*(`v' - L.`v')/L.`v'}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. connectedness r_*, model(tvpvar) prior(minnesota) gamma(0.1) asymmetric datevar(date) replace}{p_end}
{phang2}{cmd:. connectedness, part(negative)}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:connectedness} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations used{p_end}
{synopt:{cmd:e(T0)}}number of periods with measures{p_end}
{synopt:{cmd:e(k)}}number of variables{p_end}
{synopt:{cmd:e(lags)}}lag order{p_end}
{synopt:{cmd:e(horizon)}}forecast horizon{p_end}
{synopt:{cmd:e(window)}}window size (0 if none){p_end}
{synopt:{cmd:e(TCI)}}total connectedness index of the averaged decomposition{p_end}
{synopt:{cmd:e(cTCI)}}corrected total connectedness index{p_end}
{synopt:{cmd:e(kappa1)}, {cmd:e(kappa2)}}forgetting and decay factors ({cmd:tvpvar}){p_end}
{synopt:{cmd:e(gamma)}}Minnesota prior tightness ({cmd:prior(minnesota)}){p_end}
{synopt:{cmd:e(TCI_pos)}, {cmd:e(cTCI_pos)}}TCI of positive changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(TCI_neg)}, {cmd:e(cTCI_neg)}}TCI of negative changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(nbands)}}number of frequency bands ({cmd:frequency}){p_end}
{synopt:{cmd:e(TCI_f}{it:b}{cmd:)}, {cmd:e(cTCI_f}{it:b}{cmd:)}}TCI and corrected TCI of band {it:b} ({cmd:frequency}){p_end}

{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:connectedness}{p_end}
{synopt:{cmd:e(model)}}{cmd:var} or {cmd:tvpvar}{p_end}
{synopt:{cmd:e(approach)}}{cmd:time}, {cmd:frequency}, {cmd:joint} or {cmd:extjoint}{p_end}
{synopt:{cmd:e(bands)}}band labels in periods, for example {cmd:1-5 5-Inf} ({cmd:frequency}){p_end}
{synopt:{cmd:e(partition)}}band limits in radians ({cmd:frequency}){p_end}
{synopt:{cmd:e(scenario)}}{cmd:abs} or {cmd:wth} ({cmd:frequency}){p_end}
{synopt:{cmd:e(fevdtype)}}{cmd:generalized} or {cmd:orthogonalized}{p_end}
{synopt:{cmd:e(corrected)}}{cmd:corrected} or empty{p_end}
{synopt:{cmd:e(prior)}}prior ({cmd:tvpvar}){p_end}
{synopt:{cmd:e(names)}}variable names{p_end}
{synopt:{cmd:e(frame)}}name of the results frame{p_end}
{synopt:{cmd:e(timevar)}}time variable{p_end}
{synopt:{cmd:e(datevar)}}date variable, if specified{p_end}
{synopt:{cmd:e(first)}, {cmd:e(last)}}first and last period of the sample (dates if {cmd:datevar()}){p_end}
{synopt:{cmd:e(asymmetric)}}{cmd:asymmetric} or empty{p_end}
{synopt:{cmd:e(frame_pos)}, {cmd:e(frame_neg)}}results frames of positive and negative changes{p_end}

{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:e(table)}}connectedness table, (k+4) x (k+1); (k+3) x (k+1) for {cmd:joint}{p_end}
{synopt:{cmd:e(fevd)}}averaged decomposition in percent, k x k{p_end}
{synopt:{cmd:e(table_pos)}, {cmd:e(fevd_pos)}}same, positive changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(table_neg)}, {cmd:e(fevd_neg)}}same, negative changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(table_f}{it:b}{cmd:)}, {cmd:e(fevd_f}{it:b}{cmd:)}}same, frequency band {it:b} ({cmd:frequency}){p_end}

{p2col 5 18 22 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{pstd}
The results frame holds one observation per period with measures, the time variable,
and the variables {cmd:tci}, {cmd:ctci}, {cmd:to_}{it:i}, {cmd:from_}{it:i},
{cmd:net_}{it:i}, {cmd:npt_}{it:i}, {cmd:npdc_}{it:i}{cmd:_}{it:j},
{cmd:pci_}{it:i}{cmd:_}{it:j} and {cmd:inf_}{it:i}{cmd:_}{it:j}, where {it:i} and
{it:j} index the variables in {it:varlist}. Variable labels give the names.
{cmd:npdc_}{it:i}{cmd:_}{it:j} is the net transmission from variable {it:i} to
variable {it:j}. Unless {opt nofevd} is specified, the frame also holds
{cmd:fevd_}{it:i}{cmd:_}{it:j}, the share (percent) of variable {it:i}'s forecast
error variance due to variable {it:j}. The frame's characteristics ({cmd:char list})
record the variable names and settings that {helpb cngraph}, {helpb cnnetwork} and
{cmd:estat group} read.

{pstd}
With {cmd:approach(frequency)}, these variables hold the total over the bands, and
each band {it:b} adds the same set with the suffix {cmd:_f}{it:b}, for example
{cmd:tci_f1}, {cmd:net_2_f1} and {cmd:fevd_1_2_f1}. With {cmd:approach(joint)},
the frame holds {cmd:tci}, {cmd:to_}{it:i}, {cmd:from_}{it:i} and
{cmd:net_}{it:i}, the joint measures, and {cmd:fevd_}{it:i}{cmd:_}{it:j}, the
generalized decomposition. With {cmd:approach(extjoint)}, all variables refer to
the extended joint table.


{marker methods}{...}
{title:Methods and formulas}

{pstd}
Let {bf:A}_h be the moving-average matrices of the VAR and {bf:Sigma} the error
covariance. The generalized decomposition is

{p 8 8 2}
theta_ij = sigma_jj^(-1) sum_{h=0}^{H-1} (e_i' {bf:A}_h {bf:Sigma} e_j)^2,
normalized so that each row sums to one, and multiplied by 100.

{pstd}
With w_ij the normalized share of variable i's forecast error variance due to
shocks in j: TO_j = sum_{i!=j} w_ij, FROM_i = sum_{j!=i} w_ij,
NET_i = TO_i - FROM_i, NPDC_ij = w_ji - w_ij, TCI = (1/k) sum_i FROM_i,
corrected TCI = TCI k/(k-1) and
PCI_ij = 200 (w_ij + w_ji)/(w_ii + w_ij + w_ji + w_jj).

{pstd}
{it:Frequency.} Let {bf:Psi}(w) = sum_{h=0}^{H} {bf:A}_h exp(-i w h) at the
Fourier frequencies w_r = 2 pi r/(H+1), r = 0, ..., H. With generalized shocks,
the contribution of frequency w is |({bf:Psi}(w) {bf:Sigma})_ij|^2 /
(sigma_jj d_i), where d_i is the frequency-domain forecast error variance of
variable i over the frequencies used; each row is then scaled so that it sums to
one over those frequencies. A band's table sums these contributions over the
frequencies in the band. The band PCI and influence index use the total table in
their denominators. These steps follow {cmd:FrequencyConnectedness()} in the R
package.

{pstd}
{it:Joint.} FROM_i is the share of variable i's forecast error variance explained
by the joint shock of all other variables,
100 sum_h e_i' {bf:A}_h {bf:Sigma} {bf:M}_i ({bf:M}_i' {bf:Sigma} {bf:M}_i)^(-1)
{bf:M}_i' {bf:Sigma} {bf:A}_h' e_i / (sum_h e_i' {bf:A}_h {bf:Sigma} {bf:A}_h' e_i),
where {bf:M}_i removes variable i; the joint TCI is the mean of FROM_i, and TO_j
rescales the generalized TO_j by the ratio of the generalized and joint TCIs. The
extended joint table multiplies row i of the generalized table by FROM_i divided by
its off-diagonal row sum, so that all measures of the time domain are available.

{pstd}
For the TVP-VAR, the coefficient prediction covariance is the filtered covariance
divided by kappa1, and the error covariance is updated as
kappa2 {bf:Sigma}_{t-1} + (1-kappa2) u_t u_t'.


{marker references}{...}
{title:References}

{phang}
Adekoya, O. B., A. B. Akinseye, N. Antonakakis, I. Chatziantoniou, D. Gabauer, and
J. Oliyide. 2022. Crude oil and Islamic sectoral stocks: Asymmetric TVP-VAR
connectedness and investment strategies. {it:Resources Policy} 78: 102877.

{phang}
Antonakakis, N., I. Chatziantoniou, and D. Gabauer. 2020. Refined measures of dynamic
connectedness based on time-varying parameter vector autoregressions.
{it:Journal of Risk and Financial Management} 13: 84.

{phang}
Balcilar, M., D. Gabauer, and Z. Umar. 2021. Crude oil futures contracts and
commodity markets: New evidence from a TVP-VAR extended joint connectedness approach.
{it:Resources Policy} 73: 102219.

{phang}
Baruník, J., and T. Křehlík. 2018. Measuring the frequency dynamics of
financial connectedness and systemic risk. {it:Journal of Financial Econometrics}
16: 271-296.

{phang}
Chatziantoniou, I., and D. Gabauer. 2021. EMU risk-synchronisation and financial
fragility through the prism of dynamic connectedness.
{it:Quarterly Review of Economics and Finance} 79: 1-14.

{phang}
Chatziantoniou, I., D. Gabauer, and R. Gupta. 2023. Integration and risk
transmission in the market for crude oil: New evidence from a time-varying parameter
frequency connectedness approach. {it:Resources Policy} 84: 103729.

{phang}
Diebold, F. X., and K. Yilmaz. 2009. Measuring financial asset return and volatility
spillovers, with application to global equity markets. {it:Economic Journal} 119: 158-171.

{phang}
Diebold, F. X., and K. Yilmaz. 2012. Better to give than to receive: Predictive
directional measurement of volatility spillovers.
{it:International Journal of Forecasting} 28: 57-66.

{phang}
Gabauer, D. 2025. ConnectednessApproach: Connectedness approach. R package version 1.0.4.
{browse "https://CRAN.R-project.org/package=ConnectednessApproach"}.

{phang}
Koop, G., and D. Korobilis. 2013. Large time-varying parameter VARs.
{it:Journal of Econometrics} 177: 185-198.

{phang}
Koop, G., M. H. Pesaran, and S. M. Potter. 1996. Impulse response analysis in
nonlinear multivariate models. {it:Journal of Econometrics} 74: 119-147.

{phang}
Lastrapes, W. D., and T. F. P. Wiesen. 2021. The joint spillover index.
{it:Economic Modelling} 94: 681-691.

{phang}
Pesaran, M. H., and Y. Shin. 1998. Generalized impulse response analysis in linear
multivariate models. {it:Economics Letters} 58: 17-29.


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
Help:  {helpb connectedness postestimation}, {helpb cngraph}, {helpb cnnetwork}, {manhelp var TS}{p_end}
