{smcl}
{* *! version 0.9.2  24oct2026}{...}
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
{bf:connectedness} {hline 2} Connectedness measures from VAR-type, TVP-VAR and DCC-GARCH models


{marker syntax}{...}
{title:Syntax}

{p 8 22 2}
{cmd:connectedness} {varlist} {ifin} [{cmd:,} {it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt m:odel(model)}}{cmd:var} (default), {cmd:qvar}, {cmd:lad}, {cmd:lasso}, {cmd:ridge}, {cmd:elasticnet}, {cmd:tvpvar} or {cmd:dccgarch}{p_end}
{synopt:{opt l:ags(#)}}VAR lag order; default is {cmd:lags(1)}{p_end}
{synopt:{opt w:indow(#)}}rolling-window size ({cmd:var}, {cmd:qvar}, {cmd:lad}) or prior training sample ({cmd:tvpvar}); default 0{p_end}
{synopt:{opt quant:ile(numlist)}}quantile tau for {cmd:qvar}, one value or one per variable; default {cmd:0.5}{p_end}

{syntab:LASSO, ridge and elastic net}
{synopt:{opt alp:ha(#)}}elastic-net mixing parameter for {cmd:elasticnet}; default is to choose from a grid{p_end}
{synopt:{opt nalp:ha(#)}}grid of {it:#}+1 values from 0 to 1 for {cmd:elasticnet}; default {cmd:10}{p_end}
{synopt:{opt nfold:s(#)}}cross-validation folds; default {cmd:10}{p_end}
{synopt:{opt loss(loss)}}cross-validation loss, {cmd:mae} (default) or {cmd:mse}{p_end}
{synopt:{opt fold:s(folds)}}{cmd:random} (default) or {cmd:interleaved} fold assignment{p_end}
{synopt:{opt seed(#)}}random-number seed for the folds{p_end}
{synopt:{opt adapt:ive}}adaptive LASSO, ridge or elastic net, with penalty weights 1/|b_OLS|{p_end}

{syntab:TVP-VAR}
{synopt:{opt kappa1(#)}}forgetting factor for the coefficients; default {cmd:0.99}{p_end}
{synopt:{opt kappa2(#)}}decay factor for the error covariance; default {cmd:0.99}{p_end}
{synopt:{opt pr:ior(prior)}}{cmd:bayes} (default), {cmd:minnesota} or {cmd:uninformative}{p_end}
{synopt:{opt gam:ma(#)}}Minnesota prior tightness; default {cmd:0.01}{p_end}

{syntab:DCC-GARCH}
{synopt:{opt stand:ardize}}subtract the response without the shock (standardized volatility impulse responses){p_end}
{synopt:{opt dccp:arams(matname)}}use the given 4k+2 parameters instead of estimating them{p_end}

{syntab:Connectedness}
{synopt:{opt app:roach(approach)}}{cmd:time} (default), {cmd:frequency}, {cmd:joint}, {cmd:extjoint}, {cmd:r2} or {cmd:modelfree}{p_end}
{synopt:{opt h:orizon(#)}}forecast horizon {it:H}; default {cmd:horizon(10)}{p_end}
{synopt:{opt gen:eralized}}generalized FEVD; the default{p_end}
{synopt:{opt orth:ogonalized}}Cholesky (orthogonalized) FEVD{p_end}
{synopt:{opt clus:ters(clusters)}}clustered identification: shocks correlated within clusters, orthogonal across them{p_end}
{synopt:{opt clord:er(order)}}{cmd:average} (default) over all orderings of the clusters, or {cmd:given}{p_end}
{synopt:{opt norm:alize(norm)}}normalization of the generalized FEVD: {cmd:row} (default), {cmd:maxrow}, {cmd:maxcol} or {cmd:spectral}{p_end}
{synopt:{opt corr:ected}}report the corrected total connectedness index{p_end}
{synopt:{opt asym:metric}}also estimate the systems of positive and negative changes{p_end}
{synopt:{opt pos:itive(varlist)}}positive series, such as positive realized semivariances, for asymmetric connectedness; requires {opt negative()}{p_end}
{synopt:{opt neg:ative(varlist)}}negative series, such as negative realized semivariances{p_end}

{syntab:R2 and model-free}
{synopt:{opt corm:ethod(method)}}{cmd:pearson} (default), {cmd:spearman} or {cmd:kendall} correlations{p_end}

{syntab:Frequency}
{synopt:{opt band:s(numlist)}}band limits in periods, for example {cmd:bands(5 20)}{p_end}
{synopt:{opt part:ition(numlist)}}band limits in radians, in decreasing order{p_end}
{synopt:{opt scen:ario(scenario)}}{cmd:abs} (default) or {cmd:wth}{p_end}

{syntab:Results}
{synopt:{opt fr:ame(name)}}frame for the time-varying measures; default {cmd:cn_results}{p_end}
{synopt:{opt replace}}replace an existing results frame{p_end}
{synopt:{opt nofevd}}do not store the full FEVD for every period{p_end}
{synopt:{opt nocov:ariance}}do not store the covariance matrix of every period{p_end}
{synopt:{opt date:var(varname)}}copy a date variable into the results frame for graphs{p_end}

{syntab:Reporting}
{synopt:{opt notab:le}}suppress the connectedness table but show the header and TCI{p_end}
{synopt:{opt dig:its(#)}}decimals in the table; default {cmd:digits(2)}{p_end}
{synoptline}
{p 4 6 2}
You must {helpb tsset} your data before using {cmd:connectedness}; the estimation
sample must be contiguous in time. {it:varlist} may contain time-series operators.{p_end}
{p 4 6 2}
See {helpb connectedness postestimation} for {cmd:estat table}, {cmd:estat group}
(aggregated, conditional, inclusive, exclusive, internal and external connectedness),
{cmd:estat correlation} and {cmd:estat bootstrap} (bootstrap bands and event tests),
{helpb cngraph} and {helpb cnnetwork} for graphs, and {helpb cnportfolio} for
portfolio weights, hedge ratios and hedging effectiveness.{p_end}


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
The same measures are available from quantile, least absolute deviation, LASSO, ridge
and elastic-net VARs and, for volatility connectedness, from a DCC-GARCH model
(Gabauer 2020); see {opt model()}. Besides the time domain, {opt approach()} gives
frequency, joint, extended joint, R2 decomposed and model-free connectedness.

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
of the sample. {cmd:qvar} fits the quantile VAR of Chatziantoniou, Gabauer and
Stenfors (2021): each equation is a linear quantile regression at quantile tau,
estimated like R's {cmd:quantreg::rq()} with its default Frisch-Newton
interior-point method. {cmd:lad} fits a least absolute deviation (median) VAR, the
exact median-regression solution as in R's {cmd:L1pack::lad()}. For {cmd:qvar} and
{cmd:lad} the error covariance is the uncentered residual cross-product divided by
the number of observations, as in the R package.

{pstd}
{cmd:lasso}, {cmd:ridge} and {cmd:elasticnet} fit penalized VARs as in Gabauer,
Gupta, Marfatia and Miller (2024) and R's {cmd:ElasticNetVAR()}: each equation is an
elastic net (Zou and Hastie 2005) of the variable on the lags of all variables,
with the penalty chosen by cross-validation as the largest lambda within one
standard error of the minimum (lambda.1se), using the algorithm of R's
{cmd:glmnet} (Friedman, Hastie and Tibshirani 2010): coordinate descent with
standardized regressors on glmnet's lambda path. {cmd:lasso} sets alpha = 1 and
{cmd:ridge} alpha = 0; {cmd:elasticnet} fits the grid of alphas and keeps, for each
equation, the one with the smallest in-sample mean absolute error. The error
covariance is the uncentered residual cross-product divided by the number of
observations. These models can be slow with rolling windows: each window
requires (folds + 1) elastic-net paths per alpha and equation.

{phang}
{opt alpha(#)}, {opt nalpha(#)}, {opt nfolds(#)} and {opt loss()} correspond to
the {cmd:ElasticNet} settings of the R package. {opt folds(random)} assigns the
observations at random to folds of nearly equal size, as {cmd:cv.glmnet()} does,
drawing new folds for each cross-validation, so results depend on the
random-number seed; use {opt seed(#)} for reproducible results. Because R and Stata
have different random-number generators, the two programs give the same results
only with the same folds; {opt folds(interleaved)} assigns observation {it:i} to
fold 1 + mod({it:i}-1, {it:nfolds}) and is used by the certification script.

{phang}
{opt adaptive} fits the adaptive LASSO, ridge or elastic net used by Demirer,
Diebold, Liu and Yilmaz (2018): each slope of an equation is penalized with weight
1/|b_j|, where b_j is its OLS estimate in the regression of the equation on all
lags with an intercept, so that coefficients with small OLS estimates are shrunk
most. The weights are passed to the elastic net as {cmd:glmnet}'s
{cmd:penalty.factor}: they are rescaled to sum to the number of regressors and apply
to the standardized coefficients, and the cross-validation is otherwise unchanged.
With rolling windows the weights are recomputed in each window. The R package has
no adaptive option; the results are checked against {cmd:glmnet} itself.

{phang}
{cmd:dccgarch} fits the DCC-GARCH(1,1) model of Engle (2002) and gives the volatility
connectedness of Gabauer (2020), as R's {cmd:ConnectednessApproach(}..., {cmd:model="DCC-GARCH")}.
Each variable follows a GARCH(1,1) (Bollerslev 1986) with a constant mean and normal
errors; the conditional correlations follow a DCC(1,1). The model is estimated in two
steps by maximum likelihood, as R's {cmd:rmgarch::dccfit()} does: first each GARCH
margin, then the DCC parameters from the standardized residuals. The decomposition
comes from the volatility impulse responses of Gabauer (2020), computed for every
period of the sample, so {opt window()} is not allowed and {opt lags()} is not used.
{opt orthogonalized} and {opt asymmetric} are not allowed. At least 100 observations
are required. Besides {cmd:approach(time)}, {cmd:approach(r2)} gives the DCC-GARCH R2
decomposed connectedness of Cocca, Gabauer and Pomberger (2024): the R2 decomposition
is applied to the conditional correlation matrix of each period, as R's
{cmd:R2Correlations()} applied to {cmd:rmgarch::rcor()}, instead of to the
volatility impulse responses (see {it:Methods and formulas}).

{phang}
{opt quantile(numlist)} sets the quantile tau, strictly between 0 and 1, for
{cmd:model(qvar)}. With one value per variable, equation {it:i} is estimated at the
{it:i}th quantile, which gives the quantile-on-quantile connectedness of Gabauer and
Stenfors (2024) when two variables are used.

{phang}
{opt lags(#)} sets the VAR lag order {it:p}.

{phang}
{opt window(#)} has two meanings, following the R package. With {cmd:model(var)},
each rolling regression uses {it:#} - {it:p} observations, so a window of 200
with one lag uses 199 rows of data. With {cmd:model(tvpvar)} and {cmd:prior(bayes)},
the prior is an OLS VAR on the first {it:#} - {it:p} observations;
{cmd:window(0)} uses the whole sample.

{dlgtab:DCC-GARCH}

{phang}
{opt standardize} computes the volatility impulse responses as the difference between
the conditional covariance paths with and without the shock, as the
{cmd:DCC_config=list(standardize=TRUE)} setting of the R package does. By default the
path with the shock is used.

{phang}
{opt dccparams(matname)} skips the estimation and uses the parameters in the Stata
matrix {it:matname}, a row or column vector of 4k+2 values in the order of
{cmd:rmgarch}: mu, omega, alpha1 and beta1 of each variable in {it:varlist} order,
followed by the DCC parameters a and b. The parameters estimated by R can thus be
entered to reproduce its results exactly; the certification script does this.

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
{opt clusters(clusters)} gives the clustered network connectedness of Buchwalter,
Diebold and Yilmaz (2025). The variables are grouped into clusters, such as regions
or countries, listed as in {cmd:estat group}: each cluster in double quotes,
optionally preceded by a name and a colon, for example
{cmd:clusters("Americas: US BRA" "Europe: UK FRA GER")}. Variables not listed form
clusters of their own. Shocks may be correlated within a cluster, as in the
generalized decomposition, but are orthogonal across clusters, as in the Cholesky
decomposition: the residuals of each cluster are projected on those of the clusters
before it. One cluster thus gives the generalized FEVD, and one variable per cluster
the Cholesky FEVD. Besides the usual measures, {cmd:connectedness} reports, for each
cluster, the shares of its variables' forecast error variance due to their own shocks
(own), to shocks in the same cluster (co-movement) and to shocks in other clusters
(contagion), which sum to 100, and a cluster-to-cluster contagion table. Clustered
identification applies to {cmd:approach(time)} with VAR-type and TVP-VAR models, and
may not be combined with {opt orthogonalized} or {opt normalize()}.

{phang}
{opt clorder(order)} sets how the ordering of the clusters is handled. Only the
ordering across clusters matters, not the order within them. {cmd:clorder(average)},
the default, averages the row-normalized decompositions over all orderings of the
clusters, as the paper does; this allows up to 6 clusters (720 orderings).
{cmd:clorder(given)} uses the clusters in the order listed, followed by the clusters
of the variables not listed, in {it:varlist} order.

{phang}
{opt normalize(norm)} sets how the generalized FEVD is normalized.
{cmd:normalize(row)}, the default, divides each row by its sum (Diebold and Yilmaz
2012), so that every row sums to 100. Caloia, Cipollini and Muzzioli (2019) show
that row normalization can distort the ranking of net transmitters and receivers, and
propose dividing the whole unnormalized decomposition of Pesaran and Shin (1998) by a
single number: {cmd:maxrow}, its largest row sum (their preferred choice),
{cmd:maxcol}, its largest column sum, or {cmd:spectral}, its spectral radius. The
shares then sum to at most 100 by row, TO and FROM keep their relative sizes, and
TCI, the sum of the off-diagonal shares divided by k, reproduces the values in their
Table 11. {opt normalize()} applies to {cmd:approach(time)} with VAR-type and TVP-VAR
models; {cmd:estat group} requires the row-normalized decomposition. For a linear
VAR, the generalized FEVD of Lanne and Nyberg (2016) coincides with row
normalization.

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
{cmd:r2} gives the R2 decomposed connectedness of Naeem et al. (2024), extended to
contemporaneous and lagged components by Balli et al. (2023), and {cmd:modelfree}
the model-free connectedness of Gabauer, Chatziantoniou and Stenfors (2023). Both
work on the data directly, so {opt model()}, {opt orthogonalized}, {opt asymmetric}
and {opt horizon()} do not apply. With {cmd:r2}, {opt lags(#)} may be zero (only
contemporaneous relations) and {opt window(#)} rows of data give
{it:#} - 2{it:p} rows of (x_t, ..., x_{t-p}) per window, as in the R package;
{cmd:modelfree} uses no lags. The exception is {cmd:model(dccgarch) approach(r2)},
which decomposes the DCC-GARCH conditional correlations of every period (Cocca,
Gabauer and Pomberger 2024); it has no lags and takes neither {opt cormethod()} nor
{opt standardize}.

{phang}
{opt cormethod(method)} sets the correlation used by {cmd:approach(r2)} and
{cmd:approach(modelfree)}: Pearson, Spearman (rank) or Kendall (tau-b). For the
R2 decomposition, rank correlations are converted as in the R package,
2 sin(pi r/6) for Spearman and sin(pi tau/2) for Kendall.

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
three sets of measures automatically. The spillover asymmetry measure, the TCI of the
positive system minus the TCI of the negative one, is stored in {cmd:e(SAM)} and,
for every period, in the results frame (see {it:Stored results}).

{phang}
{opt positive(varlist)} and {opt negative(varlist)} give the positive and negative
series directly, one for each variable in {it:varlist}, in the same order. With
realized volatilities in {it:varlist} and the positive and negative realized
semivariances in {opt positive()} and {opt negative()}, this gives the good and bad
volatility spillovers and the spillover asymmetry measure (SAM) of Baruník, Kočenda
and Vácha (2016, 2017). The three systems are estimated and stored as with
{opt asymmetric}, which may not be combined with these options.

{phang}
{opt nofevd} stops {cmd:connectedness} from storing the full decomposition,
{cmd:fevd_}{it:i}{cmd:_}{it:j}, for every period in the results frame. The
decomposition is needed by {cmd:estat group}; omit it only when k is large and the
number of variables would exceed {help maxvar}. ({cmd:savefevd}, from version 0.1,
is still accepted and has no effect.)

{phang}
{opt nocovariance} stops {cmd:connectedness} from storing the covariance matrix of
every period, {cmd:cov_}{it:i}{cmd:_}{it:j} for {it:i} <= {it:j}, in the results
frame: the residual covariance of the VAR, of each rolling window or of each TVP-VAR
period, or the conditional covariance of the DCC-GARCH model. {helpb cnportfolio}
and {cmd:estat correlation} use these matrices. They are not stored with
{cmd:approach(r2)} or {cmd:approach(modelfree)} on the data, which estimate no
covariance matrix (they are stored with {cmd:model(dccgarch) approach(r2)}), and they
are left out with a note when they would exceed {help maxvar}.

{phang}
{opt datevar(varname)} copies a date variable into the results frame; {helpb cngraph}
then uses it as the horizontal axis, and the header and {helpb cnnetwork} report dates.
This is useful when the data are {cmd:tsset} by an observation index, as with daily
trading data. Without {opt datevar()}, a numeric variable named {cmd:date} with a date
format is used automatically when the time variable has no date format, as in the
datasets installed with the package; {cmd:datevar(}{it:timevar}{cmd:)}, naming the
time variable itself, turns this off.


{marker remarks}{...}
{title:Remarks}

{pstd}
Typing {cmd:connectedness} without arguments redisplays the results. Options
{opt digits(#)}, {opt notable}, {opt part(positive|negative)} (after
{opt asymmetric}) and {opt band(#)} (after {cmd:approach(frequency)}, to show the
table of band {it:#}) and {opt component(contemporaneous|lagged)} (after
{cmd:approach(r2)} with lags) are allowed.

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
With {cmd:model(dccgarch)}, the conditional covariances, the log likelihood and the
volatility decomposition follow {cmd:rmgarch} and {cmd:VFEVD()} exactly, including
their conventions: the GARCH recursion starts at the mean squared residual; the DCC
parameters are estimated with a presample standardized residual of zero, while the
reported correlations and the log likelihood are filtered, as in {cmd:rmgarch}, with
a presample value of one. The estimates themselves come from a different optimizer
(Newton-Raphson with analytic gradients here, {cmd:solnp} in R), so they agree with
R's to the precision of R's optimizer, about 1e-5 in the certification cases; with
the same parameters ({opt dccparams()}) the results agree to numerical precision.
Stata's estimates attain a likelihood at least as high as R's.

{pstd}
Data from 19 published applications are installed with the package as
{cmd:cn_}{it:name}{cmd:.dta}, for example {cmd:cn_dy2012.dta} and
{cmd:cn_acg2020.dta}. They are ancillary files: type {cmd:net get connectedness} to copy
them to the current directory, then load them with {cmd:use}. Each is {cmd:tsset}-ready through the observation index
{cmd:t} (monthly data also carry {cmd:month}, formatted {cmd:%tm}) and has a {cmd:date} variable.


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

{pstd}Chatziantoniou, Gabauer and Stenfors (2021): quantile connectedness{p_end}
{phang2}{cmd:. use cn_cgs2021, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness USD1Y EUR1Y GBP1Y JPY1Y, model(qvar) quantile(0.05) horizon(20) window(200) corrected datevar(date) replace}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}

{pstd}Gabauer, Gupta, Marfatia and Miller (2024): LASSO, ridge and elastic-net VAR connectedness{p_end}
{phang2}{cmd:. use cn_gghm2022, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness AK-WY, model(lasso) lags(1) seed(2024) nofevd notable replace}{p_end}
{phang2}{cmd:. cnnetwork, threshold(0.5)}{p_end}

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
{phang2}{cmd:. connectedness energy-utilities, lags(2) horizon(30) approach(joint) replace}{p_end}

{pstd}Balcilar, Gabauer and Umar (2021): extended joint connectedness{p_end}
{phang2}{cmd:. use cn_bgu2021, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness CrudeOil-Wheat, model(tvpvar) horizon(20) prior(minnesota) gamma(0.1) approach(extjoint) datevar(date) replace}{p_end}
{phang2}{cmd:. cngraph net}{p_end}

{pstd}Balli et al. (2023): contemporaneous and lagged R2 decomposed connectedness{p_end}
{phang2}{cmd:. use cn_dy2012, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, approach(r2) lags(1) window(200) datevar(date) replace}{p_end}
{phang2}{cmd:. connectedness, component(lagged)}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}

{pstd}Gabauer, Chatziantoniou and Stenfors (2023): model-free connectedness{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, approach(modelfree) window(200) cormethod(kendall) datevar(date) replace}{p_end}

{pstd}Gabauer (2020): DCC-GARCH volatility connectedness{p_end}
{phang2}{cmd:. use cn_g2020, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness EUR GBP CHF JPY, model(dccgarch) datevar(date) replace}{p_end}
{phang2}{cmd:. matrix list e(garch)}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}
{phang2}{cmd:. estat correlation}{p_end}

{pstd}Cocca, Gabauer and Pomberger (2024): DCC-GARCH R2 decomposed connectedness{p_end}
{phang2}{cmd:. use cn_cgp2024, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness GRE SOL WIN GEO, model(dccgarch) approach(r2) datevar(date) replace}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}
{phang2}{cmd:. cnnetwork}{p_end}

{pstd}Chatziantoniou, Abakah, Gabauer and Tiwari (2022): quantile frequency connectedness{p_end}
{phang2}{cmd:. use cn_gcat2022, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness S_P_GREEN_BOND_INDEX-S_P_GLOBAL_CLEAN_ENERGY, model(qvar) quantile(0.05) horizon(100) window(200) approach(frequency) bands(5) datevar(date) replace}{p_end}
{phang2}{cmd:. connectedness, band(1)}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}

{pstd}Adekoya et al. (2022): asymmetric connectedness of positive and negative changes{p_end}
{phang2}{cmd:. use cn_aaacgo2022, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. foreach v of varlist Brent-UTI {c -(}}{p_end}
{phang2}{cmd:.     generate double r_`v' = 100*(`v' - L.`v')/L.`v'}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. connectedness r_*, model(tvpvar) prior(minnesota) gamma(0.1) asymmetric datevar(date) replace}{p_end}
{phang2}{cmd:. connectedness, part(negative)}{p_end}
{phang2}{cmd:. cngraph tci}{p_end}

{pstd}Baruník, Kočenda and Vácha (2016): good and bad volatility spillovers and SAM
(squared positive and negative daily returns stand in for realized semivariances){p_end}
{phang2}{cmd:. use cn_g2020, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. foreach v in EUR GBP CHF JPY {c -(}}{p_end}
{phang2}{cmd:.     generate double rv_`v' = `v'^2}{p_end}
{phang2}{cmd:.     generate double rp_`v' = cond(`v' > 0, `v'^2, 0)}{p_end}
{phang2}{cmd:.     generate double rn_`v' = cond(`v' < 0, `v'^2, 0)}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. connectedness rv_*, positive(rp_*) negative(rn_*) window(500) datevar(date) replace}{p_end}
{phang2}{cmd:. frame cn_results: tsline sam}{p_end}

{pstd}Buchwalter, Diebold and Yilmaz (2025): clustered connectedness of equity markets by region{p_end}
{phang2}{cmd:. use cn_dy2009, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness US-TUR, lags(2) clusters("Americas: US ARG BRA CHL MEX" "Europe: UK FRA GER TUR" "AsiaPacific: HKG JPN AUS IDN KOR MYS PHL SGP TAI THA") replace}{p_end}
{phang2}{cmd:. matrix list e(cl_contagion)}{p_end}

{pstd}Caloia, Cipollini and Muzzioli (2019): scalar normalization of the Diebold and Yilmaz (2012) table{p_end}
{phang2}{cmd:. use cn_dy2012, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10)}{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) normalize(maxrow)}{p_end}

{pstd}Demirer, Diebold, Liu and Yilmaz (2018): adaptive LASSO VAR{p_end}
{phang2}{cmd:. connectedness SP500 R_10Y DJUBSCOM USDX, model(lasso) adaptive lags(2) seed(2018)}{p_end}


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
{synopt:{cmd:e(quantile)}}quantile tau, when one value is used ({cmd:qvar}){p_end}
{synopt:{cmd:e(nfolds)}}cross-validation folds ({cmd:lasso}, {cmd:ridge}, {cmd:elasticnet}){p_end}
{synopt:{cmd:e(TCI_pos)}, {cmd:e(cTCI_pos)}}TCI of positive changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(TCI_neg)}, {cmd:e(cTCI_neg)}}TCI of negative changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(nclusters)}}number of clusters ({opt clusters()}){p_end}
{synopt:{cmd:e(norderings)}}number of cluster orderings averaged ({opt clusters()}){p_end}
{synopt:{cmd:e(SAM)}}spillover asymmetry measure, {cmd:e(TCI_pos)} - {cmd:e(TCI_neg)} ({cmd:asymmetric}, {opt positive()}){p_end}
{synopt:{cmd:e(nbands)}}number of frequency bands ({cmd:frequency}){p_end}
{synopt:{cmd:e(TCI_c)}, {cmd:e(TCI_l)}}TCI of the contemporaneous and lagged components ({cmd:r2} with lags); also {cmd:e(cTCI_c)}, {cmd:e(cTCI_l)}{p_end}
{synopt:{cmd:e(TCI_f}{it:b}{cmd:)}, {cmd:e(cTCI_f}{it:b}{cmd:)}}TCI and corrected TCI of band {it:b} ({cmd:frequency}){p_end}
{synopt:{cmd:e(ll)}}log likelihood of the DCC-GARCH model ({cmd:dccgarch}){p_end}
{synopt:{cmd:e(converged)}}1 if the optimizations converged, 0 otherwise ({cmd:dccgarch}){p_end}
{synopt:{cmd:e(covariance)}}1 if the covariances are stored in the results frame, 0 otherwise{p_end}

{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:connectedness}{p_end}
{synopt:{cmd:e(model)}}{cmd:var}, {cmd:qvar}, {cmd:lad}, {cmd:lasso}, {cmd:ridge}, {cmd:elasticnet}, {cmd:tvpvar} or {cmd:dccgarch}{p_end}
{synopt:{cmd:e(quantiles)}}quantile of each equation ({cmd:qvar}){p_end}
{synopt:{cmd:e(alphas)}}alpha values considered ({cmd:lasso}, {cmd:ridge}, {cmd:elasticnet}){p_end}
{synopt:{cmd:e(loss)}, {cmd:e(folds)}, {cmd:e(seed)}}cross-validation settings{p_end}
{synopt:{cmd:e(adaptive)}}{cmd:adaptive} or empty ({cmd:lasso}, {cmd:ridge}, {cmd:elasticnet}){p_end}
{synopt:{cmd:e(approach)}}{cmd:time}, {cmd:frequency}, {cmd:joint}, {cmd:extjoint}, {cmd:r2} or {cmd:modelfree}{p_end}
{synopt:{cmd:e(cormethod)}}correlation used ({cmd:r2}, {cmd:modelfree}): {cmd:pearson}, {cmd:spearman}, {cmd:kendall}, or {cmd:dcc} for the DCC-GARCH conditional correlations{p_end}
{synopt:{cmd:e(bands)}}band labels in periods, for example {cmd:1-5 5-Inf} ({cmd:frequency}){p_end}
{synopt:{cmd:e(partition)}}band limits in radians ({cmd:frequency}){p_end}
{synopt:{cmd:e(scenario)}}{cmd:abs} or {cmd:wth} ({cmd:frequency}){p_end}
{synopt:{cmd:e(fevdtype)}}type of decomposition, for example {cmd:generalized}{p_end}
{synopt:{cmd:e(normalize)}}{cmd:row}, {cmd:maxrow}, {cmd:maxcol} or {cmd:spectral} (generalized FEVD, {cmd:approach(time)}){p_end}
{synopt:{cmd:e(standardize)}}{cmd:standardize} or empty ({cmd:dccgarch}){p_end}
{synopt:{cmd:e(dccparams)}}name of the matrix in {opt dccparams()} ({cmd:dccgarch}){p_end}
{synopt:{cmd:e(corrected)}}{cmd:corrected} or empty{p_end}
{synopt:{cmd:e(prior)}}prior ({cmd:tvpvar}){p_end}
{synopt:{cmd:e(names)}}variable names{p_end}
{synopt:{cmd:e(frame)}}name of the results frame{p_end}
{synopt:{cmd:e(timevar)}}time variable{p_end}
{synopt:{cmd:e(datevar)}}date variable, if specified{p_end}
{synopt:{cmd:e(first)}, {cmd:e(last)}}first and last period of the sample (dates if {cmd:datevar()}){p_end}
{synopt:{cmd:e(asymmetric)}}{cmd:asymmetric} or empty{p_end}
{synopt:{cmd:e(frame_pos)}, {cmd:e(frame_neg)}}results frames of positive and negative changes{p_end}
{synopt:{cmd:e(clusters)}}cluster names ({opt clusters()}){p_end}
{synopt:{cmd:e(clorder)}}{cmd:average} or {cmd:given} ({opt clusters()}){p_end}
{synopt:{cmd:e(positive)}, {cmd:e(negative)}}the variables in {opt positive()} and {opt negative()}{p_end}

{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:e(table)}}connectedness table, (k+4) x (k+1); (k+3) x (k+1) for {cmd:joint}{p_end}
{synopt:{cmd:e(fevd)}}averaged decomposition in percent, k x k{p_end}
{synopt:{cmd:e(table_pos)}, {cmd:e(fevd_pos)}}same, positive changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(table_neg)}, {cmd:e(fevd_neg)}}same, negative changes ({cmd:asymmetric}){p_end}
{synopt:{cmd:e(table_f}{it:b}{cmd:)}, {cmd:e(fevd_f}{it:b}{cmd:)}}same, frequency band {it:b} ({cmd:frequency}){p_end}
{synopt:{cmd:e(table_c)}, {cmd:e(fevd_c)}}same, contemporaneous component ({cmd:r2} with lags); {cmd:e(table_l)}, {cmd:e(fevd_l)} for the lagged one{p_end}
{synopt:{cmd:e(cl_shares)}}own, co-movement and contagion shares, TO and NET of each cluster, and their averages over clusters, (C+1) x 5 ({opt clusters()}){p_end}
{synopt:{cmd:e(cl_contagion)}}C x C table: share of the row cluster's variance due to the column cluster, averaged over the row cluster's variables{p_end}
{synopt:{cmd:e(cl_contagion_avg)}}the same divided by the number of variables in the column cluster, the definition printed in the paper{p_end}
{synopt:{cmd:e(cl_id)}}cluster of each variable, 1 x k{p_end}
{synopt:{cmd:e(garch)}}GARCH(1,1) parameters mu, omega, alpha and beta, k x 4 ({cmd:dccgarch}){p_end}
{synopt:{cmd:e(dcc)}}DCC(1,1) parameters a and b, 1 x 2 ({cmd:dccgarch}){p_end}

{p2col 5 18 22 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{pstd}
The results frame holds one observation per period with measures, the time variable,
and the variables {cmd:tci}, {cmd:ctci}, {cmd:to_}{it:i}, {cmd:from_}{it:i},
{cmd:net_}{it:i}, {cmd:npt_}{it:i}, {cmd:npdc_}{it:i}{cmd:_}{it:j},
{cmd:pci_}{it:i}{cmd:_}{it:j} and {cmd:inf_}{it:i}{cmd:_}{it:j}, where {it:i} and
{it:j} index the variables in {it:varlist}. Variable labels give the names.
{cmd:npdc_}{it:i}{cmd:_}{it:j} is the net transmission from variable {it:i} to
variable {it:j}. {cmd:inf_}{it:i}{cmd:_}{it:j} is missing in periods where both
pairwise shares are zero. Unless {opt nofevd} is specified, the frame also holds
{cmd:fevd_}{it:i}{cmd:_}{it:j}, the share (percent) of variable {it:i}'s forecast
error variance due to variable {it:j}. Unless {opt nocovariance} is specified, it
holds {cmd:cov_}{it:i}{cmd:_}{it:j} ({it:i} <= {it:j}), the covariance matrix of the
period (see {opt nocovariance}). The frame's characteristics ({cmd:char list})
record the variable names and settings that {helpb cngraph}, {helpb cnnetwork} and
{cmd:estat group} read.

{pstd}
With {cmd:approach(frequency)}, these variables hold the total over the bands, and
each band {it:b} adds the same set with the suffix {cmd:_f}{it:b}, for example
{cmd:tci_f1}, {cmd:net_2_f1} and {cmd:fevd_1_2_f1}. With {cmd:approach(joint)},
the frame holds {cmd:tci}, {cmd:to_}{it:i}, {cmd:from_}{it:i} and
{cmd:net_}{it:i}, the joint measures, and {cmd:fevd_}{it:i}{cmd:_}{it:j}, the
generalized decomposition. With {cmd:approach(extjoint)}, all variables refer to
the extended joint table. With {cmd:approach(r2)} and lags, the variables hold the
overall measures and the suffixes {cmd:_c} and {cmd:_l} the contemporaneous and
lagged components; {cmd:fevd_}{it:i}{cmd:_}{it:j} then holds the share (percent) of
variable {it:i}'s variance explained by variable {it:j}.

{pstd}
With {opt asymmetric} or {opt positive()} and {opt negative()}, the results frame
also holds the spillover asymmetry measures of every period: {cmd:sam}, the TCI of
the positive system minus that of the negative system, and {cmd:sam_to_}{it:i},
{cmd:sam_from_}{it:i} and {cmd:sam_net_}{it:i}, the same differences of TO, FROM and
NET. Graph them with, for example, {cmd:frame cn_results: tsline sam}.

{pstd}
With {opt clusters()}, the results frame also holds, for every period and each cluster
{it:c} (numbered in the order of {cmd:e(clusters)}), {cmd:cl_own_}{it:c},
{cmd:cl_comove_}{it:c} and {cmd:cl_contag_}{it:c}, the own, co-movement and contagion
shares, {cmd:cl_to_}{it:c} and {cmd:cl_net_}{it:c}, the contagion transmitted to the
other clusters and its net value, {cmd:cl_own}, {cmd:cl_comove} and {cmd:cl_contag},
their averages over clusters, and {cmd:cl_}{it:c}{cmd:_}{it:d}, the contagion to
cluster {it:c} from cluster {it:d}.


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
{it:Normalizations.} Dividing theta_ij by d_i = sum_{h=0}^{H-1} e_i' {bf:A}_h
{bf:Sigma} {bf:A}_h' e_i gives the unnormalized decomposition {bf:W} of Pesaran and
Shin (1998), whose rows need not sum to one. {cmd:normalize(row)} uses
w_ij = W_ij / sum_j W_ij; the scalar normalizations of Caloia, Cipollini and Muzzioli
(2019) use w_ij = W_ij / c, with c the largest row sum ({cmd:maxrow}), the largest
column sum ({cmd:maxcol}) or the largest eigenvalue modulus ({cmd:spectral}) of
{bf:W}. The measures above then follow from w_ij unchanged.

{pstd}
{it:Adaptive elastic net.} For equation i with regressors x_j (all lags of all
variables), the coefficients minimize
(1/2n) sum_t (y_t - b_0 - sum_j b_j x_jt)^2 + lambda sum_j v_j (alpha |b_j| + (1-alpha) b_j^2/2)
on standardized regressors, with v_j = 1 for the elastic net and v_j proportional to
1/|b_j^OLS|, rescaled to sum to the number of regressors, with {opt adaptive}, as in
{cmd:glmnet}; lambda is chosen by cross-validation as without {opt adaptive}.

{pstd}
{it:Clustered identification.} For an ordering of the clusters, let the residuals of
cluster b be projected on those of the clusters before it,
{bf:Q}^(-1)_(b,<b) = -{bf:Sigma}_(b,<b) {bf:Sigma}_(<b,<b)^(-1), with identity
blocks on the diagonal, so that {bf:Omega} = {bf:Q}^(-1) {bf:Sigma}
{bf:Q}^(-1)' is block diagonal (Buchwalter, Diebold and Yilmaz 2025, eqs. 7 and 8).
The decomposition is
omega_jj^(-1) sum_{h=0}^{H-1} (e_i' {bf:A}_h {bf:Q} {bf:Omega} e_j)^2 /
sum_{h=0}^{H-1} e_i' {bf:A}_h {bf:Sigma} {bf:A}_h' e_i, normalized by row; with
{cmd:clorder(average)} it is averaged over all orderings. With N_c the number of
variables in cluster c and w_ij the decomposition, own_c = (1/N_c) sum_{i in c}
w_ii, co-movement_c = (1/N_c) sum_{i,j in c, i!=j} w_ij, contagion_c = (1/N_c)
sum_{i in c} sum_{j not in c} w_ij, TO_c = (1/N_c) sum_{i not in c} sum_{j in c}
w_ij and NET_c = TO_c - contagion_c; the contagion table holds (1/N_c) sum_{i in c}
sum_{j in d} w_ij. The paper prints the contagion from cluster d to c with a further
division by N_d; with that definition the own, co-movement and contagion shares do
not sum to one, as the paper states they do, unless every cluster has one variable,
so the package reports the version that sums to 100 and stores the printed one in
{cmd:e(cl_contagion_avg)}.

{pstd}
{it:Spillover asymmetry.} With TCI+ and TCI- the total connectedness of the
positive and negative systems, SAM = TCI+ - TCI- (Baruník, Kočenda and Vácha 2016),
and the directional versions are the differences of TO, FROM and NET.

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
{it:R2 decomposition.} For each window and variable i, the R2 of a regression of
x_it on the other variables at time t and on all variables at lags 1 to p is
split among the regressors by the relative weights of the correlation matrix: with
{bf:R}_xx = {bf:V} {bf:L} {bf:V}' and {bf:S} = {bf:V} {bf:L}^(1/2) {bf:V}', the
shares are ({bf:S}.^2) ({bf:S}^(-1) {bf:r}_xy).^2. Contemporaneous shares give the
contemporaneous table, lag shares summed over lags the lagged table, and their sum
the overall table, as in {cmd:R2Connectedness()}. {it:Model-free.} The table is
the matrix of squared correlations with each row divided by its sum.

{pstd}
{it:DCC-GARCH.} With residuals r_it = x_it - mu_i, the conditional variances are
h_it = omega_i + alpha_i r_i,t-1^2 + beta_i h_i,t-1, starting at the mean of r_it^2.
With z_t the standardized residuals and {bf:Qbar} their sample covariance,
{bf:Q}_t = (1-a-b) {bf:Qbar} + a z_t-1 z_t-1' + b {bf:Q}_t-1, {bf:R}_t is {bf:Q}_t
scaled to unit diagonal, and {bf:H}_t = {bf:D}_t {bf:R}_t {bf:D}_t with
{bf:D}_t = diag(h_t)^(1/2). For each period, the volatility impulse response
follows {cmd:VFEVD()}: a unit shock to all variables enters the GARCH recursions at
the first step, h_1 = alpha + beta.*diag({bf:H}_t), and afterwards
h_j+1 = beta.*h_j; the correlations follow {bf:Q}_1 = (1-a-b) {bf:Qbar} +
a diag(h_1) + b {bf:H}_t and {bf:Q}_j+1 = (1-a-b) {bf:Qbar} + (a+b) {bf:Q}_j. With
{bf:V}_j the covariance path (minus the path without the shock under
{opt standardize}), the table is the transpose of sum_{j=0}^{H} {bf:V}_j.^2 with
each row divided by its sum.

{pstd}
{it:DCC-GARCH R2 decomposition.} With {cmd:approach(r2)}, the conditional correlation
matrix {bf:R}_t of each period is decomposed as in {cmd:R2Correlations()}: row {it:i}
of the table holds the relative-weight shares of the R2 of variable {it:i} regressed
on the others, computed from {bf:R}_t as for the data-based R2 decomposition above,
and the own share is 1, so the table is not row-normalized. A period in which the
shares of variable {it:i} sum to more than 1 (an indefinite {bf:R}_t submatrix)
takes the previous period's shares for that row, as R does. In the first period, R
stops with an error; Stata keeps the shares and reports a note. The pairwise
connectedness index of {cmd:R2Correlations()} is on a scale of 0 to 2; the results
frame holds it in percent, like the other approaches.

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
Balli, F., H. O. Balli, T. H. N. Dang, and D. Gabauer. 2023. Contemporaneous and
lagged R2 decomposed connectedness approach: New evidence from the energy futures
market. {it:Finance Research Letters} 57: 104168.

{phang}
Baruník, J., and T. Křehlík. 2018. Measuring the frequency dynamics of
financial connectedness and systemic risk. {it:Journal of Financial Econometrics}
16: 271-296.

{phang}
Baruník, J., E. Kočenda, and L. Vácha. 2016. Asymmetric connectedness on the U.S.
stock market: Bad and good volatility spillovers. {it:Journal of Financial Markets}
27: 55-78.

{phang}
Baruník, J., E. Kočenda, and L. Vácha. 2017. Asymmetric volatility connectedness on
the forex market. {it:Journal of International Money and Finance} 77: 39-56.

{phang}
Bollerslev, T. 1986. Generalized autoregressive conditional heteroskedasticity.
{it:Journal of Econometrics} 31: 307-327.

{phang}
Buchwalter, B., F. X. Diebold, and K. Yilmaz. 2025. Clustered network connectedness: A
new measurement framework, with application to global equity markets. NBER Working
Paper 34796.

{phang}
Caloia, F. G., A. Cipollini, and S. Muzzioli. 2019. How do normalization schemes
affect net spillovers? A replication of the Diebold and Yilmaz (2012) study.
{it:Energy Economics} 84: 104536.

{phang}
Chatziantoniou, I., E. J. A. Abakah, D. Gabauer, and A. K. Tiwari. 2022. Quantile
time-frequency price connectedness between green bond, green equity, sustainable
investments and clean energy markets. {it:Journal of Cleaner Production} 361: 132088.

{phang}
Chatziantoniou, I., and D. Gabauer. 2021. EMU risk-synchronisation and financial
fragility through the prism of dynamic connectedness.
{it:Quarterly Review of Economics and Finance} 79: 1-14.

{phang}
Chatziantoniou, I., D. Gabauer, and R. Gupta. 2023. Integration and risk
transmission in the market for crude oil: New evidence from a time-varying parameter
frequency connectedness approach. {it:Resources Policy} 84: 103729.

{phang}
Chatziantoniou, I., D. Gabauer, and A. Stenfors. 2021. Interest rate swaps and the
transmission mechanism of monetary policy: A quantile connectedness approach.
{it:Economics Letters} 204: 109891.

{phang}
Cocca, T., D. Gabauer, and S. Pomberger. 2024. Clean energy market connectedness and
investment strategies: New evidence from DCC-GARCH R2 decomposed connectedness
measures. {it:Energy Economics} 136: 107680.

{phang}
Demirer, M., F. X. Diebold, L. Liu, and K. Yilmaz. 2018. Estimating global bank
network connectedness. {it:Journal of Applied Econometrics} 33: 1-15.

{phang}
Diebold, F. X., and K. Yilmaz. 2009. Measuring financial asset return and volatility
spillovers, with application to global equity markets. {it:Economic Journal} 119: 158-171.

{phang}
Diebold, F. X., and K. Yilmaz. 2012. Better to give than to receive: Predictive
directional measurement of volatility spillovers.
{it:International Journal of Forecasting} 28: 57-66.

{phang}
Engle, R. 2002. Dynamic conditional correlation: A simple class of multivariate
generalized autoregressive conditional heteroskedasticity models.
{it:Journal of Business and Economic Statistics} 20: 339-350.

{phang}
Friedman, J., T. Hastie, and R. Tibshirani. 2010. Regularization paths for
generalized linear models via coordinate descent. {it:Journal of Statistical Software}
33(1): 1-22.

{phang}
Gabauer, D. 2020. Volatility impulse response analysis for DCC-GARCH models: The role
of volatility transmission mechanisms. {it:Journal of Forecasting} 39: 788-796.

{phang}
Gabauer, D. 2025. ConnectednessApproach: Connectedness approach. R package version 1.0.4.
{browse "https://CRAN.R-project.org/package=ConnectednessApproach"}.

{phang}
Gabauer, D., I. Chatziantoniou, and A. Stenfors. 2023. Model-free connectedness
measures. {it:Finance Research Letters} 54: 103804.

{phang}
Gabauer, D., R. Gupta, H. A. Marfatia, and S. M. Miller. 2024. Estimating U.S.
housing price network connectedness: Evidence from dynamic elastic net, lasso, and
ridge vector autoregressive models. {it:International Review of Economics and Finance}
89: 349-362.

{phang}
Gabauer, D., and A. Stenfors. 2024. Quantile-on-quantile connectedness measures:
Evidence from the US treasury yield curve. {it:Finance Research Letters} 60: 104852.

{phang}
Ghalanos, A. 2026. rmgarch: Multivariate GARCH models. R package version 1.4-3.
{browse "https://CRAN.R-project.org/package=rmgarch"}.

{phang}
Koenker, R. 2005. {it:Quantile Regression}. Cambridge: Cambridge University Press.

{phang}
Koop, G., and D. Korobilis. 2013. Large time-varying parameter VARs.
{it:Journal of Econometrics} 177: 185-198.

{phang}
Koop, G., M. H. Pesaran, and S. M. Potter. 1996. Impulse response analysis in
nonlinear multivariate models. {it:Journal of Econometrics} 74: 119-147.

{phang}
Lanne, M., and H. Nyberg. 2016. Generalized forecast error variance decomposition for
linear and nonlinear multivariate models. {it:Oxford Bulletin of Economics and
Statistics} 78: 595-603.

{phang}
Lastrapes, W. D., and T. F. P. Wiesen. 2021. The joint spillover index.
{it:Economic Modelling} 94: 681-691.

{phang}
Naeem, M. A., I. Chatziantoniou, D. Gabauer, and S. Karim. 2024. Measuring the G20
stock market return transmission mechanism: Evidence from the R2 connectedness
approach. {it:International Review of Financial Analysis} 91: 102986.

{phang}
Pesaran, M. H., and Y. Shin. 1998. Generalized impulse response analysis in linear
multivariate models. {it:Economics Letters} 58: 17-29.

{phang}
Portnoy, S., and R. Koenker. 1997. The Gaussian hare and the Laplacian tortoise:
Computability of squared-error versus absolute-error estimators.
{it:Statistical Science} 12: 279-300.

{phang}
Zou, H., and T. Hastie. 2005. Regularization and variable selection via the elastic
net. {it:Journal of the Royal Statistical Society, Series B} 67: 301-320.


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
Help:  {helpb connectedness postestimation}, {helpb cngraph}, {helpb cnnetwork}, {helpb cnportfolio}, {manhelp var TS}{p_end}
