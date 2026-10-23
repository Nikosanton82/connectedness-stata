{smcl}
{* *! version 0.8.0  23oct2026}{...}
{title:Title}

{phang}
{bf:cnportfolio} {hline 2} Portfolio weights, hedge ratios and hedging effectiveness after connectedness


{title:Syntax}

{p 8 16 2}
{cmd:cnportfolio} {it:returns} {ifin}{cmd:,} {opt m:ethod(method)} [{it:options}]

{p 4 4 2}
{it:returns} are the return series of the k variables of the last {cmd:connectedness}
estimates, in the same order, as decimal fractions (for example, a percentage change
divided by 100).

{synoptset 24 tabbed}{...}
{synopthdr:method}
{synoptline}
{synopt:{cmd:mcp}}minimum connectedness portfolio (Broadstock, Chatziantoniou and Gabauer 2022){p_end}
{synopt:{cmd:mhp}}multivariate hedging portfolio (Cocca, Gabauer and Pomberger 2024){p_end}
{synopt:{cmd:rpp}}risk-parity portfolio{p_end}
{synopt:{cmd:ewp}}equally weighted portfolio{p_end}
{synopt:{cmd:bivariate}}bivariate portfolio weights of every pair (Kroner and Ng 1998){p_end}
{synopt:{cmd:hedge}}hedge ratios of every pair (Kroner and Sultan 1993){p_end}
{synoptline}

{synopthdr}
{synoptline}
{synopt:{opt mat:rix(matrix)}}{cmd:pci}, {cmd:cov} or {cmd:cor}; default {cmd:pci} for {cmd:mcp} and {cmd:cov} otherwise{p_end}
{synopt:{opt fr:ame(name)}}results frame of {cmd:connectedness}; default the frame of the last estimates{p_end}
{synopt:{opt test(test)}}test of equal variances: {cmd:fisher} (default), {cmd:bartlett}, {cmd:fk}, {cmd:levene} or {cmd:bf}{p_end}
{synopt:{opt cum:ulative(type)}}cumulative returns: {cmd:sum} (default) or {cmd:prod}{p_end}
{synopt:{opt nolong}}allow short positions ({cmd:mcp}, {cmd:bivariate}){p_end}
{synopt:{opt sc:ale(#)}}periods per year for annualization; default from the dates{p_end}
{synopt:{opt sav:ing(name[, replace])}}frame for the series; default {cmd:cn_portfolio}, replaced{p_end}
{synopt:{opt dig:its(#)}}decimals in the table; default 2{p_end}
{synopt:{opt notab:le}}do not display the table{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:cnportfolio} turns connectedness estimates into portfolio and hedging strategies,
following the portfolio functions of the R package {it:ConnectednessApproach}
(Gabauer 2025). In every period with estimates it computes portfolio weights, hedge
weights or hedge ratios from a k x k matrix: the pairwise connectedness index (PCI),
the covariance matrix or the correlation matrix of that period. It then reports, for
each asset or pair, the mean, standard deviation and 5% and 95% quantiles of the
weights or ratios, the hedging effectiveness of the portfolio, a test of equal
variances of the asset and the portfolio, and the portfolio's annualized Sharpe ratio.

{pstd}
The matrices come from the results frame of {helpb connectedness}. Covariances are
stored there by default: the residual covariance of the VAR, of each rolling window
or of each TVP-VAR period, and the conditional covariance of the DCC-GARCH model
({cmd:model(dccgarch)}), the usual choice for hedging. With static estimates
(no {opt window()}), one matrix is used for all observations of {it:returns}; with
time-varying estimates, the periods of the results frame are matched to the data on
the time variable.

{pstd}
The series of weights, ratios and portfolio returns are saved in a frame, by default
{cmd:cn_portfolio}, for graphs or further analysis.


{title:Options}

{phang}
{opt method(method)} is required.

{phang2}
{cmd:mcp} gives the weights w = V^-1 1 / (1'V^-1 1) for the matrix V of each period,
the global minimum-variance weights when V is a covariance matrix and the minimum
connectedness weights when V is the PCI matrix (with 100 on the diagonal), as in
Broadstock, Chatziantoniou and Gabauer (2022). By default the weights are restricted
to [0, 1] and rescaled to sum to 1, as in the R package; {opt nolong} keeps the
unrestricted weights.

{phang2}
{cmd:mhp} hedges each asset i with all the others: the hedged portfolio is
r_i - sum_j b_ij r_j, where b_i are the coefficients of a regression of asset i on the
other assets implied by the covariance matrix of the period (Cocca, Gabauer and
Pomberger 2024). The table reports one row per pair i/j with the summary of b_ij.

{phang2}
{cmd:rpp} gives the long-only risk-parity weights, under which every asset contributes
the same share of the portfolio variance. These weights are unique; they are computed
by Newton's method on the convex formulation of Spinu (2013) and equal those of
{cmd:riskParityPortfolio()} in R.

{phang2}
{cmd:ewp} gives the weights 1/k.

{phang2}
{cmd:bivariate} gives, for each pair, the weight of asset i in the minimum-variance
portfolio of i and j, w = (h_jj - h_ij)/(h_ii - 2h_ij + h_jj) (Kroner and Ng 1998),
restricted to [0, 1] unless {opt nolong} is specified; the row i/j summarizes w and
the portfolio w r_i + (1 - w) r_j.

{phang2}
{cmd:hedge} gives, for each pair, the hedge ratio h_ij/h_jj of a long position in
asset i and a short position in asset j (Kroner and Sultan 1993); the row i/j
summarizes the ratio and the hedged portfolio r_i - (h_ij/h_jj) r_j.

{phang}
{opt matrix(matrix)} chooses the matrix: {cmd:pci} reads the pairwise connectedness
indices {cmd:pci_}{it:i}{cmd:_}{it:j}, {cmd:cov} the covariances
{cmd:cov_}{it:i}{cmd:_}{it:j}, and {cmd:cor} the correlations implied by them. The
covariances are not available after {cmd:approach(r2)} or {cmd:approach(modelfree)},
nor after {cmd:connectedness, nocovariance}.

{phang}
{opt frame(name)} names the results frame; by default the frame of the last
{cmd:connectedness} estimates or, failing that, {cmd:cn_results}.

{phang}
{opt test(test)} chooses the test of equal variances of the asset and the portfolio
whose p-value is shown: {cmd:fisher} (F test, as R's {cmd:var.test()}),
{cmd:bartlett}, {cmd:fk} (Fligner-Killeen, as {cmd:fligner.test()}), {cmd:levene}
(absolute deviations from the means) or {cmd:bf} (Brown-Forsythe: absolute deviations
from the medians, as {cmd:car::leveneTest()}). The p-values of all five tests are
stored in {cmd:r(pvalues)}.

{phang}
{opt cumulative(type)} sets the cumulative returns saved with the series: the
cumulative sum ({cmd:sum}) or the compounded return, prod(1 + r) - 1 ({cmd:prod}).

{phang}
{opt scale(#)} sets the number of periods per year used to annualize the return,
risk and Sharpe ratio. By default it follows the periodicity of the dates, as the R
package PerformanceAnalytics does: 252 for daily, 52 for weekly, 12 for monthly, 4 for
quarterly and 1 for yearly data. Without a date format, the measures are not
annualized ({cmd:scale(1)}).

{phang}
{opt saving(name[, replace])} names the frame for the series; the default frame
{cmd:cn_portfolio} is replaced each time.


{title:Remarks}

{pstd}
Hedging effectiveness is HE = 1 - Var(portfolio)/Var(asset) (Ederington 1979). The
annualized return is the geometric return, prod(1 + r)^(scale/n) - 1, the risk is the
standard deviation times sqrt(scale), and the Sharpe ratio is their ratio, with a
risk-free rate of zero.

{pstd}
The results reproduce those of R's {cmd:MinimumConnectednessPortfolio()},
{cmd:MultivariateHedgingPortfolio()}, {cmd:RiskParityPortfolio()},
{cmd:EquallyWeightedPortfolio()}, {cmd:BivariatePortfolio()} and {cmd:HedgeRatio()}
(ConnectednessApproach 1.0.4, PerformanceAnalytics 2.0.8, riskParityPortfolio 0.2.2)
to about 3e-14, and the risk-parity weights to the precision of R's solver (6e-9); see
{cmd:test/r/make_fixtures_portfolio.R}. Three choices differ from the R package. First,
each row's p-value compares that row's asset and portfolio; {cmd:HedgeRatio()} keeps
only the last pair's p-value in each row. Second, the option R calls "Brown-Forsythe"
tests equality of means, so {cmd:test(bf)} is the Brown-Forsythe test of variances,
the one R calls "Levene", and {cmd:test(levene)} centers on the means. Third, a pair
whose mean bivariate weight is exactly 0.5 is kept in the table, while R drops it
together with the diagonal.


{title:Examples}

{pstd}Minimum connectedness and bivariate portfolios from a static VAR{p_end}
{phang2}{cmd:. use cn_g2020, clear}{p_end}
{phang2}{cmd:. tsset t}{p_end}
{phang2}{cmd:. foreach v in EUR GBP CHF JPY {c -(}}{p_end}
{phang2}{cmd:.     generate double r_`v' = `v'/100}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. connectedness EUR GBP CHF JPY, lags(1) horizon(10) replace}{p_end}
{phang2}{cmd:. cnportfolio r_EUR r_GBP r_CHF r_JPY, method(mcp)}{p_end}
{phang2}{cmd:. cnportfolio r_EUR r_GBP r_CHF r_JPY, method(bivariate)}{p_end}

{pstd}Hedging with the conditional covariances of a DCC-GARCH model{p_end}
{phang2}{cmd:. connectedness EUR GBP CHF JPY, model(dccgarch) replace notable}{p_end}
{phang2}{cmd:. cnportfolio r_EUR r_GBP r_CHF r_JPY, method(mhp)}{p_end}
{phang2}{cmd:. cnportfolio r_EUR r_GBP r_CHF r_JPY, method(hedge) test(bf) saving(hr, replace)}{p_end}
{phang2}{cmd:. frame hr: tsline hr_1_2}{p_end}

{pstd}Time-varying minimum connectedness weights from a TVP-VAR{p_end}
{phang2}{cmd:. use cn_acg2020, clear}{p_end}
{phang2}{cmd:. tsset month}{p_end}
{phang2}{cmd:. foreach v in EUR GBP JPY CHF {c -(}}{p_end}
{phang2}{cmd:.     generate double r_`v' = `v'/100}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) replace notable}{p_end}
{phang2}{cmd:. cnportfolio r_EUR r_GBP r_JPY r_CHF, method(mcp) cumulative(prod)}{p_end}
{phang2}{cmd:. frame cn_portfolio: tsline w_1 w_2 w_3 w_4}{p_end}


{title:Stored results}

{pstd}
{cmd:cnportfolio} stores the following in {cmd:r()}:

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of periods{p_end}
{synopt:{cmd:r(k)}}number of assets{p_end}
{synopt:{cmd:r(scale)}}periods per year used for annualization{p_end}

{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:r(method)}}method{p_end}
{synopt:{cmd:r(matrix)}}{cmd:pci}, {cmd:cov}, {cmd:cor} or {cmd:none}{p_end}
{synopt:{cmd:r(test)}}test whose p-value is in {cmd:r(table)}{p_end}
{synopt:{cmd:r(cumulative)}}{cmd:sum} or {cmd:prod}{p_end}
{synopt:{cmd:r(long)}}{cmd:long} or {cmd:nolong}{p_end}
{synopt:{cmd:r(frame)}}results frame used{p_end}
{synopt:{cmd:r(saving)}}frame holding the series{p_end}
{synopt:{cmd:r(rows)}}row labels, for example {cmd:EUR/GBP}{p_end}

{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:r(table)}}one row per asset or pair: Mean, Std_Dev, p5 and p95 of the weights or ratios, HE, p_value, Return, Risk and SR{p_end}
{synopt:{cmd:r(pvalues)}}p-values of the five tests ({cmd:fisher}, {cmd:bartlett}, {cmd:fk}, {cmd:levene}, {cmd:bf}){p_end}

{pstd}
The frame of series holds the time variable (and the date variable of the results),
the weights {cmd:w_}{it:i} and the portfolio return {cmd:ret} and cumulative return
{cmd:cumret} ({cmd:mcp}, {cmd:rpp}, {cmd:ewp}); the hedge weights
{cmd:beta_}{it:i}{cmd:_}{it:j} and the hedged returns {cmd:ret_}{it:i} and
{cmd:cumret_}{it:i} ({cmd:mhp}); or, for each pair, the weight {cmd:w_}{it:i}{cmd:_}{it:j}
or ratio {cmd:hr_}{it:i}{cmd:_}{it:j} with {cmd:ret_}{it:i}{cmd:_}{it:j} and
{cmd:cumret_}{it:i}{cmd:_}{it:j} ({cmd:bivariate}, {cmd:hedge}).


{title:References}

{phang}
Broadstock, D. C., I. Chatziantoniou, and D. Gabauer. 2022. Minimum connectedness
portfolios and the market for green bonds: Advocating socially responsible investment
(SRI) activity. In {it:Applications in Energy Finance}, 217-253. Cham: Palgrave
Macmillan.

{phang}
Cocca, T., D. Gabauer, and S. Pomberger. 2024. Clean energy market connectedness and
investment strategies: New evidence from DCC-GARCH R2 decomposed connectedness
measures. {it:Energy Economics} 136: 107680.

{phang}
Ederington, L. H. 1979. The hedging performance of the new futures markets.
{it:Journal of Finance} 34: 157-170.

{phang}
Gabauer, D. 2025. {it:ConnectednessApproach: Connectedness Approach}. R package
version 1.0.4.

{phang}
Kroner, K. F., and V. K. Ng. 1998. Modeling asymmetric comovements of asset returns.
{it:Review of Financial Studies} 11: 817-844.

{phang}
Kroner, K. F., and J. Sultan. 1993. Time-varying distributions and dynamic hedging
with foreign currency futures. {it:Journal of Financial and Quantitative Analysis} 28:
535-551.

{phang}
Spinu, F. 2013. An algorithm for computing risk parity weights. SSRN working paper,
OMERS Capital Markets.


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

{p 7 14 2}
Help:  {helpb connectedness}, {helpb connectedness postestimation}, {helpb cngraph}, {helpb cnnetwork}{p_end}
