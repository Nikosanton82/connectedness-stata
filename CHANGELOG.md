# Changelog

## 0.8.2 (2026-10-24)
- Scalar normalizations of the generalized FEVD (Caloia, Cipollini and Muzzioli 2019):
  `normalize(maxrow)`, `normalize(maxcol)` and `normalize(spectral)` divide the
  unnormalized decomposition of Pesaran and Shin (1998) by its largest row sum, largest
  column sum or spectral radius instead of normalizing each row (`normalize(row)`, the
  default). They apply to `approach(time)` with VAR-type and TVP-VAR models;
  `estat group` requires the row normalization. `e(normalize)` records the choice.
- Adaptive LASSO, ridge and elastic-net VARs (Demirer, Diebold, Liu and Yilmaz 2018):
  option `adaptive` penalizes each slope with weight 1/|b_OLS|, passed to the elastic
  net as `glmnet`'s `penalty.factor` (rescaled to sum to the number of regressors).
  The glmnet port now supports penalty factors.
- Spillover asymmetry measure (Barunik, Kocenda and Vacha 2016, 2017): options
  `positive(varlist)` and `negative(varlist)` take the positive and negative series
  directly, such as realized semivariances, and estimate the three systems as
  `asymmetric` does. With either, `e(SAM)` holds the TCI of the positive system minus
  that of the negative system, and the results frame holds the per-period `sam`,
  `sam_to_i`, `sam_from_i` and `sam_net_i`.
- Validation: 9 R cases from `test/r/make_fixtures_082.R` (ConnectednessApproach 1.0.4,
  glmnet 4.1.8): the three normalizations of the Diebold and Yilmaz (2012) VAR and the
  maximum-row normalization of a TVP-VAR, computed from R's VAR() and TVPVAR() and
  passed to R's ConnectednessTable(), and five adaptive LASSO, ridge and elastic-net
  cases against glmnet with penalty factors. The R cases agree with R to 3e-13. TO,
  FROM, NET and TCI also match Table 11 of Caloia, Cipollini and Muzzioli (2019) under
  row and maximum-row normalization, to the three decimals published. SAM is checked
  against separate estimations of the positive and negative systems and against
  `asymmetric`. 214 of 214 certification tests pass in Stata 19
  (`docs/certification-log-v0.8.2.txt`).
- Replication do-file `cn_rep_2019ccm.do` (Caloia, Cipollini and Muzzioli 2019).

## 0.8.1 (2026-10-23)
- DCC-GARCH R2 decomposed connectedness (Cocca, Gabauer and Pomberger 2024):
  `model(dccgarch) approach(r2)` applies the R2 decomposition to the conditional
  correlation matrix of every period, as R's `R2Correlations()` applied to
  `rmgarch::rcor()`. As in R, a period whose shares for a variable sum to more than 1
  takes the previous period's shares; in the first period, where R stops with an
  error, Stata keeps the shares and reports a note. The pairwise connectedness index
  is stored in percent (R reports it on a scale of 0 to 2). The conditional
  covariances are stored in the results frame, so `cnportfolio` and
  `estat correlation` work after it.
- New `estat correlation [, partial]`: the conditional correlations of the covariance
  matrices stored in the results frame (VAR, rolling windows, TVP-VAR or DCC-GARCH),
  as R's `ConditionalCorrelation()`, or the partial correlations, as
  `PartialCorrelations()`. The series `cor_i_j` or `pcor_i_j` are saved in a frame
  (default `cn_corr`) and the average matrix is displayed and returned in `r(mean)`.
- Help: examples of DCC-GARCH R2 connectedness, conditional correlations and quantile
  frequency connectedness (Chatziantoniou, Abakah, Gabauer and Tiwari 2022), and the
  formulas of the DCC-GARCH R2 decomposition.
- Validation: 5 R cases from `test/r/make_fixtures_dccr2.R` (ConnectednessApproach
  1.0.4, rugarch 1.5.3, rmgarch 1.3.9): DCC-GARCH R2 connectedness of g2020 and
  cgp2024, and conditional and partial correlations after DCC-GARCH, static VAR and
  TVP-VAR. With rmgarch's estimates the results agree with R to 6e-14, and with
  Stata's own estimates to 2e-5. 187 of 187 certification tests pass in Stata 19
  (`docs/certification-log-v0.8.1.txt`).
- Replication do-files: `cn_rep_2022cagt.do` (Chatziantoniou, Abakah, Gabauer and
  Tiwari 2022: quantile frequency connectedness of green bonds, green equities,
  sustainable investment and clean energy at the 5%, 50% and 95% quantiles), and
  `cn_rep_2024cgp.do` now starts with the DCC-GARCH R2 connectedness of the paper.

## 0.8.0 (2026-10-23)
- New command `cnportfolio`: portfolio and hedging strategies from the matrices of
  each period, as the portfolio functions of the R package. `method(mcp)` gives the
  minimum connectedness portfolio (Broadstock, Chatziantoniou and Gabauer 2022),
  `method(mhp)` the multivariate hedging portfolio (Cocca, Gabauer and Pomberger 2024),
  `method(rpp)` the risk-parity portfolio, `method(ewp)` equal weights,
  `method(bivariate)` the bivariate portfolio weights of Kroner and Ng (1998) and
  `method(hedge)` the hedge ratios of Kroner and Sultan (1993). It reports the mean,
  standard deviation and 5% and 95% quantiles of the weights or ratios, hedging
  effectiveness, a test of equal variances (`test()`: F, Bartlett, Fligner-Killeen,
  Levene or Brown-Forsythe; all five p-values in `r(pvalues)`) and the annualized
  return, risk and Sharpe ratio, and saves the series in a frame. Options `matrix(pci|cov|cor)`,
  `cumulative(sum|prod)`, `nolong`, `scale()` and `saving()`.
- `connectedness` stores the covariance matrix of every period in the results frame,
  `cov_i_j` (i <= j): the residual covariance of the VAR, of each rolling window or of
  each TVP-VAR period, or the DCC-GARCH conditional covariance. New option
  `nocovariance`; the covariances are left out, with a note, when they would exceed
  maxvar, and are not stored with `approach(r2)` or `approach(modelfree)`.
  `e(covariance)` records whether they were stored.
- Differences from the R package, documented in `help cnportfolio`: each row's p-value
  compares that row's asset and portfolio (R's `HedgeRatio()` keeps only the last
  pair of each row); `test(bf)` is the Brown-Forsythe test of variances (R's
  `VarianceTest()` uses that name for a test of means and calls the variance test
  "Levene"), and `test(levene)` centers on the means; a pair whose mean bivariate
  weight is exactly 0.5 stays in the table.
- Validation: 17 R cases from `test/r/make_fixtures_portfolio.R` (ConnectednessApproach
  1.0.4, PerformanceAnalytics 2.0.8, riskParityPortfolio 0.2.2): all six methods with
  static VAR, TVP-VAR and DCC-GARCH matrices (PCI, covariance and correlation). The
  tables, p-values and series agree with R to 3e-14, and the risk-parity weights to
  6e-9, the precision of R's solver (`docs/r-fixtures-log-v0.8.0.txt`).
  175 of 175 certification tests pass in Stata 19 (`docs/certification-log-v0.8.0.txt`).
- Replication do-file `cn_rep_2024cgp.do` (Cocca, Gabauer and Pomberger 2024: multivariate
  hedging, minimum connectedness, risk-parity and bivariate portfolios from a DCC-GARCH
  model of clean energy markets).

## 0.7.2 (2026-10-20)
- `cngraph` takes the look of `cnnetwork`: light background, muted axes and gridlines,
  horizontal axis labels and a fixed palette (dark blue, mustard, teal, rust, grey,
  plum). NET and NPDC panels are areas in dark blue above zero and mustard below zero,
  matching net transmitters and receivers in `cnnetwork`; overlays from `compare()` or
  `asymmetric` take the next colors as well as dashed patterns. Date axes are labeled
  at the start of years in steps of 1, 2, 5 or 10 years, with fewer labels in panel
  graphs. Option `nostyle` keeps the graph scheme.
  151 of 151 certification tests pass in Stata 19 (`docs/certification-log-v0.7.2.txt`).

## 0.7.1 (2026-10-20)
- `cnnetwork` redesigned: curved edges whose width and shade grow with strength in
  six classes, arrowheads that end at the receiving node, net transmitters in dark
  blue and net receivers in mustard (drawn as true circles, slightly larger for larger
  NET), names inside the nodes for up to 12 short names and outside the circle
  otherwise, a legend below the circle, a subtitle with the period ("average over
  <first> - <last>") and a light background that does not depend on the graph scheme.
  Nodes, edges and labels shrink with the number of variables, so that large networks
  (for example 51 variables) stay readable. New options `straight`, `curvature()`,
  `colors()`, `labels()` and `nonote`.
- Dates: without `datevar()`, `connectedness` uses a numeric variable named `date` with
  a date format when the time variable has no date format, so that the results frame,
  the header, `cngraph` and `cnnetwork` show dates rather than observation numbers;
  `datevar()` naming the time variable turns this off.
- `cngraph` labels monthly and quarterly date axes by year, as it already did for daily
  dates.
- Data: `month` in `cn_acg2020`, `cn_gg2018` and `cn_gghm2022` now has the `%tm` format
  (values unchanged).
  151 of 151 certification tests pass in Stata 19 (`docs/certification-log-v0.7.1.txt`);
  all help-file examples and replication do-files run.

## 0.7.0 (2026-10-20)
- `model(dccgarch)`: DCC-GARCH(1,1) volatility connectedness (Gabauer 2020), as R's
  `ConnectednessApproach(..., model = "DCC-GARCH")`. GARCH(1,1) margins with a
  constant mean and normal errors and a DCC(1,1) correlation model are estimated in
  two steps by maximum likelihood (Newton-Raphson with analytic gradients), following
  the conventions of `rugarch` and `rmgarch::dccfit()`; the volatility FEVD of every
  period is a port of R's `VFEVD()`.
- Options `standardize` (R's `DCC_config = list(standardize = TRUE)`) and
  `dccparams()`, which uses given parameters instead of estimating them.
- Stored results `e(garch)`, `e(dcc)`, `e(ll)` and `e(converged)`.
- `cngraph, compare()` labels LASSO, ridge, elastic-net and DCC-GARCH results.
- Validation: 3 R cases (the Gabauer 2020 data with the default settings and with
  `standardize`, `corrected` and `horizon(20)`, and the data of Cocca, Gabauer and
  Pomberger 2024). With R's parameters the measures agree
  to 1.9e-12 and the log likelihood to 1e-15; Stata's own estimates agree with R's
  to 5.3e-6 (the precision of R's optimizer) and the measures to 8.0e-4.
  151 of 151 certification tests pass (`docs/certification-log-v0.7.0.txt`).
- Replication do-file `cn_rep_2020g.do`.

## 0.6.0 (2026-10-20)
- `model(lasso)`, `model(ridge)` and `model(elasticnet)`: penalized VARs as in R's
  `ElasticNetVAR()` (Gabauer, Gupta, Marfatia and Miller 2024). Each equation is
  fitted by a Mata port of glmnet's Gaussian coordinate descent (covariance updates),
  with the penalty chosen by 10-fold cross-validation at lambda.1se exactly as
  `cv.glmnet()` does, including its per-fold lambda paths and interpolation.
  `elasticnet` searches a grid of alphas (`nalpha()`, default 10 steps) or uses
  `alpha()`, and keeps the alpha with the smallest in-sample mean absolute error.
- Options `nfolds()`, `loss(mae|mse)`, `folds(random|interleaved)` and `seed()`.
- Validation: 7 new R cases with the folds that `folds(interleaved)` uses (LASSO,
  ridge, elastic net with a grid and with a fixed alpha, MSE loss, rolling windows,
  and a 51-variable LASSO VAR); largest difference from R 1.3e-11 (the 51-variable
  case) and below 1e-13 otherwise.
  132 of 132 certification tests pass (`docs/certification-log-v0.6.0.txt`).
- Certification script: the `seed()` reproducibility test keeps the first TCI in a
  scalar; Stata 19 rounds the local macro it used before, which made the exact
  comparison fail although the results are identical.

## 0.5.0 (2026-10-20)
- `approach(r2)`: R2 decomposed connectedness (Naeem et al. 2024), with contemporaneous
  and lagged components when `lags()` > 0 (Balli et al. 2023). The overall measures
  use the base names; the components get the suffixes `_c` and `_l`.
- `approach(modelfree)`: model-free connectedness from squared correlations
  (Gabauer, Chatziantoniou and Stenfors 2023); `estat group` works on these results.
- `cormethod(pearson|spearman|kendall)` for both approaches.
- `connectedness, component()`, `cngraph, component()` and `cnnetwork, component()`
  show the contemporaneous or lagged parts; `cngraph tci` draws all three TCIs.
- The pairwise influence index is stored as missing, instead of stopping the
  estimation, in periods where both shares of a pair are zero.
- Validation: 10 new R cases (Pearson, Spearman and Kendall; 0, 1 and 2 lags; static
  and rolling; 18 variables); largest relative difference from R below 1e-12;
  119 of 119 certification tests pass (`docs/certification-log-v0.5.0.txt`).

## 0.4.0 (2026-10-20)
- `model(qvar)`: quantile VAR connectedness (Chatziantoniou, Gabauer and Stenfors
  2021), with `quantile()` set to one value or one value per variable
  (quantile-on-quantile connectedness, Gabauer and Stenfors 2024). Each equation is a
  linear quantile regression solved by a Mata port of the Frisch-Newton interior-point
  method of R's `quantreg::rq()` (the "fn" method used by the R package).
- `model(lad)`: least absolute deviation VAR, matching `L1pack::lad()`.
- Both models work with rolling windows and with the time, frequency, joint and
  extended joint approaches.
- Header shows the model (QVAR or LAD-VAR) and the quantile; `e(quantile)` and
  `e(quantiles)` store it.
- Options `kappa1()`, `kappa2()`, `gamma()` and `prior()` are now rejected with
  models other than `tvpvar` instead of being ignored.
- Validation: 7 new R cases (static, rolling, quantile-on-quantile, LAD, QVAR
  frequency and QVAR extended joint); largest relative difference from R 6.0e-12;
  102 of 102 certification tests pass (`docs/certification-log-v0.4.0.txt`).
- Fixed: `cngraph` and `cnnetwork` printed a spurious "invalid syntax" message
  (without failing) whenever `band()` was not specified. Found when testing the
  installed package in Stata 19.
- Help: the joint-connectedness example of `connectedness` now specifies `replace`,
  so that the examples run in sequence.
- `cngraph`: the legend of single-panel graphs (the TCI with bands, comparisons or
  asymmetric parts) is placed below the graph, so that the plot keeps its width
  under the default scheme of Stata 18 and 19.

## 0.3.0 (2026-10-20)
- `approach()` option: `frequency` (Baruník and Křehlík 2018; with `model(tvpvar)`,
  Chatziantoniou, Gabauer and Gupta 2023), `joint` (Lastrapes and Wiesen 2021) and
  `extjoint` (Balcilar, Gabauer and Umar 2021), for static, rolling and TVP-VAR models.
- Frequency bands set in periods with `bands()` or in radians with `partition()`;
  `scenario(abs|wth)` for absolute or within-band measures. The results frame holds the
  total and, for each band, the same measures with the suffix `_f<b>`.
- `connectedness, band(#)` and `estat table, band(#)` show a band's table; the header
  lists the band TCIs.
- `cngraph`: `band()` option; `cngraph tci` draws the total and every band after
  `approach(frequency)`. `cnnetwork`: `band()` option.
- `estat group` requires `approach(time)`, as in the R package.
- Validation: 13 new R cases from `test/r/make_fixtures_ext.R` (ConnectednessApproach
  1.0.4), covering generalized and orthogonalized, ABS and WTH, static, rolling and
  TVP-VAR, and corrected frequency results; NumPy reference port in
  `test/python/ca_reference_ext.py`. Largest relative difference from R: 1.4e-12 for
  static and rolling models, 2.6e-8 for the TVP-VAR (Minnesota prior, k = 12);
  90 of 90 certification tests pass (`docs/certification-log-v0.3.0.txt`).
- Replication do-files for Baruník and Křehlík (2018), Lastrapes and Wiesen (2021),
  Balcilar, Gabauer and Umar (2021) and Chatziantoniou, Gabauer and Gupta (2023).

## 0.2.1 (2026-10-18), prepared for Stata Journal submission
- `notable` (in `connectedness` and `estat group`) now suppresses only the table; the
  header and the TCIs are still shown.
- Header reports the Minnesota prior tightness.
- `cngraph, compare()`: the legend labels the current results by model (for example,
  TVP-VAR) instead of "current".
- `cnnetwork`: transmitters and receivers are told apart by filled and hollow markers,
  and colors follow the graph scheme (grayscale with `set scheme sj`).
- Help files: American spelling, author and "Also see" sections in journal format.

## 0.2.0 (2026-10-18), development
- `asymmetric` option: estimates the systems of all, positive and negative changes
  (Adekoya et al. 2022) into three frames; `connectedness, part()` shows each table.
- `estat group`: aggregated, conditional, inclusive, exclusive, internal and external
  connectedness (Gabauer and Gupta 2018; Stenfors et al. 2022; Chatziantoniou et al.
  2023), with group TCIs; results in their own frame.
- `cnnetwork`: network plots of net pairwise directional connectedness or the
  pairwise connectedness index, averaged or for one period.
- `cngraph`: `frame()` option, `gtci` and `tciext` measures, overlays for `compare()`
  and asymmetric results on every measure, year labels for daily dates.
- The full FEVD is now stored in the results frame by default (`nofevd` to skip);
  `savefevd` is accepted and has no effect.
- Header shows dates when `datevar()` is given.
- Validation extended to 14 R replication cases, 8 group-measure cases and 3
  asymmetric systems; reference port reproduces 8 more published tables.

## 0.1.0 (2026-10-18), development
- `connectedness`: static and rolling-window VAR (OLS) and TVP-VAR (Kalman filter with
  forgetting factors; Bayes, Minnesota and uninformative priors); generalized and
  orthogonalized FEVD; TCI, corrected TCI, TO, FROM, NET, NPT, NPDC, PCI, influence.
- Results frame with the time-varying measures; optional full FEVD (`savefevd`).
- `estat table`, `cngraph` (tci, to, from, net, npt, npdc, pci, inf; `compare()` for TCI).
- 18 datasets from the ConnectednessApproach R package.
- Certification suite against fixtures from a validated Python port of the R code
  and from R itself (ConnectednessApproach 1.0.4, R 4.5.2): 31 tests, all passed on
  Stata 18 (largest relative difference 9e-9, TVP-VAR with 3,915 periods). Logs in docs/.
