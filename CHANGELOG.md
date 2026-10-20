# Changelog

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
