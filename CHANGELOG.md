# Changelog

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
