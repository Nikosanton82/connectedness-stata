# Changelog

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
