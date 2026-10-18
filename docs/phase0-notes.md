# Phase 0 notes (18 October 2026)

Findings from reading the ConnectednessApproach R source (GitHub, version 1.0.4,
licence GPL-3). These settle the open questions in the blueprint.

- **Horizon.** `FEVD()` sums the moving-average matrices for h = 0, ..., nfore - 1.
- **Rolling windows.** `ConnectednessApproach()` sets `window.size = window.size - nlag`,
  so each rolling VAR uses window - nlag rows of data and window - 2*nlag regression rows.
- **TVP-VAR and `window.size`.** Not a rolling window. It is passed (after subtracting
  nlag) to `BayesPrior()`, which subtracts nlag again: the prior is an OLS VAR with
  intercept on window - 2*nlag regression observations. With `window.size = NULL` the
  whole sample is used. The TVP-VAR then runs over the full sample and measures are
  reported for every period.
- **TVP-VAR details.** Data are demeaned; the first error covariance is the prior's;
  for periods 2 to nlag+1 the covariance update uses 0.1 x_t x_t'; unstable draws
  (companion eigenvalue modulus > 1) keep the previous coefficients and shrink the
  state by 0.99. Matrix inverses use MASS::ginv.
- **VAR.** OLS with intercept; residual covariance divided by n (no degrees-of-freedom
  correction).
- **Sign convention.** R stores NPDC[i,j] = CT[i,j] - CT[j,i] and plots it as "j - i".
  The Stata results frame stores npdc_i_j (i < j) as the net transmission from i to j,
  which is the same quantity R plots under the title "i - j".
- **TABLE.** TO row ends with sum(TO); NET row ends with the TCI; NPT counts, for each
  variable, the series to which it is a net transmitter.
- **Licence.** GPL-3, so the Stata port is released under GPL-3. The 18 datasets ship
  with the CRAN package and are redistributed under the same licence.

Validation of the reference port (test/python/ca_reference.py) against the tables
printed on the replication pages: identical to every printed digit for Diebold and
Yilmaz (2009) static, Diebold and Yilmaz (2012) static, the VAR(4) error covariance
of dy2012 (six decimals), Gabauer and Gupta (2018) TVP-VAR and Chatziantoniou and
Gabauer (2021) TVP-VAR with corrected TCI.
