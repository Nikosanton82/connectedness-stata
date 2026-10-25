connectedness: VAR-type, TVP-VAR and DCC-GARCH connectedness measures for Stata
================================================================================
Version 0.9.4, 25 October 2026 (the Stata Journal submission describes version 0.2.1)
Author:  Nikolaos Antonakakis, University of Nicosia, UNIC Athens, antonakakis.n@unic.ac.cy
License: GNU General Public License version 3 (see LICENSE)

Requirements
  Stata 16.1 or newer. No community-contributed packages are needed.

Installation
  From a local copy of this folder:
      net install connectedness, from("<path to this folder>") replace
  From GitHub (once the repository is public):
      net install connectedness, from("https://raw.githubusercontent.com/Nikosanton82/connectedness-stata/main") replace

Contents
  ado/            commands, help files and the Mata engine (cn_core.mata)
  data/           19 datasets from published applications (cn_*.dta); after installing,
                  type  net get connectedness  to copy them to the current directory
  replication/    do-files reproducing the replication pages of the R package
  test/           certification suite (run_all.do), reference fixtures and the
                  Python and R scripts that generate them

Commands
  connectedness   estimation (static and rolling VAR, QVAR and LAD-VAR, LASSO, ridge
                  and elastic-net VARs, TVP-VAR, DCC-GARCH; time,
                  frequency, joint, extended joint, R2 decomposed (also of the
                  DCC-GARCH correlations) and model-free connectedness; row or
                  scalar normalization; clustered identification; adaptive
                  elastic nets; asymmetric connectedness and the spillover
                  asymmetry measure)
  estat table     redisplay the connectedness table
  estat group     aggregated, conditional, inclusive, exclusive, internal and
                  external connectedness
  estat correlation  conditional or partial correlations of each period
  estat bootstrap    bootstrap bands (also of the spillover asymmetry measure),
                     event probabilities and detection of changes
  cngraph         graphs of TCI, TO, FROM, NET, NPT, NPDC, PCI, influence, group TCIs
  cnnetwork       network plots of pairwise connectedness
  cnportfolio     portfolio weights, hedge ratios and hedging effectiveness (minimum
                  connectedness, multivariate hedging, risk-parity, equally weighted
                  and bivariate portfolios; hedge ratios)

Data
  The datasets are those distributed with the ConnectednessApproach R package
  (D. Gabauer, version 1.0.4, CRAN, GPL-3) and are redistributed under the same
  license. Each file records its source paper in its data label. Variables whose
  R names were not valid Stata names were renamed (for example US.FPU -> US_FPU);
  the original names are kept as variable labels. cn_gkn2024.dta holds the daily
  volatilities of Diebold and Yilmaz (2009) from the replication files of
  Greenwood-Nimmo, Kocenda and Nguyen (2024), Mendeley Data, doi:10.17632/cvp7v9nxkf.1,
  CC BY 4.0.

Reproducing the checks
  In Stata:  cd "<this folder>/test"
             do run_all.do
  The log is written to test/logs/cert_connectedness.log. test/run_install.do
  installs the package into a temporary directory and runs every replication. To check against R itself,
  run test/r/make_fixtures.R first (on Windows, test/r/run_r_fixtures.bat; edit the R
  path inside it if needed; the script installs ConnectednessApproach if missing);
  the certification script then also compares with the R_* fixtures. The fixtures of
  the extensions (version 0.3.0 onwards) come from test/r/make_fixtures_ext.R; its
  DCC-GARCH cases (version 0.7.0) use the rugarch and rmgarch packages, which are
  installed with ConnectednessApproach. The fixtures of cnportfolio (version 0.8.0)
  come from test/r/make_fixtures_portfolio.R (test/r/run_r_portfolio.bat), which uses
  the PerformanceAnalytics and riskParityPortfolio packages, also installed with
  ConnectednessApproach. The fixtures of DCC-GARCH R2 decomposed connectedness and
  estat correlation (version 0.8.1) come from test/r/make_fixtures_dccr2.R
  (test/r/run_r_dccr2.bat), and those of the normalizations and adaptive elastic nets
  (version 0.8.2) from test/r/make_fixtures_082.R (test/r/run_r_082.bat). Clustered
  connectedness (version 0.9.0) is checked against an independent R implementation in
  test/r/make_fixtures_090.R (test/r/run_r_090.bat), and estat bootstrap (version 0.9.1)
  against one in test/r/make_fixtures_091.R (test/r/run_r_091.bat), which shares its
  resampled rows with Stata; its pairwise and NPDC bands and stat() (version 0.9.2)
  against test/r/make_fixtures_092.R (test/r/run_r_092.bat); the cluster table and the
  within/cross split of the TCI (version 0.9.3) against test/r/make_fixtures_093.R
  (test/r/run_r_093.bat); the bootstrap of the spillover asymmetry measure (version
  0.9.4) against test/r/make_fixtures_094.R (test/r/run_r_094.bat), with a Monte Carlo
  study of its coverage in test/r/sam_coverage_094.R.
