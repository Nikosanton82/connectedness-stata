connectedness: VAR-type, TVP-VAR and DCC-GARCH connectedness measures for Stata
================================================================================
Version 0.7.0, 20 October 2026 (the Stata Journal submission describes version 0.2.1)
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
  data/           18 datasets from published applications (cn_*.dta); after installing,
                  type  net get connectedness  to copy them to the current directory
  replication/    do-files reproducing the replication pages of the R package
  test/           certification suite (run_all.do), reference fixtures and the
                  Python and R scripts that generate them

Commands
  connectedness   estimation (static and rolling VAR, QVAR and LAD-VAR, LASSO, ridge
                  and elastic-net VARs, TVP-VAR, DCC-GARCH; time,
                  frequency, joint, extended joint, R2 decomposed and model-free
                  connectedness; asymmetric option)
  estat table     redisplay the connectedness table
  estat group     aggregated, conditional, inclusive, exclusive, internal and
                  external connectedness
  cngraph         graphs of TCI, TO, FROM, NET, NPT, NPDC, PCI, influence, group TCIs
  cnnetwork       network plots of pairwise connectedness

Data
  The datasets are those distributed with the ConnectednessApproach R package
  (D. Gabauer, version 1.0.4, CRAN, GPL-3) and are redistributed under the same
  license. Each file records its source paper in its data label. Variables whose
  R names were not valid Stata names were renamed (for example US.FPU -> US_FPU);
  the original names are kept as variable labels.

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
  installed with ConnectednessApproach.
