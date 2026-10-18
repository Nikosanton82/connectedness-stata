# connectedness

Stata package for the connectedness approach: Diebold and Yilmaz (2009, 2012) VAR
connectedness, rolling-window measures, and the TVP-VAR connectedness measures of
Antonakakis, Chatziantoniou and Gabauer (2020), with the corrected total, pairwise and
influence indices of Chatziantoniou and Gabauer (2021) and Gabauer (2021).

The numerical routines port the [ConnectednessApproach](https://gabauerdavid.github.io/ConnectednessApproach/)
R package by David Gabauer (GPL-3) and reproduce its results to numerical precision;
see `test/`.

```stata
net install connectedness, from("https://raw.githubusercontent.com/Nikosanton82/connectedness-stata/main") replace
net get connectedness          // copies the example datasets to the current directory
use cn_acg2020, clear
tsset month
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96)
cngraph tci
cngraph net
```

Status: development (version 0.2.1). Scope of the first release and of the Stata Journal
article: static and rolling VAR, TVP-VAR, corrected TCI, PCI, influence, and the
aggregated, conditional and decomposed group measures. Later releases add quantile,
LASSO-type, frequency, joint, R2 and DCC-GARCH connectedness.

Requires Stata 16.1 or newer. License: GPL-3.
