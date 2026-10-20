# connectedness

Stata package for the connectedness approach: Diebold and Yilmaz (2009, 2012) VAR
connectedness, rolling-window measures, and the TVP-VAR connectedness measures of
Antonakakis, Chatziantoniou and Gabauer (2020), with the corrected total, pairwise and
influence indices of Chatziantoniou and Gabauer (2021) and Gabauer (2021). Since
version 0.3.0 it also computes frequency connectedness (Baruník and Křehlík 2018),
the joint spillover index (Lastrapes and Wiesen 2021) and extended joint
connectedness (Balcilar, Gabauer and Umar 2021), and since version 0.4.0 quantile
connectedness (Chatziantoniou, Gabauer and Stenfors 2021) from quantile and LAD VARs.

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

* frequency connectedness: bands of 1-5 and more than 5 months
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(100) window(200) approach(frequency) bands(5)
cngraph tci
```

Status: development (version 0.4.0). The Stata Journal article describes version 0.2.1:
static and rolling VAR, TVP-VAR, corrected TCI, PCI, influence, and the aggregated,
conditional and decomposed group measures. Version 0.3.0 adds frequency, joint and
extended joint connectedness, and version 0.4.0 quantile and LAD VARs. Later releases
add LASSO-type, R2 and DCC-GARCH connectedness.

Requires Stata 16.1 or newer. License: GPL-3.
