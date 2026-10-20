# connectedness

Stata package for the connectedness approach: Diebold and Yilmaz (2009, 2012) VAR
connectedness, rolling-window measures, and the TVP-VAR connectedness measures of
Antonakakis, Chatziantoniou and Gabauer (2020), with the corrected total, pairwise and
influence indices of Chatziantoniou and Gabauer (2021) and Gabauer (2021). Since
version 0.3.0 it also computes frequency connectedness (Baruník and Křehlík 2018),
the joint spillover index (Lastrapes and Wiesen 2021) and extended joint
connectedness (Balcilar, Gabauer and Umar 2021), since version 0.4.0 quantile
connectedness (Chatziantoniou, Gabauer and Stenfors 2021) from quantile and LAD VARs,
since version 0.5.0 R2 decomposed connectedness (Naeem et al. 2024; Balli et al.
2023) and model-free connectedness (Gabauer, Chatziantoniou and Stenfors 2023), and
since version 0.6.0 LASSO, ridge and elastic-net VARs (Gabauer et al. 2024), and
since version 0.7.0 DCC-GARCH volatility connectedness (Gabauer 2020).

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

Status: development (version 0.7.1). The Stata Journal article describes version 0.2.1:
static and rolling VAR, TVP-VAR, corrected TCI, PCI, influence, and the aggregated,
conditional and decomposed group measures. Version 0.3.0 adds frequency, joint and
extended joint connectedness, version 0.4.0 quantile and LAD VARs, and version 0.5.0
R2 decomposed and model-free connectedness, version 0.6.0 LASSO, ridge and
elastic-net VARs, and version 0.7.0 DCC-GARCH volatility connectedness.

Requires Stata 16.1 or newer. License: GPL-3.
