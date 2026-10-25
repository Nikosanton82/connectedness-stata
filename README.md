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
since version 0.6.0 LASSO, ridge and elastic-net VARs (Gabauer et al. 2024),
since version 0.7.0 DCC-GARCH volatility connectedness (Gabauer 2020), since
version 0.8.0 portfolio and hedging strategies with `cnportfolio`: minimum
connectedness, multivariate hedging, risk-parity, equally weighted and bivariate
portfolios and hedge ratios, and since version 0.8.1 DCC-GARCH R2 decomposed
connectedness (Cocca, Gabauer and Pomberger 2024) and conditional and partial
correlations (`estat correlation`), and since version 0.8.2 scalar normalizations
(Caloia, Cipollini and Muzzioli 2019), adaptive LASSO, ridge and elastic-net VARs and
the spillover asymmetry measure of good and bad volatility (Baruník, Kočenda and Vácha
2016), since version 0.9.0 clustered network connectedness (Buchwalter, Diebold and
Yilmaz 2026, *Journal of Econometrics*), and since version 0.9.1 bootstrap inference and event tests
(Greenwood-Nimmo, Kočenda and Nguyen 2024), extended in version 0.9.2 to the pairwise
and net pairwise (NPDC) measures and in version 0.9.4 to the spillover asymmetry measure,
and since version 0.9.5 control clusters, the ordering-sensitivity diagnostic and
simultaneous clustered orthogonalization (Buchwalter, Diebold and Yilmaz 2026b, 2026c).

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

* minimum connectedness portfolio of the four currencies (decimal returns)
foreach v in EUR GBP JPY CHF {
    generate double r_`v' = `v'/100
}
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) replace
cnportfolio r_EUR r_GBP r_JPY r_CHF, method(mcp)
```

Status: development (version 0.9.5). The Stata Journal article describes version 0.2.1:
static and rolling VAR, TVP-VAR, corrected TCI, PCI, influence, and the aggregated,
conditional and decomposed group measures. Version 0.3.0 adds frequency, joint and
extended joint connectedness, version 0.4.0 quantile and LAD VARs, and version 0.5.0
R2 decomposed and model-free connectedness, version 0.6.0 LASSO, ridge and
elastic-net VARs, version 0.7.0 DCC-GARCH volatility connectedness, version 0.8.0
portfolio and hedging strategies (`cnportfolio`), version 0.8.1 DCC-GARCH R2
decomposed connectedness and `estat correlation`, and version 0.8.2 scalar
normalizations, adaptive elastic nets and the spillover asymmetry measure, version
0.9.0 clustered network connectedness, version 0.9.1 bootstrap inference
(`estat bootstrap`), and version 0.9.2 bootstrap bands for the pairwise and NPDC
measures and event tests for any measure, and version 0.9.3 the cluster table and the
within/cross split of the TCI of the published clustered connectedness paper, and
version 0.9.4 bootstrap bands for the spillover asymmetry measure (`estat bootstrap, sam`),
and version 0.9.5 a control cluster (`clcontrols()`), the ordering-sensitivity
diagnostic (`estat clorder`) and simultaneous clustered orthogonalization
(`clorder(simultaneous)`).

Requires Stata 16.1 or newer. License: GPL-3.
