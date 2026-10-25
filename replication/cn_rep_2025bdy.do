* Replication in the spirit of Buchwalter, Diebold and Yilmaz (2026, Journal of
* Econometrics), clustered network connectedness of global equity markets. The
* paper's weekly data (16 markets,
* 2002-2021) come from WRDS and are not distributed; this do-file applies the same
* design to the weekly equity returns of Diebold and Yilmaz (2009), 19 markets in three
* regions: a VAR(3) estimated by the adaptive elastic net (alpha 0.5, weights from
* OLS, 10-fold cross-validation), horizon 12, shocks correlated within regions and
* orthogonal across them, averaged over the 6 orderings of the regions; then rolling
* 104-week windows with an OLS VAR(2).
version 16.1
use cn_dy2009, clear
tsset t
local reg `""Americas: US ARG BRA CHL MEX" "Europe: UK FRA GER TUR" "AsiaPacific: HKG JPN AUS IDN KOR MYS PHL SGP TAI THA""'
quietly ds t date, not
local v `r(varlist)'
connectedness `v', model(elasticnet) alpha(0.5) adaptive lags(3) horizon(12) seed(2025) ///
    clusters(`reg') frame(bdy_static) replace notable
matrix list e(cl_shares), format(%8.2f)
matrix list e(cl_contagion), format(%8.2f)
connectedness `v', lags(2) horizon(12) window(104) clusters(`reg') datevar(date) ///
    frame(bdy_roll) replace notable
frame bdy_roll: tsline cl_contag_1 cl_contag_2 cl_contag_3, ///
    title("Contagion received by region") legend(order(1 "Americas" 2 "Europe" 3 "Asia-Pacific"))
