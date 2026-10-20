* Replication: Chatziantoniou, Gabauer and Gupta (2023), TVP-VAR frequency connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2022GabauerChatziantoniouGupta
version 16.1
use cn_cgg2022, clear
tsset t
connectedness Europe-OPEC, model(tvpvar) lags(1) horizon(100) window(200) ///
    approach(frequency) bands(5) datevar(date) frame(cgg2022) replace
connectedness, band(1)
connectedness, band(2)
cngraph tci
cngraph net, band(1)
cngraph net, band(2)
