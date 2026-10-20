* Replication: Balcilar, Gabauer and Umar (2021), TVP-VAR extended joint connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2021BalcilarGabauerUmar
version 16.1
use cn_bgu2021, clear
tsset t
connectedness CrudeOil-Wheat, model(tvpvar) lags(1) horizon(20) prior(minnesota) gamma(0.1) ///
    approach(extjoint) datevar(date) frame(bgu2021) replace
cngraph tci
cngraph net
cnnetwork
