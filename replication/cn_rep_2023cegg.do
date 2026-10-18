* Replication: Chatziantoniou, Elsayed, Gabauer and Gozgor (2023), decomposed and partial
* connectedness of oil price shocks and exchange rates
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2022ChatziantoniouElsayedGabauerGozgor
version 16.1
use cn_cegg2022, clear
tsset t
connectedness Demand_Shock-Japan, model(tvpvar) lags(1) horizon(20) prior(minnesota) ///
    gamma(0.1) corrected datevar(date) frame(cegg) replace
local G `""Oil: Demand_Shock Supply_Shock Risk_Shock" "Exporters: Brazil Canada Norway Russia" "Importers: China Europe India Japan""'
estat group, type(internal) groups(`G') frame(cegg_int) replace
estat group, type(external) groups(`G') frame(cegg_ext) replace
estat group, type(inclusive) groups(Demand_Shock Supply_Shock Risk_Shock) frame(cegg_inc) replace
estat group, type(exclusive) groups(Demand_Shock Supply_Shock Risk_Shock) frame(cegg_exc) replace
cngraph tci
cngraph tci, frame(cegg_inc)
cngraph tci, frame(cegg_exc)
cngraph net, frame(cegg_inc)
