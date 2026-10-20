* Replication: Gabauer, Chatziantoniou and Stenfors (2023), model-free connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2023GabauerChatziantoniouStenfors
version 16.1
use cn_dy2012, clear
tsset t
connectedness SP500 R_10Y DJUBSCOM USDX, approach(modelfree) frame(mf) replace
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) frame(mb) replace
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) window(200) ///
    datevar(date) frame(dmb) replace notable
estimates store dmb
foreach m in pearson spearman kendall {
    connectedness SP500 R_10Y DJUBSCOM USDX, approach(modelfree) window(200) ///
        cormethod(`m') datevar(date) frame(mf_`m') replace notable
    estimates store mf_`m'
}
cngraph tci, compare(dmb)
