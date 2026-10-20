* Replication: Chatziantoniou, Gabauer and Stenfors (2021), quantile connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2021ChatziantoniouGabauerStenfors
version 16.1
use cn_cgs2021, clear
tsset t
foreach q in 0.05 0.5 0.95 {
    local s = 100*`q'
    connectedness USD1Y EUR1Y GBP1Y JPY1Y, model(qvar) quantile(`q') lags(1) horizon(20) ///
        window(200) corrected datevar(date) frame(cgs2021_q`s') replace
    estimates store q`s'
}
cngraph tci, compare(q5 q95)
