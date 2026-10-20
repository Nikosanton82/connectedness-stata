* Replication: Naeem et al. (2024) and Balli et al. (2023), R2 decomposed connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2023BalliBalliDangGabauer
* and https://gabauerdavid.github.io/ConnectednessApproach/2023NaeemChatziantoniouGabauerKarim
version 16.1
use cn_dy2012, clear
tsset t
* contemporaneous only (Naeem et al. 2024)
connectedness SP500 R_10Y DJUBSCOM USDX, approach(r2) lags(0) window(200) ///
    datevar(date) frame(naeem) replace
* contemporaneous and lagged (Balli et al. 2023), three correlation measures
foreach m in pearson spearman kendall {
    connectedness SP500 R_10Y DJUBSCOM USDX, approach(r2) lags(1) window(200) ///
        cormethod(`m') datevar(date) frame(balli_`m') replace
}
connectedness, component(contemporaneous)
connectedness, component(lagged)
cngraph tci
