* Replication: Diebold and Yilmaz (2012), generalized VAR connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2012DieboldYilmaz
version 16.1
use cn_dy2012, clear
tsset t
* Full-sample table (nlag = 4, nfore = 10)
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) frame(dy2012) replace
* Rolling windows of 200 days
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) window(200) ///
    datevar(date) frame(dy2012) replace
cngraph tci, yscale(range(5 35))
cngraph to
cngraph from
cngraph net
cngraph npdc
cngraph pci
