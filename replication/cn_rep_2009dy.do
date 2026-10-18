* Replication: Diebold and Yilmaz (2009), orthogonalized VAR connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2009DieboldYilmaz
version 16.1
use cn_dy2009, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
* Full-sample table (nlag = 2, nfore = 10, Cholesky decomposition)
connectedness `v', lags(2) horizon(10) orthogonalized frame(dy2009) replace
* Rolling windows of 200 weeks
connectedness `v', lags(2) horizon(10) orthogonalized window(200) datevar(date) ///
    frame(dy2009) replace notable
cngraph tci, yscale(range(40 60)) ylabel(40(5)60)
