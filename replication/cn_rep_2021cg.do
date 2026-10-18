* Replication: Chatziantoniou and Gabauer (2021), corrected TCI and pairwise connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2021ChatziantoniouGabauer
version 16.1
use cn_cg2021, clear
tsset t
connectedness Austria-Spain, model(tvpvar) lags(1) horizon(10) window(200) corrected ///
    datevar(date) frame(cg2021) replace
cngraph tci
cngraph net
cngraph pci, select(Germany)
