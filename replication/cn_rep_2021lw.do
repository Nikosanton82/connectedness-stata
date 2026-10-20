* Replication: Lastrapes and Wiesen (2021), joint spillover index
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2021LastrapesWiesen
version 16.1
use cn_lw2021, clear
tsset t
connectedness energy-utilities, lags(2) horizon(30) approach(joint) frame(lw2021) replace
connectedness energy-utilities, lags(2) horizon(30) window(200) approach(joint) ///
    datevar(date) frame(lw2021_roll) replace notable
cngraph tci, yscale(range(40 100))
