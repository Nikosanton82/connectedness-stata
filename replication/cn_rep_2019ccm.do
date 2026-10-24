* Replication: Caloia, Cipollini and Muzzioli (2019), normalization schemes and net
* spillovers, a replication of Diebold and Yilmaz (2012). Static VAR(4) of the four
* range-based volatilities, horizon 10, under row normalization and under the scalar
* normalizations (maximum row sum, maximum column sum, spectral radius). Table 11 of
* the paper reports TO, FROM, NET and the TCI as shares; the tables here are in
* percent. Under the maximum row sum, bonds become larger net receivers than
* commodities, the ranking change the paper discusses.
version 16.1
use cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
connectedness `v', lags(4) horizon(10) frame(ccm_row) replace
foreach n in maxrow maxcol spectral {
    connectedness `v', lags(4) horizon(10) normalize(`n') frame(ccm_`n') replace
}
* rolling 200-day windows, as in Diebold and Yilmaz (2012), with the maximum row sum
connectedness `v', lags(4) horizon(10) window(200) normalize(maxrow) frame(ccm_roll) ///
    datevar(date) replace notable
cngraph net
