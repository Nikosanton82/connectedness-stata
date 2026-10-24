* Replication: Chatziantoniou, Abakah, Gabauer and Tiwari (2022), quantile
* time-frequency connectedness of green bond, green equity, sustainable investment
* and clean energy markets (data of the R package's gcat2022). Rolling quantile VAR(1)
* with a 200-day window, horizon 100 and two bands: up to 5 days (short run) and
* longer than 5 days (long run), at both tails and at the median, as
* ConnectednessApproach(..., model="QVAR", connectedness="Frequency").
version 16.1
use cn_gcat2022, clear
tsset t
local v S_P_GREEN_BOND_INDEX-S_P_GLOBAL_CLEAN_ENERGY
foreach q in 05 95 50 {
    connectedness `v', model(qvar) quantile(0.`q') lags(1) horizon(100) window(200) ///
        approach(frequency) bands(5) datevar(date) frame(gcat2022_q`q') replace
    connectedness, band(1)
    connectedness, band(2)
    estimates store q`q'
}
cngraph tci, compare(q05 q95)
cngraph net, band(1)
