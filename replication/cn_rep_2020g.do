* Replication: Gabauer (2020), DCC-GARCH volatility impulse responses and volatility
* connectedness of the EUR, GBP, CHF and JPY exchange rates (data of the R package's
* g2020), as ConnectednessApproach(g2020, model = "DCC-GARCH", nfore = 10)
version 16.1
use cn_g2020, clear
tsset t
connectedness EUR GBP CHF JPY, model(dccgarch) horizon(10) datevar(date) frame(g2020_dcc) replace
matrix list e(garch)
matrix list e(dcc)
cngraph tci
cngraph net
* standardized volatility impulse responses
connectedness EUR GBP CHF JPY, model(dccgarch) horizon(10) standardize datevar(date) ///
    frame(g2020_dcc_std) replace
