* Replication: Barunik and Krehlik (2018), frequency connectedness (static VAR)
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2018BarunikKrehlik
* The R page uses simulated data; here the frequency bands are applied to the
* Diebold and Yilmaz (2012) volatilities instead.
version 16.1
use cn_dy2012, clear
tsset t
connectedness SP500 R_10Y DJUBSCOM USDX, lags(2) horizon(100) approach(frequency) ///
    bands(4) frame(bk2018) replace
connectedness, band(1)
connectedness, band(2)
