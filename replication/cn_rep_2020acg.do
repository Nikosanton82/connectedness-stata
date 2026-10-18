* Replication: Antonakakis, Chatziantoniou and Gabauer (2020), TVP-VAR connectedness
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2020AntonakakisChatziantoniouGabauer
version 16.1
use cn_acg2020, clear
tsset month
* Rolling-window VARs of 50, 100 and 200 months for comparison
foreach w in 50 100 200 {
    connectedness EUR GBP JPY CHF, lags(1) horizon(12) window(`w') ///
        frame(acg_roll`w') replace notable
    estimates store roll`w'
}
* TVP-VAR: kappa1 = 0.99, kappa2 = 0.96, Bayes prior on the first 200 months
connectedness EUR GBP JPY CHF, model(tvpvar) lags(1) horizon(12) window(200) ///
    kappa1(0.99) kappa2(0.96) prior(bayes) frame(acg_tvp) replace
cngraph tci, compare(roll50 roll100 roll200) yscale(range(20 80))
cngraph net, yscale(range(-20 20))
cngraph npdc, yscale(range(-20 20))
