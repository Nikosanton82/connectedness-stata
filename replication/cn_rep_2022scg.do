* Replication: Stenfors, Chatziantoniou and Gabauer (2022), conditional and aggregated
* connectedness across European yield curves
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2022StenforsChatziantoniouGabauer
version 16.1
use cn_cgs2022, clear
tsset t
connectedness X2YCHF-X10YSEK, model(tvpvar) lags(1) horizon(10) window(250) corrected ///
    datevar(date) frame(cgs2022) replace
cngraph tci, yscale(range(70 100))
* Conditional connectedness: NOK curve, then 2-year rates across currencies
estat group, type(conditional) groups(X2YNOK X5YNOK X10YNOK) frame(cgs_nok) replace
cngraph tci, frame(cgs_nok)
estat group, type(conditional) groups(X2YCHF X2YDKK X2YEUR X2YGBP X2YNOK X2YSEK) ///
    frame(cgs_2y) replace
cngraph net, frame(cgs_2y)
* Aggregated connectedness by currency and by maturity
estat group, type(aggregated) frame(cgs_cur) replace ///
    groups("CHF: X2YCHF X5YCHF X10YCHF" "DKK: X2YDKK X5YDKK X10YDKK" ///
    "EUR: X2YEUR X5YEUR X10YEUR" "GBP: X2YGBP X5YGBP X10YGBP" ///
    "NOK: X2YNOK X5YNOK X10YNOK" "SEK: X2YSEK X5YSEK X10YSEK")
cngraph net, frame(cgs_cur)
cnnetwork, frame(cgs_cur)
estat group, type(aggregated) frame(cgs_mat) replace ///
    groups("Y2: X2YCHF X2YDKK X2YEUR X2YGBP X2YNOK X2YSEK" ///
    "Y5: X5YCHF X5YDKK X5YEUR X5YGBP X5YNOK X5YSEK" ///
    "Y10: X10YCHF X10YDKK X10YEUR X10YGBP X10YNOK X10YSEK")
cngraph net, frame(cgs_mat)
