* Replication: Adekoya, Akinseye, Antonakakis, Chatziantoniou, Gabauer and Oliyide (2022),
* asymmetric TVP-VAR connectedness of crude oil and Islamic sectoral stocks
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2022Adekoya
version 16.1
use cn_aaacgo2022, clear
tsset t
foreach v of varlist Brent-UTI {
    generate double r_`v' = 100*(`v' - L.`v')/L.`v'
    label variable r_`v' "`v' (percentage change)"
}
* window(200) mirrors the R replication; it is not used with the Minnesota prior
connectedness r_*, model(tvpvar) lags(1) horizon(10) window(200) prior(minnesota) ///
    gamma(0.1) asymmetric datevar(date) frame(aaacgo) replace
connectedness, part(positive)
connectedness, part(negative)
cngraph tci, yscale(range(50 100))
cngraph net, yscale(range(-50 50))
cngraph npdc, select(r_Brent) yscale(range(-20 20))
cngraph pci, select(r_Brent)
