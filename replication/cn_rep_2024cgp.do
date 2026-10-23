* Replication: Cocca, Gabauer and Pomberger (2024), investment strategies for clean
* energy markets (data of the R package's cgp2024): the multivariate hedging portfolio
* from the DCC-GARCH conditional covariances, as MultivariateHedgingPortfolio(), and
* minimum connectedness and risk-parity portfolios for comparison. The portfolio
* functions take decimal returns, so the percentage returns are divided by 100.
version 16.1
use cn_cgp2024, clear
tsset t
foreach v in GRE SOL WIN GEO {
    generate double r_`v' = `v'/100
}
connectedness GRE SOL WIN GEO, model(dccgarch) horizon(10) datevar(date) frame(cgp2024_dcc) replace
cnportfolio r_GRE r_SOL r_WIN r_GEO, method(mhp) saving(cgp2024_mhp, replace)
cnportfolio r_GRE r_SOL r_WIN r_GEO, method(mcp) saving(cgp2024_mcp, replace)
cnportfolio r_GRE r_SOL r_WIN r_GEO, method(rpp) saving(cgp2024_rpp, replace)
cnportfolio r_GRE r_SOL r_WIN r_GEO, method(bivariate) saving(cgp2024_bpp, replace)
frame cgp2024_mcp: tsline cumret, title("Minimum connectedness portfolio, cumulative return")
