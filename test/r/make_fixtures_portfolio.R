# Reference fixtures for cnportfolio (version 0.8.0), computed with the portfolio
# functions of the ConnectednessApproach R package (MinimumConnectednessPortfolio,
# MultivariateHedgingPortfolio, RiskParityPortfolio, EquallyWeightedPortfolio,
# BivariatePortfolio, HedgeRatio) and PerformanceAnalytics.
#
# Run from the test/r folder:   Rscript make_fixtures_portfolio.R
# Requires ConnectednessApproach 1.0.4 with its imports (PerformanceAnalytics,
# riskParityPortfolio, car), rugarch and rmgarch.
#
# Output, per case <c>, in ../fixtures:
#   R_pf_<c>.csv      one row per asset (single portfolios) or per ordered pair i/j,
#                     i != j (pair methods): mean, sd, q5, q95 of the weights, hedge
#                     weights or hedge ratios; HE; the p-values of five tests of equal
#                     variances of the asset and the portfolio returns; the annualized
#                     return, standard deviation and Sharpe ratio of the portfolio
#   R_pf_<c>_ret.csv  portfolio returns per period (thinned), same order as the rows
#   R_pf_<c>_w.csv    weights, hedge weights or hedge ratios per period (thinned)
# Weights and portfolio returns are taken from the R output objects. Hedging
# effectiveness is 1 - var(portfolio)/var(asset), as in the R functions. The p-values
# are computed here for every row from the same two series, because HedgeRatio()
# fills its p-values by row and keeps only the last pair of each row, and because
# VarianceTest()'s "Brown-Forsythe" option tests means: p_levene centres on the
# group means and p_bf on the group medians (car::leveneTest's default, which
# VarianceTest() labels "Levene").

if (!exists("ConnectednessApproach")) {
  lib = Sys.getenv("R_LIBS_USER")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(ConnectednessApproach))
  cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
}
suppressPackageStartupMessages({library(zoo); library(PerformanceAnalytics); library(car)})
cat("PerformanceAnalytics version:", as.character(packageVersion("PerformanceAnalytics")),
    " riskParityPortfolio version:", as.character(packageVersion("riskParityPortfolio")), "\n")

out = file.path("..", "fixtures")
dir.create(out, showWarnings = FALSE)

getx = function(d) {
  x = get(data(list = d))
  if (!is(x, "zoo")) x = zoo::zoo(as.matrix(x[, -1]), as.Date(x[, 1]))
  x
}
q7 = function(v, p) as.numeric(stats::quantile(v, p, type = 7))
wsum = function(w) c(mean(w), sd(w), q7(w, 0.05), q7(w, 0.95))
pvals = function(a, b) {
  y = c(a, b)
  g = factor(rep(c("A", "B"), c(length(a), length(b))))
  c(var.test(a, b)$p.value,
    bartlett.test(y, g)$p.value,
    fligner.test(y, g)$p.value,
    car::leveneTest(y, g, center = mean)[["Pr(>F)"]][1],
    car::leveneTest(y, g, center = median)[["Pr(>F)"]][1])
}
perf = function(pr, idx) {
  z = zoo::zoo(as.numeric(pr), order.by = idx)
  c(as.numeric(PerformanceAnalytics::Return.annualized(z)),
    as.numeric(PerformanceAnalytics::StdDev.annualized(z)),
    as.numeric(PerformanceAnalytics::SharpeRatio(z, FUN = "StdDev", annualize = TRUE)))
}
row1 = function(w, a, pr, idx) c(wsum(w), 1 - var(pr)/var(a), pvals(a, pr), perf(pr, idx))
thin = function(m) {
  nt = nrow(m)
  m = cbind(period = seq_len(nt), m)
  step = if (nt > 2000) 50 else 10
  if (nt > 300) m = m[(m[, "period"] %% step == 1) | (m[, "period"] == nt), , drop = FALSE]
  m
}
cols = c("mean", "sd", "q5", "q95", "he", "p_fisher", "p_bartlett", "p_fk", "p_levene",
         "p_bf", "ret", "risk", "sr")
writepf = function(nm, tab, ret, w) {
  colnames(tab) = cols
  write.csv(tab, file.path(out, paste0("R_pf_", nm, ".csv")), row.names = FALSE)
  write.csv(thin(ret), file.path(out, paste0("R_pf_", nm, "_ret.csv")), row.names = FALSE)
  write.csv(thin(w), file.path(out, paste0("R_pf_", nm, "_w.csv")), row.names = FALSE)
  cat(sprintf("%-22s rows=%3d periods=%5d  mean HE=%.6f  SR[1]=%.6f\n", nm, nrow(tab),
              nrow(ret), mean(tab[, "he"]), tab[1, "sr"]))
}

# single portfolios: MCP, RPP, EWP (one return series)
pf_single = function(nm, fit, x) {
  k = ncol(x); idx = zoo::index(x); X = as.matrix(x)
  W = as.matrix(fit$portfolio_weights)
  pr = as.numeric(fit$portfolio_return)
  tab = t(sapply(1:k, function(i) row1(W[, i], X[, i], pr, idx)))
  writepf(nm, tab, matrix(pr, ncol = 1), W)
}
# multivariate hedging: one hedged series per asset, k - 1 betas each
pf_mhp = function(nm, fit, x) {
  k = ncol(x); idx = zoo::index(x); X = as.matrix(x)
  tab = c(); W = c()
  for (i in 1:k) {
    pr = as.numeric(fit$portfolio_return[, i])
    for (j in (1:k)[-i]) {
      b = as.numeric(fit$Beta[i, j, ])
      tab = rbind(tab, row1(b, X[, i], pr, idx)); W = cbind(W, b)
    }
  }
  writepf(nm, tab, as.matrix(fit$portfolio_return), W)
}
# pair methods: weight of i in the i/j portfolio (BivariatePortfolio stores it as
# portfolio_weights[j, i, ]) or hedge ratio of i on j (hedge_ratio[i, j, ])
pf_pairs = function(nm, fit, x, hedge) {
  k = ncol(x); idx = zoo::index(x); X = as.matrix(x)
  tab = c(); W = c(); R = c()
  for (i in 1:k) for (j in (1:k)[-i]) {
    if (hedge) { w = as.numeric(fit$hedge_ratio[i, j, ]); pr = as.numeric(fit$portfolio_return[i, j, ]) }
    else       { w = as.numeric(fit$portfolio_weights[j, i, ]); pr = as.numeric(fit$portfolio_return[j, i, ]) }
    tab = rbind(tab, row1(w, X[, i], pr, idx)); W = cbind(W, w); R = cbind(R, pr)
  }
  writepf(nm, tab, R, W)
}

# ---------------------------------------------------------------------------
# A. Static VAR(1) of the g2020 exchange rates (returns in percent; the portfolio
#    functions get decimal returns, x/100, as in the R package's examples)
x = getx("g2020")
xr = x/100
fit = VAR(x, configuration = list(nlag = 1))
dca = TimeConnectedness(Phi = fit$B, Sigma = fit$Q, nfore = 10, generalized = TRUE)
pci = dca$PCI[, , 1]
Q = matrix(as.numeric(fit$Q), ncol(x), dimnames = list(colnames(x), colnames(x)))
pf_single("var_mcp_pci", MinimumConnectednessPortfolio(xr, pci), xr)
pf_single("var_mcp_cov_short", MinimumConnectednessPortfolio(xr, Q, long = FALSE), xr)
pf_single("var_mcp_cor", MinimumConnectednessPortfolio(xr, cov2cor(Q)), xr)
pf_single("var_rpp_cov", RiskParityPortfolio(xr, Q), xr)
pf_single("var_ewp", EquallyWeightedPortfolio(xr), xr)
pf_pairs("var_bivariate_cov", BivariatePortfolio(xr, Q), xr, FALSE)
pf_pairs("var_hedge_cov", HedgeRatio(xr, Q, method = "cumprod"), xr, TRUE)
pf_mhp("var_mhp_cov", MultivariateHedgingPortfolio(xr, Q), xr)

# ---------------------------------------------------------------------------
# B. TVP-VAR(1) of the acg2020 monthly exchange rates (kappa1 0.99, kappa2 0.96,
#    Bayes prior from 200 observations, horizon 12): time-varying matrices
x = getx("acg2020")
k = ncol(x)
nlag = 1
prior = BayesPrior(x = x, size = 200 - nlag, nlag = nlag)
fit = TVPVAR(x, configuration = list(l = c(0.99, 0.96), nlag = nlag, prior = prior))
dca = TimeConnectedness(Phi = fit$B_t, Sigma = fit$Q_t, nfore = 12, generalized = TRUE)
# R's TVPVAR() returns one covariance matrix per observation (the first from the
# initial values), so the returns are the whole sample
xr = x/100
Qt = fit$Q_t
dimnames(Qt) = list(colnames(x), colnames(x), as.character(zoo::index(xr)))
PCIt = dca$PCI
dimnames(PCIt) = dimnames(Qt)
cat("acg2020 TVP-VAR periods:", dim(Qt)[3], " returns:", nrow(xr), "\n")
pf_single("tvp_mcp_pci", MinimumConnectednessPortfolio(xr, PCIt), xr)
pf_single("tvp_rpp_cov", RiskParityPortfolio(xr, Qt), xr)
pf_pairs("tvp_bivariate_cov", BivariatePortfolio(xr, Qt, method = "cumprod"), xr, FALSE)
pf_pairs("tvp_hedge_cov", HedgeRatio(xr, Qt), xr, TRUE)
pf_mhp("tvp_mhp_cov", MultivariateHedgingPortfolio(xr, Qt), xr)

# ---------------------------------------------------------------------------
# C. DCC-GARCH(1,1) of the g2020 exchange rates: conditional covariances H_t of
#    rmgarch::dccfit() with the specification of ConnectednessApproach(); the
#    estimates are written to R_pf_dcc_par.csv for dccparams()
suppressPackageStartupMessages({library(rugarch); library(rmgarch)})
x = getx("g2020")
k = ncol(x)
us = rugarch::ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                         variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                         distribution.model = "norm")
ms = rmgarch::dccspec(uspec = rugarch::multispec(replicate(k, us)), dccOrder = c(1, 1),
                      distribution = "mvnorm")
dfit = rmgarch::dccfit(ms, data = x)
write.table(t(dfit@mfit$matcoef[, 1]), file.path(out, "R_pf_dcc_par.csv"), sep = ",",
            row.names = FALSE, col.names = TRUE)
Ht = rmgarch::rcov(dfit)
dimnames(Ht) = list(colnames(x), colnames(x), as.character(zoo::index(x)))
xr = x/100
pf_mhp("dcc_mhp_cov", MultivariateHedgingPortfolio(xr, Ht), xr)
pf_pairs("dcc_bivariate_cov", BivariatePortfolio(xr, Ht), xr, FALSE)
pf_pairs("dcc_hedge_cov", HedgeRatio(xr, Ht), xr, TRUE)
pf_single("dcc_mcp_cor", MinimumConnectednessPortfolio(xr, rmgarch::rcor(dfit)), xr)
cat("rugarch version:", as.character(packageVersion("rugarch")),
    " rmgarch version:", as.character(packageVersion("rmgarch")), "\n")
