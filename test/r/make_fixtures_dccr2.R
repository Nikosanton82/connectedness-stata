# Reference fixtures for version 0.8.1: R2 decomposed connectedness of DCC-GARCH
# conditional correlations (Cocca, Gabauer and Pomberger 2024) with R2Correlations(),
# and the conditional and partial correlations of estat correlation with
# ConditionalCorrelation() and PartialCorrelations().
#
# Run from the test/r folder:   Rscript make_fixtures_dccr2.R
# Requires ConnectednessApproach 1.0.4, rugarch and rmgarch.
#
# Output in ../fixtures:
#   R_dccr2_<c>_par.csv    rmgarch estimates (for dccparams()) and log likelihood
#   R_dccr2_<c>_ts.csv     per-period tci, to_i, from_i, net_i, npdc_i_j, pci_i_j
#   R_dccr2_<c>_fevd.csv   time average of the R2 decomposition (percent, own = 100)
#   R_corr_<c>_ts.csv      per-period cor_i_j and pcor_i_j (i < j)
#   R_corr_<c>_mean.csv    time averages: correlations, then partial correlations
#                          (k rows each, diagonal 1)
# R2Correlations() reports the PCI as 2(a + b)/(a + b + 1 + 1), on a scale of 0 to 2;
# the fixtures multiply it by 100, the percent scale of ConnectednessTable() and of
# the Stata results. PartialCorrelations() puts -1 on the diagonal; the averages here
# use 1, as estat correlation does.

if (!exists("ConnectednessApproach")) {
  lib = Sys.getenv("R_LIBS_USER")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(ConnectednessApproach))
  cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
}
suppressPackageStartupMessages({library(zoo); library(rugarch); library(rmgarch)})
cat("rugarch version:", as.character(packageVersion("rugarch")),
    " rmgarch version:", as.character(packageVersion("rmgarch")), "\n")

out = file.path("..", "fixtures")
dir.create(out, showWarnings = FALSE)

getx = function(d) {
  x = get(data(list = d))
  if (!is(x, "zoo")) x = zoo::zoo(as.matrix(x[, -1]), as.Date(x[, 1]))
  x
}
pairs = function(k) {
  pr = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) pr = c(pr, paste0(i, "_", j))
  pr
}
pairser = function(A, lower) {
  k = dim(A)[1]
  out = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) out = cbind(out, if (lower) A[j, i, ] else A[i, j, ])
  out
}
thin = function(ts) {
  nt = nrow(ts)
  step = if (nt > 2000) 50 else 10
  if (nt > 300) ts = ts[(ts[, "period"] %% step == 1) | (ts[, "period"] == nt), , drop = FALSE]
  ts
}
dccfit = function(x) {
  k = ncol(x)
  us = rugarch::ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                           variance.model = list(garchOrder = c(1, 1), model = "sGARCH"),
                           distribution.model = "norm")
  ms = rmgarch::dccspec(uspec = rugarch::multispec(replicate(k, us)), dccOrder = c(1, 1),
                        distribution = "mvnorm")
  rmgarch::dccfit(ms, data = x)
}

# ---------------------------------------------------------------------------
# R2 decomposed connectedness of the DCC-GARCH conditional correlations
dccr2 = function(nm, d) {
  x = getx(d)
  k = ncol(x)
  fit = dccfit(x)
  par = c(fit@mfit$matcoef[, 1], llh = fit@mfit$llh)
  write.table(t(par), file.path(out, paste0("R_dccr2_", nm, "_par.csv")), sep = ",",
              row.names = FALSE, col.names = TRUE)
  R = rmgarch::rcor(fit)
  dimnames(R) = list(colnames(x), colnames(x), as.character(zoo::index(x)))
  dca = R2Correlations(R)
  nt = dim(dca$CT)[3]
  ts = cbind(dca$TCI, dca$TO, dca$FROM, dca$NET, pairser(dca$NPDC, TRUE),
             100*pairser(dca$PCI, FALSE))
  colnames(ts) = c("tci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k),
                   paste0("npdc_", pairs(k)), paste0("pci_", pairs(k)))
  ts = thin(cbind(period = seq_len(nt), ts))
  write.csv(as.data.frame(ts), file.path(out, paste0("R_dccr2_", nm, "_ts.csv")), row.names = FALSE)
  write.table(apply(100*dca$CT, 1:2, mean), file.path(out, paste0("R_dccr2_", nm, "_fevd.csv")),
              sep = ",", row.names = FALSE, col.names = FALSE)
  cat(sprintf("%-22s periods=%5d mean TCI=%.6f\n", paste0("dccr2_", nm), nt, mean(dca$TCI)))
  fit
}

# ---------------------------------------------------------------------------
# Conditional and partial correlations of a k x k x t covariance array
corrfx = function(nm, Q) {
  k = dim(Q)[1]
  nt = dim(Q)[3]
  C = ConditionalCorrelation(Q)
  P = PartialCorrelations(Q)
  ts = cbind(pairser(C, FALSE), pairser(P, FALSE))
  ts = matrix(ts, nrow = nt)
  colnames(ts) = c(paste0("cor_", pairs(k)), paste0("pcor_", pairs(k)))
  write.csv(as.data.frame(thin(cbind(period = seq_len(nt), ts))),
            file.path(out, paste0("R_corr_", nm, "_ts.csv")), row.names = FALSE)
  Cm = apply(C, 1:2, mean); diag(Cm) = 1
  Pm = apply(P, 1:2, mean); diag(Pm) = 1
  write.table(rbind(Cm, Pm), file.path(out, paste0("R_corr_", nm, "_mean.csv")), sep = ",",
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("%-22s periods=%5d mean cor_1_2=%.6f mean pcor_1_2=%.6f\n", paste0("corr_", nm),
              nt, Cm[1, 2], Pm[1, 2]))
}

fit = dccr2("g2020", "g2020")
H = rmgarch::rcov(fit)
corrfx("dcc_g2020", H)
invisible(dccr2("cgp2024", "cgp2024"))

# static VAR(1) of g2020: one residual covariance matrix
x = getx("g2020")
v = VAR(x, configuration = list(nlag = 1))
Q = array(as.numeric(v$Q), c(ncol(x), ncol(x), 1))
corrfx("var_g2020", Q)

# TVP-VAR(1) of acg2020 (kappa2 0.96, Bayes prior from 200 observations)
x = getx("acg2020")
prior = BayesPrior(x = x, size = 199, nlag = 1)
tv = TVPVAR(x, configuration = list(l = c(0.99, 0.96), nlag = 1, prior = prior))
corrfx("tvp_acg2020", tv$Q_t)
