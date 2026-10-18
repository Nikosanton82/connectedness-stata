# Generate reference fixtures with the ConnectednessApproach R package itself.
#
# Run from the test/r folder:   Rscript make_fixtures.R
# Requires: install.packages("ConnectednessApproach")   (CRAN version 1.0.4)
# Output:   ../fixtures/R_<case>_fevd.csv and ../fixtures/R_<case>_ts.csv, in exactly
#           the layout written by ../python/make_fixtures.py, so the Stata
#           certification script can test against R directly.

lib = Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
}
if (!requireNamespace("ConnectednessApproach", quietly = TRUE)) {
  cat("Installing ConnectednessApproach from CRAN...\n")
  install.packages("ConnectednessApproach", repos = "https://cloud.r-project.org", lib = .libPaths()[1])
}
suppressPackageStartupMessages(library(ConnectednessApproach))
cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")

out = file.path("..", "fixtures")
dir.create(out, showWarnings = FALSE)

cases = list(
  dy2009_orth_static  = list(d = "dy2009", nlag = 2, nfore = 10, ws = NULL, model = "VAR", gen = FALSE),
  dy2009_orth_roll200 = list(d = "dy2009", nlag = 2, nfore = 10, ws = 200, model = "VAR", gen = FALSE),
  dy2012_gen_static   = list(d = "dy2012", nlag = 4, nfore = 10, ws = NULL, model = "VAR"),
  dy2012_gen_roll200  = list(d = "dy2012", nlag = 4, nfore = 10, ws = 200, model = "VAR"),
  acg2020_tvp_bayes   = list(d = "acg2020", nlag = 1, nfore = 12, ws = 200, model = "TVP-VAR",
                             tvp = list(kappa1 = 0.99, kappa2 = 0.96, prior = "BayesPrior")),
  acg2020_tvp_minn    = list(d = "acg2020", nlag = 1, nfore = 12, ws = NULL, model = "TVP-VAR",
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "MinnesotaPrior", gamma = 0.1)),
  acg2020_tvp_unin    = list(d = "acg2020", nlag = 1, nfore = 12, ws = NULL, model = "TVP-VAR",
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "UninformativePrior")),
  acg2020_tvp_lag2    = list(d = "acg2020", nlag = 2, nfore = 12, ws = 200, model = "TVP-VAR",
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior")),
  acg2020_var_roll100 = list(d = "acg2020", nlag = 1, nfore = 12, ws = 100, model = "VAR"),
  gg2018_tvp_bayes    = list(d = "gg2018", nlag = 1, nfore = 10, ws = 200, model = "TVP-VAR",
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior")),
  cg2021_tvp_corr     = list(d = "cg2021", nlag = 1, nfore = 10, ws = 200, model = "TVP-VAR", corr = TRUE,
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior")),
  cgs2022_tvp_corr    = list(d = "cgs2022", nlag = 1, nfore = 10, ws = 250, model = "TVP-VAR", corr = TRUE,
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior")),
  cegg2022_tvp_minn   = list(d = "cegg2022", nlag = 1, nfore = 20, ws = NULL, model = "TVP-VAR", corr = TRUE,
                             tvp = list(kappa1 = 0.99, kappa2 = 0.99, prior = "MinnesotaPrior", gamma = 0.1))
)

per_period = function(ct) {
  k = nrow(ct)
  TO = colSums(ct) - diag(ct)
  FROM = rowSums(ct) - diag(ct)
  NET = TO - FROM
  TCI = mean(FROM)
  NPT = rowSums((ct - t(ct)) < 0)
  npdc = pci = inf = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    npdc = c(npdc, ct[j, i] - ct[i, j])
    pci  = c(pci, 200 * (ct[i, j] + ct[j, i]) / (ct[i, i] + ct[i, j] + ct[j, i] + ct[j, j]))
    inf  = c(inf, 100 * abs(ct[j, i] - ct[i, j]) / (ct[i, j] + ct[j, i]))
  }
  c(TCI, TCI * k / (k - 1), TO, FROM, NET, NPT, npdc, pci, inf)
}

colnames_ts = function(k) {
  pr = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) pr = c(pr, paste0(i, "_", j))
  c("tci", "ctci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k), paste0("npt_", 1:k),
    paste0("npdc_", pr), paste0("pci_", pr), paste0("inf_", pr))
}

thin = function(ts) {
  nt = nrow(ts)
  step = if (nt > 2000) 50 else 10
  if (nt > 300) ts = ts[(ts[, "period"] %% step == 1) | (ts[, "period"] == nt), , drop = FALSE]
  ts
}

DCA = list()
for (nm in names(cases)) {
  cs = cases[[nm]]
  x = get(data(list = cs$d))
  vc = list(TVPVAR = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior", gamma = 0.01))
  if (!is.null(cs$tvp)) vc$TVPVAR = modifyList(vc$TVPVAR, cs$tvp)
  gen = if (is.null(cs$gen)) TRUE else cs$gen
  dca = suppressMessages(ConnectednessApproach(x, nlag = cs$nlag, nfore = cs$nfore, window.size = cs$ws,
          model = cs$model, connectedness = "Time", corrected = isTRUE(cs$corr), VAR_config = vc,
          Connectedness_config = list(TimeConnectedness = list(generalized = gen))))
  CT = 100 * dca$CT
  k = dim(CT)[1]
  nt = dim(CT)[3]
  ts = t(sapply(1:nt, function(i) per_period(CT[, , i])))
  if (nt == 1) ts = matrix(ts, nrow = 1)
  colnames(ts) = colnames_ts(k)
  ts = thin(cbind(period = 1:nt, ts))
  DCA[[nm]] = dca
  avg = apply(CT, 1:2, mean)
  write.table(avg, file.path(out, paste0("R_", nm, "_fevd.csv")),
              sep = ",", row.names = FALSE, col.names = FALSE)
  write.csv(as.data.frame(ts), file.path(out, paste0("R_", nm, "_ts.csv")), row.names = FALSE)
  cat(sprintf("%-22s k=%2d periods=%5d TCI=%.4f\n", nm, k, nt, mean(rowSums(avg) - diag(avg))))
}

# ---------------------------------------------------------------------------
# Group measures computed by the R functions themselves (series only: the R
# tables are rounded text, and R's group PCI and influence are not finite for
# time-domain results)
group_ts = function(res, kind, k, corrected) {
  TO = res$TO; FROM = res$FROM; NET = res$NET; NPT = res$NPT; NPDC = res$NPDC
  n = ncol(TO); nt = nrow(TO)
  tci = rowMeans(FROM)
  out = cbind(tci, tci * n / (n - 1), TO, FROM, NET, NPT)
  pr = c()
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    out = cbind(out, NPDC[j, i, ])
    pr = c(pr, paste0(i, "_", j))
  }
  cn = c("tci", "ctci", paste0("to_", 1:n), paste0("from_", 1:n), paste0("net_", 1:n),
         paste0("npt_", 1:n), paste0("npdc_", pr))
  if (kind == "aggregated") {
    e = as.numeric(res$TCI_ext)
    if (corrected) out = cbind(out, e * (k - 1) / k, e) else out = cbind(out, e, e * k / (k - 1))
    cn = c(cn, "tciext", "ctciext")
  }
  if (kind %in% c("internal", "external")) {
    g = matrix(res$gTCI, nrow = nt)
    m = ncol(g)
    if (corrected) out = cbind(out, g * (k - 1) / k, g) else out = cbind(out, g, g * k / (k - 1))
    cn = c(cn, paste0("gtci_", 1:m), paste0("cgtci_", 1:m))
  }
  colnames(out) = cn
  thin(cbind(period = 1:nt, out))
}

gcases = list(
  gg2018_internal     = list(base = "gg2018_tvp_bayes", kind = "internal", g = list(1:4, 5:8)),
  gg2018_external     = list(base = "gg2018_tvp_bayes", kind = "external", g = list(1:4, 5:8)),
  cgs2022_conditional = list(base = "cgs2022_tvp_corr", kind = "conditional", g = seq(1, 18, 3)),
  cgs2022_aggregated  = list(base = "cgs2022_tvp_corr", kind = "aggregated",
                             g = list(1:3, 4:6, 7:9, 10:12, 13:15, 16:18)),
  cegg2022_internal   = list(base = "cegg2022_tvp_minn", kind = "internal", g = list(1:3, 4:7, 8:11)),
  cegg2022_external   = list(base = "cegg2022_tvp_minn", kind = "external", g = list(1:3, 4:7, 8:11)),
  cegg2022_inclusive  = list(base = "cegg2022_tvp_minn", kind = "inclusive", g = 1:3),
  cegg2022_exclusive  = list(base = "cegg2022_tvp_minn", kind = "exclusive", g = 1:3)
)
for (nm in names(gcases)) {
  gc = gcases[[nm]]
  dca = DCA[[gc$base]]
  k = ncol(dca$NET)
  res = switch(gc$kind,
    aggregated  = AggregatedConnectedness(dca, groups = gc$g),
    conditional = ConditionalConnectedness(dca, group = gc$g),
    inclusive   = InclusiveConnectedness(dca, group = gc$g),
    exclusive   = ExclusiveConnectedness(dca, group = gc$g),
    internal    = InternalConnectedness(dca, groups = gc$g),
    external    = ExternalConnectedness(dca, groups = gc$g))
  ts = group_ts(res, gc$kind, k, isTRUE(dca$config$corrected))
  write.csv(as.data.frame(ts), file.path(out, paste0("R_", nm, "_ts.csv")), row.names = FALSE)
  cat(sprintf("%-22s n=%2d mean TCI=%.4f\n", nm, ncol(res$TO), mean(rowMeans(res$FROM))))
}

# ---------------------------------------------------------------------------
# Asymmetric connectedness (Adekoya et al. 2022 replication page)
data("aaacgo2022")
Y = Yp = Yn = aaacgo2022[-1, ]
for (i in 1:ncol(Y)) {
  z = embed(as.numeric(aaacgo2022[, i]), 2)
  Y[, i] = Yp[, i] = Yn[, i] = 100 * (z[, 1] - z[, 2]) / z[, 2]
  Yp[which(Y[, i] < 0), i] = 0
  Yn[which(Y[, i] > 0), i] = 0
}
for (part in c("all", "pos", "neg")) {
  yy = switch(part, all = Y, pos = Yp, neg = Yn)
  dca = suppressMessages(ConnectednessApproach(yy, model = "TVP-VAR", connectedness = "Time",
          nlag = 1, nfore = 10, window.size = 200,
          VAR_config = list(TVPVAR = list(kappa1 = 0.99, kappa2 = 0.99, prior = "MinnesotaPrior", gamma = 0.1))))
  CT = 100 * dca$CT
  nt = dim(CT)[3]
  ts = t(sapply(1:nt, function(i) per_period(CT[, , i])))
  colnames(ts) = colnames_ts(dim(CT)[1])
  ts = thin(cbind(period = 1:nt, ts))
  avg = apply(CT, 1:2, mean)
  nm = paste0("aaacgo2022_asym_", part)
  write.table(avg, file.path(out, paste0("R_", nm, "_fevd.csv")), sep = ",", row.names = FALSE, col.names = FALSE)
  write.csv(as.data.frame(ts), file.path(out, paste0("R_", nm, "_ts.csv")), row.names = FALSE)
  cat(sprintf("%-22s k=%2d periods=%5d TCI=%.4f\n", nm, dim(CT)[1], nt, mean(rowSums(avg) - diag(avg))))
}
