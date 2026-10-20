# Reference fixtures for the connectedness extensions (version 0.3.0 onwards),
# computed with the ConnectednessApproach R package itself.
#
# Run from the test/r folder:   Rscript make_fixtures_ext.R
# Requires the CRAN package ConnectednessApproach 1.0.4 (see make_fixtures.R).
# Output: ../fixtures/R_<case>_ts.csv (per-period measures, thinned) and
#         ../fixtures/R_<case>_fevd.csv (time-averaged connectedness table body),
#         plus R_<case>_fevd_f<b>.csv for each frequency band.
# Unlike make_fixtures.R, the series are taken from the R output objects (TCI, TO,
# FROM, NET, NPT, NPDC, PCI, INFLUENCE) rather than recomputed from the FEVD.

if (!exists("ConnectednessApproach")) {
  lib = Sys.getenv("R_LIBS_USER")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(ConnectednessApproach))
  cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
}

out = file.path("..", "fixtures")
dir.create(out, showWarnings = FALSE)

thin = function(ts) {
  nt = nrow(ts)
  step = if (nt > 2000) 50 else 10
  if (nt > 300) ts = ts[(ts[, "period"] %% step == 1) | (ts[, "period"] == nt), , drop = FALSE]
  ts
}
pairs = function(k) {
  pr = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) pr = c(pr, paste0(i, "_", j))
  pr
}
# lower-triangle pair series: our npdc_i_j (i < j) is R's NPDC[j, i] (net from i to j)
pairser = function(A, lower) {
  k = dim(A)[1]
  out = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) out = cbind(out, if (lower) A[j, i, ] else A[i, j, ])
  matrix(out, nrow = dim(A)[3])
  out
}
writefx = function(nm, ts, avg, extra = list()) {
  ts = thin(cbind(period = seq_len(nrow(ts)), ts))
  write.csv(as.data.frame(ts), file.path(out, paste0("R_", nm, "_ts.csv")), row.names = FALSE)
  write.table(avg, file.path(out, paste0("R_", nm, "_fevd.csv")), sep = ",",
              row.names = FALSE, col.names = FALSE)
  for (e in names(extra)) {
    write.table(extra[[e]], file.path(out, paste0("R_", nm, "_", e, ".csv")), sep = ",",
                row.names = FALSE, col.names = FALSE)
  }
  cat(sprintf("%-26s periods=%5d mean TCI=%.6f\n", nm, nrow(ts), mean(ts[, 2])))
}

run = function(d, nlag, nfore, ws = NULL, model = "VAR", conn = "Time", corr = FALSE,
               tvp = NULL, freq = NULL) {
  x = get(data(list = d))
  vc = list(TVPVAR = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior", gamma = 0.01))
  if (!is.null(tvp)) vc$TVPVAR = modifyList(vc$TVPVAR, tvp)
  cc = list(TimeConnectedness = list(generalized = TRUE),
            FrequencyConnectedness = list(partition = c(pi, pi/2, 0), generalized = TRUE, scenario = "ABS"))
  if (!is.null(freq)) cc$FrequencyConnectedness = modifyList(cc$FrequencyConnectedness, freq)
  suppressMessages(ConnectednessApproach(x, nlag = nlag, nfore = nfore, window.size = ws, model = model,
                   connectedness = conn, corrected = corr, VAR_config = vc, Connectedness_config = cc))
}

# ---------------------------------------------------------------------------
# Frequency connectedness (Barunik and Krehlik 2018; Chatziantoniou, Gabauer and
# Gupta 2023). Band 0 (base names) is the total over the bands; band b gets the
# suffix _f<b>.
freq_fx = function(nm, dca) {
  CT = 100 * dca$CT                        # k x k x t x (1 + bands)
  k = dim(CT)[1]
  nb = dim(CT)[4]
  ts = NULL
  extra = list()
  for (b in 1:nb) {
    suf = if (b == 1) "" else paste0("_f", b - 1)
    nt = dim(CT)[3]
    tk = function(A) matrix(A[, , b], nrow = nt)
    kk = function(A) array(A[, , , b], c(k, k, nt))
    blk = cbind(dca$TCI[, b], tk(dca$TO), tk(dca$FROM), tk(dca$NET), tk(dca$NPT),
                pairser(kk(dca$NPDC), TRUE), pairser(kk(dca$PCI), FALSE),
                pairser(kk(dca$INFLUENCE), FALSE))
    colnames(blk) = paste0(c("tci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k),
                             paste0("npt_", 1:k), paste0("npdc_", pairs(k)), paste0("pci_", pairs(k)),
                             paste0("inf_", pairs(k))), suf)
    ts = cbind(ts, blk)
    if (b > 1) extra[[paste0("fevd_f", b - 1)]] = apply(CT[, , , b, drop = FALSE], 1:2, mean)
  }
  # with corrected, R reports the corrected TCI in TCI; store it under ctci
  if (isTRUE(dca$config$corrected)) colnames(ts) = sub("^tci", "ctci", colnames(ts))
  writefx(nm, ts, apply(CT[, , , 1, drop = FALSE], 1:2, mean), extra)
}

fcases = list(
  dy2012_freq_static  = list(d = "dy2012", nlag = 4, nfore = 100,
                             freq = list(partition = c(pi + 0.00001, pi/5, pi/20, 0))),
  dy2012_freq_default = list(d = "dy2012", nlag = 4, nfore = 10),
  dy2012_freq_roll200 = list(d = "dy2012", nlag = 4, nfore = 50, ws = 200,
                             freq = list(partition = c(pi + 0.00001, pi/5, 0))),
  dy2012_freq_wth     = list(d = "dy2012", nlag = 4, nfore = 100,
                             freq = list(partition = c(pi + 0.00001, pi/5, pi/20, 0), scenario = "WTH")),
  dy2009_freq_orth    = list(d = "dy2009", nlag = 2, nfore = 100,
                             freq = list(partition = c(pi + 0.00001, pi/5, 0), generalized = FALSE)),
  cgg2022_freq_tvp    = list(d = "cgg2022", nlag = 1, nfore = 100, ws = 200, model = "TVP-VAR",
                             freq = list(partition = c(pi + 0.00001, pi/5, 0))),
  acg2020_freq_tvp_corr = list(d = "acg2020", nlag = 1, nfore = 12, ws = 200, model = "TVP-VAR", corr = TRUE,
                             tvp = list(kappa2 = 0.96), freq = list(partition = c(pi + 0.00001, pi/3, pi/12, 0)))
)
for (nm in names(fcases)) {
  cs = fcases[[nm]]
  dca = run(cs$d, cs$nlag, cs$nfore, cs$ws, if (is.null(cs$model)) "VAR" else cs$model, "Frequency",
            isTRUE(cs$corr), cs$tvp, cs$freq)
  freq_fx(nm, dca)
}

# ---------------------------------------------------------------------------
# Joint connectedness (Lastrapes and Wiesen 2021)
joint_fx = function(nm, dca) {
  CT = 100 * dca$CT
  k = dim(CT)[1]
  ts = cbind(dca$TCI[, 1], dca$TO, dca$FROM, dca$NET)
  colnames(ts) = c("tci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k))
  writefx(nm, ts, apply(CT, 1:2, mean))
}
jcases = list(
  lw2021_joint_static  = list(d = "lw2021", nlag = 2, nfore = 30),
  lw2021_joint_roll200 = list(d = "lw2021", nlag = 2, nfore = 30, ws = 200),
  bgu2021_joint_tvp    = list(d = "bgu2021", nlag = 1, nfore = 20, model = "TVP-VAR",
                              tvp = list(prior = "MinnesotaPrior", gamma = 0.1))
)
for (nm in names(jcases)) {
  cs = jcases[[nm]]
  joint_fx(nm, run(cs$d, cs$nlag, cs$nfore, cs$ws, if (is.null(cs$model)) "VAR" else cs$model,
                   "Joint", FALSE, cs$tvp))
}

# ---------------------------------------------------------------------------
# Extended joint connectedness (Balcilar, Gabauer and Umar 2021)
ejoint_fx = function(nm, dca) {
  CT = 100 * dca$CT
  k = dim(CT)[1]
  ts = cbind(dca$TCI[, 1], dca$TO, dca$FROM, dca$NET, dca$NPT,
             pairser(dca$NPDC, TRUE), pairser(dca$PCI, FALSE), pairser(dca$INFLUENCE, FALSE))
  colnames(ts) = c("tci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k),
                   paste0("npt_", 1:k), paste0("npdc_", pairs(k)), paste0("pci_", pairs(k)),
                   paste0("inf_", pairs(k)))
  writefx(nm, ts, apply(CT, 1:2, mean))
}
ecases = list(
  bgu2021_extjoint_tvp = list(d = "bgu2021", nlag = 1, nfore = 20, model = "TVP-VAR",
                              tvp = list(prior = "MinnesotaPrior", gamma = 0.1)),
  lw2021_extjoint_static = list(d = "lw2021", nlag = 2, nfore = 30),
  dy2012_extjoint_roll200 = list(d = "dy2012", nlag = 4, nfore = 10, ws = 200)
)
for (nm in names(ecases)) {
  cs = ecases[[nm]]
  ejoint_fx(nm, run(cs$d, cs$nlag, cs$nfore, cs$ws, if (is.null(cs$model)) "VAR" else cs$model,
                    "Extended Joint", FALSE, cs$tvp))
}
