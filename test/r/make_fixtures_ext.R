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
               tvp = NULL, freq = NULL, tau = 0.5, rows = NULL) {
  x = get(data(list = d))
  if (!is.null(rows)) x = x[rows, ]
  vc = list(TVPVAR = list(kappa1 = 0.99, kappa2 = 0.99, prior = "BayesPrior", gamma = 0.01),
            QVAR = list(tau = tau, method = "fn"))
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

# ---------------------------------------------------------------------------
# Quantile VAR (Chatziantoniou, Gabauer and Stenfors 2021) and LAD VAR, from
# version 0.4.0. Rolling cases use the first rows of the data to keep R's run time
# short (each window is estimated independently).
time_fx = function(nm, dca) {
  CT = 100 * dca$CT
  k = dim(CT)[1]
  nt = dim(CT)[3]
  ts = cbind(dca$TCI, dca$TO, dca$FROM, dca$NET, dca$NPT,
             pairser(dca$NPDC, TRUE), pairser(dca$PCI, FALSE), pairser(dca$INFLUENCE, FALSE))
  colnames(ts) = c(if (isTRUE(dca$config$corrected)) "ctci" else "tci", paste0("to_", 1:k),
                   paste0("from_", 1:k), paste0("net_", 1:k), paste0("npt_", 1:k),
                   paste0("npdc_", pairs(k)), paste0("pci_", pairs(k)), paste0("inf_", pairs(k)))
  writefx(nm, ts, apply(CT, 1:2, mean))
}
qcases = list(
  cgs2021_qvar_roll200  = list(d = "cgs2021", nlag = 1, nfore = 20, ws = 200, corr = TRUE, rows = 1:700),
  cgs2021_qvar_q05      = list(d = "cgs2021", nlag = 1, nfore = 20, corr = TRUE, tau = 0.05),
  cgs2021_qvar_qoq      = list(d = "cgs2021", nlag = 2, nfore = 20, tau = c(0.2, 0.8, 0.5, 0.95)),
  dy2012_lad_static     = list(d = "dy2012", nlag = 4, nfore = 10, model = "LAD"),
  dy2012_lad_roll200    = list(d = "dy2012", nlag = 4, nfore = 10, ws = 200, model = "LAD", rows = 1:600)
)
for (nm in names(qcases)) {
  cs = qcases[[nm]]
  tau = if (is.null(cs$tau)) 0.5 else cs$tau
  time_fx(nm, run(cs$d, cs$nlag, cs$nfore, cs$ws, if (is.null(cs$model)) "QVAR" else cs$model,
                  "Time", isTRUE(cs$corr), tau = tau, rows = cs$rows))
}
gcat = run("gcat2022", 1, 100, 200, "QVAR", "Frequency", freq = list(partition = c(pi + 0.00001, pi/5, 0)),
           rows = 1:600)
freq_fx("gcat2022_qvar_freq_roll200", gcat)
ejoint_fx("jcggh2022_qvar_extjoint_q90", run("jcggh2022", 1, 10, NULL, "QVAR", "Extended Joint", tau = 0.9))

# ---------------------------------------------------------------------------
# R2 decomposed connectedness (Naeem et al. 2023; Balli et al. 2023) and model-free
# (unconditional) connectedness (Gabauer, Chatziantoniou and Stenfors 2023), from
# version 0.5.0. Base names hold the overall measures; with lags, the suffixes _c
# and _l hold the contemporaneous and lagged parts.
r2run = function(d, nlag, ws = NULL, method = "pearson", decomposition = TRUE, corr = FALSE) {
  x = get(data(list = d))
  suppressMessages(ConnectednessApproach(x, nlag = nlag, window.size = ws, connectedness = "R2",
    corrected = corr, Connectedness_config = list(R2Connectedness = list(method = method,
    decomposition = decomposition, relative = FALSE))))
}
r2_fx = function(nm, dca, nlag) {
  CT = 100 * dca$CT                        # k x k x t x (nlag+1 > 1 ? 2 : 1)
  k = dim(CT)[1]
  nt = dim(CT)[3]
  TCI = matrix(dca$TCI, nrow = nt)
  nd = ncol(TCI)
  sufs = if (nd == 3) c("", "_c", "_l") else ""
  ts = NULL
  for (b in 1:nd) {
    tk = function(A) if (nd == 1) matrix(A, nrow = nt) else matrix(A[, , b], nrow = nt)
    kk = function(A) if (nd == 1) array(A, c(k, k, nt)) else array(A[, , , b], c(k, k, nt))
    blk = cbind(TCI[, b], tk(dca$TO), tk(dca$FROM), tk(dca$NET), pairser(kk(dca$NPDC), TRUE))
    colnames(blk) = paste0(c(if (isTRUE(dca$config$corrected)) "ctci" else "tci",
                             paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k),
                             paste0("npdc_", pairs(k))), sufs[b])
    ts = cbind(ts, blk)
  }
  if (nd == 3) {
    avg = apply(CT[, , , 1, drop = FALSE] + CT[, , , 2, drop = FALSE], 1:2, mean)
    extra = list(fevd_c = apply(CT[, , , 1, drop = FALSE], 1:2, mean),
                 fevd_l = apply(CT[, , , 2, drop = FALSE], 1:2, mean))
  } else {
    avg = apply(CT[, , , 1, drop = FALSE], 1:2, mean)
    extra = list()
  }
  writefx(nm, ts, avg, extra)
}
r2cases = list(
  dy2012_r2_roll200          = list(d = "dy2012", nlag = 1, ws = 200),
  dy2012_r2_spearman_roll200 = list(d = "dy2012", nlag = 1, ws = 200, method = "spearman"),
  dy2012_r2_kendall_roll200  = list(d = "dy2012", nlag = 1, ws = 200, method = "kendall"),
  dy2012_r2_lag0_roll200     = list(d = "dy2012", nlag = 0, ws = 200),
  dy2012_r2_lag2_static_corr = list(d = "dy2012", nlag = 2, corr = TRUE),
  cgs2022_r2_static          = list(d = "cgs2022", nlag = 1)
)
for (nm in names(r2cases)) {
  cs = r2cases[[nm]]
  r2_fx(nm, r2run(cs$d, cs$nlag, cs$ws, if (is.null(cs$method)) "pearson" else cs$method,
                  corr = isTRUE(cs$corr)), cs$nlag)
}
mfcases = list(
  dy2012_mf_static          = list(d = "dy2012"),
  dy2012_mf_roll200         = list(d = "dy2012", ws = 200),
  dy2012_mf_kendall_roll200 = list(d = "dy2012", ws = 200, method = "kendall"),
  dy2012_mf_spearman_static = list(d = "dy2012", method = "spearman")
)
for (nm in names(mfcases)) {
  cs = mfcases[[nm]]
  time_fx(nm, r2run(cs$d, 1, cs$ws, if (is.null(cs$method)) "pearson" else cs$method,
                    decomposition = FALSE))
}

# ---------------------------------------------------------------------------
# LASSO, ridge and elastic-net VARs (Gabauer, Gupta, Marfatia and Miller 2024), from
# version 0.6.0. cv.glmnet() draws its folds at random, so R's ElasticNetVAR() is
# replaced, for these fixtures only, by an identical function that passes
# foldid = rep(seq(nfolds), length = n) to cv.glmnet(): the folds that Stata uses with
# folds(interleaved).
ElasticNetVAR_interleaved = function(x, configuration=list(nlag=1, nfolds=10, loss="mae", alpha=NULL, n_alpha=10)) {
  nlag = configuration$nlag
  alpha_ = configuration$alpha
  nfolds = configuration$nfolds
  n_alpha = configuration$n_alpha
  loss = configuration$loss
  NAMES = colnames(x)
  k = ncol(x)
  if (is.null(alpha_)) {
    alpha_ = seq(0, 1, 1/n_alpha)
  }
  alpha = Res = B = NULL
  for (i in 1:k) {
    MAE = NULL
    B_ = array(NA, c(nlag*k, length(alpha_)))
    for (j in 1:length(alpha_)) {
      z = embed(x, nlag+1)
      X = z[,-c(1:k)]
      y = z[,i]
      fit = glmnet::cv.glmnet(X, y, alpha=alpha_[j], type.measure=loss, nfolds=nfolds,
                              foldid=rep(seq(nfolds), length=nrow(X)))
      B_[,j] = coef(fit)[-1]
      y_pred = cbind(1, X) %*% coef(fit)
      MAE[j] = mean(abs(y - y_pred[,1]))
      if (MAE[j] <= min(MAE)) {
        Y_pred = y_pred[,1]
      }
    }
    ind = which(MAE==min(MAE))[1]
    B = rbind(B, B_[,ind])
    Res = cbind(Res, y-Y_pred)
    alpha[i] = alpha_[ind]
  }
  Q = array(t(Res)%*%Res/nrow(Res), c(k, k, 1), dimnames=list(NAMES, NAMES, tail(as.character(zoo::index(x)),1)))
  results = list(B=B, Q=Q, alpha=alpha)
}
enrun = function(d, nlag, nfore, model, ws = NULL, rows = NULL, en = list(), conn = "Time") {
  x = get(data(list = d))
  if (!is.null(rows)) x = x[rows, ]
  cfg = modifyList(list(nfolds = 10, alpha = NULL, loss = "mae", n_alpha = 10), en)
  ns = if ("ConnectednessApproach" %in% loadedNamespaces()) asNamespace("ConnectednessApproach")
    else globalenv()
  old = get("ElasticNetVAR", envir = ns)
  if (identical(ns, globalenv())) assign("ElasticNetVAR", ElasticNetVAR_interleaved, envir = ns)
  else utils::assignInNamespace("ElasticNetVAR", ElasticNetVAR_interleaved, "ConnectednessApproach")
  on.exit(if (identical(ns, globalenv())) assign("ElasticNetVAR", old, envir = ns)
          else utils::assignInNamespace("ElasticNetVAR", old, "ConnectednessApproach"))
  suppressMessages(ConnectednessApproach(x, nlag = nlag, nfore = nfore, window.size = ws, model = model,
                   connectedness = conn, VAR_config = list(ElasticNet = cfg)))
}
encases = list(
  dy2012_lasso_static     = list(d = "dy2012", nlag = 2, nfore = 10, model = "LASSO"),
  dy2012_ridge_static     = list(d = "dy2012", nlag = 2, nfore = 10, model = "Ridge"),
  dy2012_elastic_static   = list(d = "dy2012", nlag = 2, nfore = 10, model = "Elastic"),
  dy2012_elastic_a05_mse  = list(d = "dy2012", nlag = 1, nfore = 10, model = "Elastic",
                                 en = list(alpha = 0.5, loss = "mse")),
  dy2012_lasso_roll200    = list(d = "dy2012", nlag = 1, nfore = 10, model = "LASSO", ws = 200, rows = 1:400),
  dy2012_elastic4_roll200 = list(d = "dy2012", nlag = 1, nfore = 10, model = "Elastic", ws = 200,
                                 rows = 1:300, en = list(n_alpha = 4)),
  gghm2022_lasso_static   = list(d = "gghm2022", nlag = 1, nfore = 10, model = "LASSO")
)
for (nm in names(encases)) {
  cs = encases[[nm]]
  time_fx(nm, enrun(cs$d, cs$nlag, cs$nfore, cs$model, cs$ws, cs$rows,
                    if (is.null(cs$en)) list() else cs$en))
}
