# Reference fixtures for version 0.8.2: scalar normalizations of the generalized
# FEVD (Caloia, Cipollini and Muzzioli 2019) and adaptive LASSO, ridge and
# elastic-net VARs (Demirer, Diebold, Liu and Yilmaz 2018).
#
# Run from the test/r folder:   Rscript make_fixtures_082.R
# Requires ConnectednessApproach 1.0.4 and glmnet (installed with it).
#
# Normalizations. ConnectednessApproach has no option for them, so the unnormalized
# generalized FEVD of Pesaran and Shin (1998) is computed here from the VAR() and
# TVPVAR() estimates. Its row-normalized version is first checked against the
# package's TimeConnectedness(); the scalar versions are then passed to the package's
# ConnectednessTable(), which computes TCI, TO, FROM, NET, NPT, NPDC, PCI and
# influence from any decomposition.
#   maxrow:   W / max(rowSums(W))
#   maxcol:   W / max(colSums(W))
#   spectral: W / max(Mod(eigen(W)$values))
#
# Adaptive elastic net. As for make_fixtures_ext.R, ElasticNetVAR() is replaced by
# a copy with interleaved folds (folds(interleaved) in Stata); here it also passes
# penalty.factor = 1/|b_OLS| to cv.glmnet(), b_OLS being the OLS slopes of the
# equation with an intercept.
#
# Output in ../fixtures: R_<case>_ts.csv (per-period measures, thinned) and
# R_<case>_fevd.csv (time-averaged decomposition, percent).

if (!exists("ConnectednessApproach")) {
  lib = Sys.getenv("R_LIBS_USER")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(ConnectednessApproach))
  cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
}
suppressPackageStartupMessages(library(zoo))
cat("glmnet version:", as.character(packageVersion("glmnet")), "\n")

out = file.path("..", "fixtures")
dir.create(out, showWarnings = FALSE)

getx = function(d) {
  x = get(data(list = d))
  if (!is(x, "zoo")) x = zoo::zoo(as.matrix(x[, -1]), as.Date(x[, 1]))
  x
}
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
pairmat = function(A, lower) {
  k = nrow(A)
  out = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) out = c(out, if (lower) A[j, i] else A[i, j])
  out
}
# measures of one period from ConnectednessTable(); our npdc_i_j is R's NPDC[j, i]
permeas = function(W) {
  k = nrow(W)
  ct = ConnectednessTable(array(W, c(k, k, 1)))
  c(ct$TCI, ct$TO, ct$FROM, ct$NET, ct$NPT, pairmat(ct$NPDC, TRUE), pairmat(ct$PCI, FALSE),
    pairmat(ct$INFLUENCE, FALSE))
}
measnames = function(k) c("tci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k),
                          paste0("npt_", 1:k), paste0("npdc_", pairs(k)), paste0("pci_", pairs(k)),
                          paste0("inf_", pairs(k)))
writefx = function(nm, ts, avg) {
  ts = thin(cbind(period = seq_len(nrow(ts)), ts))
  write.csv(as.data.frame(ts), file.path(out, paste0("R_", nm, "_ts.csv")), row.names = FALSE)
  write.table(avg, file.path(out, paste0("R_", nm, "_fevd.csv")), sep = ",",
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("%-26s periods=%5d mean TCI=%.6f\n", nm, nrow(ts), mean(ts[, 2])))
}

# ---------------------------------------------------------------------------
# Unnormalized generalized FEVD: sigma_jj^-1 sum_h (Psi_h Sigma)_ij^2 / sum_h (Psi_h Sigma Psi_h')_ii
ma = function(Phi, H) {
  k = nrow(Phi)
  p = ncol(Phi)/k
  Psi = array(0, c(k, k, H))
  Psi[, , 1] = diag(k)
  if (H > 1) for (h in 2:H) {
    P = matrix(0, k, k)
    for (j in 1:min(h - 1, p)) P = P + Phi[, ((j - 1)*k + 1):(j*k)] %*% Psi[, , h - j]
    Psi[, , h] = P
  }
  Psi
}
gfevd_raw = function(Phi, Sigma, H) {
  Psi = ma(Phi, H)
  k = nrow(Sigma)
  num = matrix(0, k, k)
  den = rep(0, k)
  for (h in 1:H) {
    num = num + sweep((Psi[, , h] %*% Sigma)^2, 2, diag(Sigma), "/")
    den = den + diag(Psi[, , h] %*% Sigma %*% t(Psi[, , h]))
  }
  num/den
}
normalize = function(W, how) {
  switch(how,
         row = W/rowSums(W),
         maxrow = W/max(rowSums(W)),
         maxcol = W/max(colSums(W)),
         spectral = W/max(Mod(eigen(W, only.values = TRUE)$values)))
}
# per-period decompositions (k x k x t) of B_t, Q_t; checks the row-normalized
# version against the package's TimeConnectedness()
normfx = function(nm, Bt, Qt, H, how) {
  k = dim(Qt)[1]
  nt = dim(Qt)[3]
  W = array(NA, c(k, k, nt))
  ts = NULL
  dmax = 0
  for (t in 1:nt) {
    B = if (length(dim(Bt)) == 3) Bt[, , t] else Bt
    Q = matrix(Qt[, , t], k, k)
    raw = gfevd_raw(B, Q, H)
    ref = TimeConnectedness(Phi = array(B, c(dim(B), 1)), Sigma = array(Q, c(k, k, 1)),
                            nfore = H, generalized = TRUE)$CT[, , 1]
    dmax = max(dmax, abs(normalize(raw, "row") - ref))
    W[, , t] = normalize(raw, how)
    ts = rbind(ts, permeas(W[, , t]))
  }
  colnames(ts) = measnames(k)
  cat(sprintf("  row-normalized check against TimeConnectedness(): max abs difference %.2e\n", dmax))
  writefx(nm, ts, apply(100*W, 1:2, mean))
}

# A. Static VAR(4) of the Diebold and Yilmaz (2012) volatilities, horizon 10: the
#    replication of Caloia, Cipollini and Muzzioli (2019, Table 11)
x = getx("dy2012")
v = VAR(x, configuration = list(nlag = 4))
Q = array(as.numeric(v$Q), c(ncol(x), ncol(x), 1))
for (how in c("maxrow", "maxcol", "spectral")) normfx(paste0("dy2012_norm_", how), v$B, Q, 10, how)

# B. TVP-VAR(1) of acg2020 (kappa2 0.96, Bayes prior from 200 observations,
#    horizon 12): one decomposition per period
x = getx("acg2020")
prior = BayesPrior(x = x, size = 199, nlag = 1)
tv = TVPVAR(x, configuration = list(l = c(0.99, 0.96), nlag = 1, prior = prior))
normfx("acg2020_tvp_norm_maxrow", tv$B_t, tv$Q_t, 12, "maxrow")

# ---------------------------------------------------------------------------
# Adaptive LASSO, ridge and elastic-net VARs
ElasticNetVAR_adaptive = function(x, configuration=list(nlag=1, nfolds=10, loss="mae", alpha=NULL, n_alpha=10)) {
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
    z = embed(x, nlag+1)
    X = z[,-c(1:k)]
    y = z[,i]
    w = 1/abs(coef(lm(y ~ X))[-1])
    for (j in 1:length(alpha_)) {
      fit = glmnet::cv.glmnet(X, y, alpha=alpha_[j], type.measure=loss, nfolds=nfolds,
                              foldid=rep(seq(nfolds), length=nrow(X)), penalty.factor=w)
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
pairser = function(A, lower) {
  k = dim(A)[1]
  out = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) out = cbind(out, if (lower) A[j, i, ] else A[i, j, ])
  out
}
time_fx = function(nm, dca) {
  CT = 100 * dca$CT
  k = dim(CT)[1]
  ts = cbind(dca$TCI, dca$TO, dca$FROM, dca$NET, dca$NPT,
             pairser(dca$NPDC, TRUE), pairser(dca$PCI, FALSE), pairser(dca$INFLUENCE, FALSE))
  ts = matrix(ts, nrow = dim(CT)[3])
  colnames(ts) = measnames(k)
  writefx(nm, ts, apply(CT, 1:2, mean))
}
aenrun = function(d, nlag, nfore, model, ws = NULL, rows = NULL, en = list()) {
  x = get(data(list = d))
  if (!is.null(rows)) x = x[rows, ]
  cfg = modifyList(list(nfolds = 10, alpha = NULL, loss = "mae", n_alpha = 10), en)
  ns = if ("ConnectednessApproach" %in% loadedNamespaces()) asNamespace("ConnectednessApproach")
    else globalenv()
  old = get("ElasticNetVAR", envir = ns)
  if (identical(ns, globalenv())) assign("ElasticNetVAR", ElasticNetVAR_adaptive, envir = ns)
  else utils::assignInNamespace("ElasticNetVAR", ElasticNetVAR_adaptive, "ConnectednessApproach")
  on.exit(if (identical(ns, globalenv())) assign("ElasticNetVAR", old, envir = ns)
          else utils::assignInNamespace("ElasticNetVAR", old, "ConnectednessApproach"))
  suppressMessages(ConnectednessApproach(x, nlag = nlag, nfore = nfore, window.size = ws, model = model,
                   connectedness = "Time", VAR_config = list(ElasticNet = cfg)))
}
aencases = list(
  dy2012_alasso_static    = list(d = "dy2012", nlag = 2, nfore = 10, model = "LASSO"),
  dy2012_aridge_static    = list(d = "dy2012", nlag = 2, nfore = 10, model = "Ridge"),
  dy2012_aenet_a05        = list(d = "dy2012", nlag = 2, nfore = 10, model = "Elastic",
                                 en = list(alpha = 0.5)),
  dy2012_aenet4_mse       = list(d = "dy2012", nlag = 1, nfore = 10, model = "Elastic",
                                 en = list(n_alpha = 4, loss = "mse")),
  dy2012_alasso_roll200   = list(d = "dy2012", nlag = 1, nfore = 10, model = "LASSO", ws = 200,
                                 rows = 1:300)
)
for (nm in names(aencases)) {
  cs = aencases[[nm]]
  time_fx(nm, aenrun(cs$d, cs$nlag, cs$nfore, cs$model, cs$ws, cs$rows,
                     if (is.null(cs$en)) list() else cs$en))
}
