# Reference fixtures for version 0.9.3: the system-wide within-cluster and cross-cluster
# parts of the TCI (Buchwalter, Diebold and Yilmaz 2026a, Fig. 5; 2026b, Section 2.3)
# and the cluster table with its TO and NET rows (2026b, Tables 1 to 3), for the cases
# of make_fixtures_090.R. The clustered decomposition is the independent implementation
# of make_fixtures_090.R (sequential projection, eqs. 7 to 9 of the published paper,
# averaged over the cluster orderings); the new measures are computed here directly from
# their definitions:
#   C_within = (1/N) sum_i sum_{j != i, same cluster} theta_ij,
#   C_cross  = (1/N) sum_i sum_{j in another cluster} theta_ij,   C_within + C_cross = TCI;
#   T[c, d]  = (1/N_c) sum_{i in c} sum_{j in d} theta_ij (rows sum to 100);
#   TO_c     = sum_{d != c} T[d, c],  FROM_c = sum_{d != c} T[c, d],  NET_c = TO_c - FROM_c.
#
# Run from the test/r folder:   Rscript make_fixtures_093.R
# Output in ../fixtures:
#   R_<case>_clts093.csv   cl_within, cl_cross, cl_tto_c and cl_tnet_c per period (thinned)
#   R_<case>_cltab.csv     cluster table of the averaged decomposition: C rows, TO, NET;
#                          columns the C clusters and FROM ("." in the two empty cells)
#   R_<case>_clsys.csv     within and cross parts of the averaged decomposition


if (!exists("ConnectednessApproach")) {
  lib = Sys.getenv("R_LIBS_USER")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(ConnectednessApproach))
  cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
}
suppressPackageStartupMessages(library(zoo))

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
permeas = function(W) {
  k = nrow(W)
  ct = ConnectednessTable(array(W, c(k, k, 1)))
  c(ct$TCI, ct$TO, ct$FROM, ct$NET, ct$NPT, pairmat(ct$NPDC, TRUE), pairmat(ct$PCI, FALSE),
    pairmat(ct$INFLUENCE, FALSE))
}
measnames = function(k) c("tci", paste0("to_", 1:k), paste0("from_", 1:k), paste0("net_", 1:k),
                          paste0("npt_", 1:k), paste0("npdc_", pairs(k)), paste0("pci_", pairs(k)),
                          paste0("inf_", pairs(k)))

# ---------------------------------------------------------------------------
# moving-average matrices A_0, ..., A_{H-1}
ma = function(Phi, H) {
  k = nrow(Phi)
  p = ncol(Phi)/k
  A = array(0, c(k, k, H))
  A[, , 1] = diag(k)
  if (H > 1) for (h in 2:H) {
    P = matrix(0, k, k)
    for (j in 1:min(h - 1, p)) P = P + Phi[, ((j - 1)*k + 1):(j*k)] %*% A[, , h - j]
    A[, , h] = P
  }
  A
}
# Q_C^-1 for the clusters in the order of the variables (cl: cluster of each variable,
# clusters appear in contiguous blocks): row block b is (-Sigma_{b,<b} Sigma_{<b,<b}^-1, I)
qinv = function(Sigma, cl) {
  k = nrow(Sigma)
  Qi = diag(k)
  blocks = unique(cl)
  for (m in seq_along(blocks)[-1]) {
    b = which(cl == blocks[m])
    prev = which(cl %in% blocks[1:(m - 1)])
    Qi[b, prev] = -Sigma[b, prev, drop = FALSE] %*% solve(Sigma[prev, prev, drop = FALSE])
  }
  Qi
}
clustvd_order = function(A, Sigma, cl, ord) {
  k = nrow(Sigma)
  H = dim(A)[3]
  # variables of the clusters in the order ord, keeping their order within clusters
  idx = unlist(lapply(ord, function(c) which(cl == c)))
  S = Sigma[idx, idx]
  Qi = qinv(S, cl[idx])
  Q = solve(Qi)
  Om = Qi %*% S %*% t(Qi)
  # back to the original order of the variables
  P = diag(k)[idx, ]                      # x_perm = P x
  QO = t(P) %*% (Q %*% Om) %*% P
  om = diag(t(P) %*% Om %*% P)
  num = matrix(0, k, k)
  den = rep(0, k)
  for (h in 1:H) {
    num = num + (A[, , h] %*% QO)^2
    den = den + diag(A[, , h] %*% Sigma %*% t(A[, , h]))
  }
  W = sweep(num, 2, om, "/")/den
  W/rowSums(W)
}
allperms = function(v) {
  if (length(v) <= 1) return(list(v))
  out = list()
  for (i in seq_along(v)) for (p in allperms(v[-i])) out[[length(out) + 1]] = c(v[i], p)
  out
}
clustvd = function(Phi, Sigma, H, cl, average = TRUE) {
  A = ma(Phi, H)
  ords = if (average) allperms(sort(unique(cl))) else list(sort(unique(cl)))
  Reduce(`+`, lapply(ords, function(o) clustvd_order(A, Sigma, cl, o)))/length(ords)
}

# the new measures of one decomposition W (percent)
newmeas = function(W, cl) {
  k = nrow(W)
  C = max(cl)
  same = outer(cl, cl, "==")
  off = !diag(k)
  wth = sum(W[same & off])/k
  crs = sum(W[!same])/k
  T = matrix(NA, C, C)
  for (c in 1:C) for (d in 1:C) T[c, d] = sum(W[cl == c, cl == d, drop = FALSE])/sum(cl == c)
  from = rowSums(T) - diag(T)
  to = colSums(T) - diag(T)
  list(wth = wth, crs = crs, T = T, from = from, to = to, net = to - from)
}
clustfx = function(nm, Bt, Qt, H, cl, average = TRUE) {
  k = dim(Qt)[1]
  nt = dim(Qt)[3]
  C = max(cl)
  W = array(NA, c(k, k, nt))
  cts = NULL
  for (t in 1:nt) {
    B = if (length(dim(Bt)) == 3) Bt[, , t] else Bt
    Q = matrix(Qt[, , t], k, k)
    W[, , t] = 100*clustvd(B, Q, H, cl, average)
    m = newmeas(W[, , t], cl)
    tci = sum(W[, , t][!diag(k)])/k
    stopifnot(abs(m$wth + m$crs - tci) < 1e-10)
    cts = rbind(cts, c(m$wth, m$crs, m$to, m$net))
  }
  colnames(cts) = c("cl_within", "cl_cross", paste0("cl_tto_", 1:C), paste0("cl_tnet_", 1:C))
  avg = apply(W, 1:2, mean)
  m = newmeas(avg, cl)
  stopifnot(max(abs(rowSums(m$T) - 100)) < 1e-10)
  tab = rbind(cbind(m$T, m$from), c(m$to, NA), c(m$net, NA))
  write.csv(as.data.frame(thin(cbind(period = seq_len(nt), cts))),
            file.path(out, paste0("R_", nm, "_clts093.csv")), row.names = FALSE)
  write.table(tab, file.path(out, paste0("R_", nm, "_cltab.csv")), sep = ",", row.names = FALSE,
              col.names = FALSE, na = ".")
  write.table(t(c(m$wth, m$crs)), file.path(out, paste0("R_", nm, "_clsys.csv")), sep = ",",
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("%-26s periods=%5d within/cross (averaged)=%.4f/%.4f net: %s\n", nm, nt, m$wth, m$crs,
              paste(sprintf("%.3f", m$net), collapse = " ")))
}

# rolling-window VAR as ConnectednessApproach(): window.size - nlag rows per window
rollvar = function(x, nlag, window.size) {
  ws = window.size - nlag
  t0 = nrow(x) - ws + 1
  k = ncol(x)
  Bt = array(NA, c(k, k*nlag, t0))
  Qt = array(NA, c(k, k, t0))
  for (i in 1:t0) {
    v = VAR(x[i:(i + ws - 1), ], configuration = list(nlag = nlag))
    Bt[, , i] = v$B
    Qt[, , i] = v$Q[, , 1]
  }
  list(B = Bt, Q = Qt)
}
staticQ = function(v, k) array(as.numeric(v$Q), c(k, k, 1))

# sanity check: one cluster = generalized FEVD of TimeConnectedness()
x = getx("dy2012")
v = VAR(x, configuration = list(nlag = 4))
g1 = clustvd(v$B, matrix(v$Q[, , 1], 4, 4), 10, rep(1, 4))
gr = TimeConnectedness(Phi = array(v$B, c(dim(v$B), 1)), Sigma = v$Q, nfore = 10, generalized = TRUE)$CT[, , 1]
cat(sprintf("one cluster against TimeConnectedness(): max abs difference %.2e\n", max(abs(g1 - gr))))

# A. Diebold and Yilmaz (2012) volatilities, static VAR(4), H = 10: two clusters
clustfx("dy2012_clust2", v$B, staticQ(v, 4), 10, c(1, 1, 2, 2))
# B. Diebold and Yilmaz (2009) weekly equity returns, static VAR(2), H = 10: three
#    regions (Americas, Europe, Asia-Pacific), averaged over the 6 orderings and in
#    the order given
x = getx("dy2009")
reg = c(US = 1, UK = 2, FRA = 2, GER = 2, HKG = 3, JPN = 3, AUS = 3, IDN = 3, KOR = 3, MYS = 3,
        PHL = 3, SGP = 3, TAI = 3, THA = 3, ARG = 1, BRA = 1, CHL = 1, MEX = 1, TUR = 2)
cl = unname(reg[colnames(x)])
v = VAR(x, configuration = list(nlag = 2))
clustfx("dy2009_regions", v$B, staticQ(v, ncol(x)), 10, cl)
clustfx("dy2009_regions_given", v$B, staticQ(v, ncol(x)), 10, cl, average = FALSE)
# C. rolling VAR(2) of dy2012, 200-observation windows, first 400 rows: clusters of
#    sizes 1 and 3
x = getx("dy2012")[1:400, ]
rv = rollvar(x, 2, 200)
clustfx("dy2012_clust_roll200", rv$B, rv$Q, 10, c(1, 2, 2, 2))
# D. TVP-VAR(1) of acg2020 (kappa2 0.96, Bayes prior from 200 observations, H = 12):
#    European currencies and the yen
x = getx("acg2020")
prior = BayesPrior(x = x, size = 199, nlag = 1)
tv = TVPVAR(x, configuration = list(l = c(0.99, 0.96), nlag = 1, prior = prior))
clustfx("acg2020_tvp_clust", tv$B_t, tv$Q_t, 12, c(1, 1, 2, 1))
# E. six clusters of three yields each (cgs2022 swap curves), static VAR(1), H = 10:
#    720 orderings
x = getx("cgs2022")
v = VAR(x, configuration = list(nlag = 1))
clustfx("cgs2022_countries", v$B, staticQ(v, ncol(x)), 10, rep(1:6, each = 3))

