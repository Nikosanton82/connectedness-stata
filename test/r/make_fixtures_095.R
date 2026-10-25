# Reference fixtures for version 0.9.5: a control cluster fixed first and the
# ordering-sensitivity diagnostic (Buchwalter, Diebold and Yilmaz 2026b, Sections 2.5
# and 2.6, arXiv 2609.05792), and simultaneous clustered orthogonalization
# (Buchwalter, Diebold and Yilmaz 2026c, arXiv 2609.12280). Implemented here from the
# papers' equations, independently of the Stata code:
#   - control cluster first (2026b, eqs. 5 to 7): the bank innovations are residualized
#     on the control cluster (Schur complement Sigma_BB.1 = Sigma_BB - Sigma_B1
#     Sigma_11^-1 Sigma_1B), the BDY block orthogonalization is applied to Sigma_BB.1
#     for each ordering of the other clusters, and Q = [I 0; Sigma_B1 Sigma_11^-1 Q_B],
#     Omega = diag(Sigma_11, Omega_B);
#   - net connectedness of each cluster in the cluster table of each ordering (2026b,
#     eq. 2), its distribution over the orderings, its mean by position, and the paired
#     first-versus-last difference holding the relative order of the other clusters
#     fixed (eq. 4), paired here by text keys of the other clusters' order;
#   - simultaneous clustered orthogonalization (2026c, eq. 24): Q^-1 = D_C R1^-1/2
#     (R1^1/2 R_C R1^1/2)^1/2 R1^-1/2 D1^-1, with the decomposition computed from Q and
#     Omega = Q^-1 Sigma Q^-1' for an arbitrary D_C (the result does not depend on it).
# The model estimates and the standard measures come from ConnectednessApproach (VAR(),
# TVPVAR(), ConnectednessTable()); TimeConnectedness() checks the generalized corner.
#
# Run from the test/r folder:   Rscript make_fixtures_095.R
# Output in ../fixtures, for each case (the files of make_fixtures_090.R and
# make_fixtures_093.R):
#   R_<case>_ts.csv, _fevd.csv, _clts.csv, _cl.csv, _ctg.csv, _cltab.csv, _clsys.csv,
#   _clts093.csv
# and for the cases averaged over orderings:
#   R_<case>_clord.csv    static: one row per ordering, pos_1..pos_C, net_1..net_C
#   R_<case>_cldiag.csv   static: estat clorder tables: C rows (averaged net, mean, sd,
#                         min, max), then the permuted clusters by position, then the
#                         paired first-versus-last (mean, sd, min, max, pairs)
#   R_<case>_clts095.csv  rolling: period, cl_dfl_c, cl_nsd_c

lib = Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ConnectednessApproach))
cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
suppressPackageStartupMessages(library(zoo))

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
wr = function(m, f) write.table(m, file.path(out, f), sep = ",", row.names = FALSE, col.names = FALSE, na = ".")

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
# decomposition from Q (u = Q eps) and Omega = V(eps), normalized by row
vdq = function(A, Sigma, QO, om) {
  k = nrow(Sigma)
  num = matrix(0, k, k)
  den = rep(0, k)
  for (h in 1:dim(A)[3]) {
    num = num + (A[, , h] %*% QO)^2
    den = den + diag(A[, , h] %*% Sigma %*% t(A[, , h]))
  }
  W = sweep(num, 2, om, "/")/den
  W/rowSums(W)
}
# BDY block orthogonalization of S for the clusters cl taken in the order ord
# (sequential projection on the earlier clusters); Q and Omega in the order of S
bdy = function(S, cl, ord) {
  k = nrow(S)
  idx = unlist(lapply(ord, function(c) which(cl == c)))
  Sp = S[idx, idx]
  clp = cl[idx]
  Qi = diag(k)
  for (m in seq_along(ord)[-1]) {
    b = which(clp == ord[m])
    prev = which(clp %in% ord[1:(m - 1)])
    Qi[b, prev] = -Sp[b, prev, drop = FALSE] %*% solve(Sp[prev, prev, drop = FALSE])
  }
  P = diag(k)[idx, , drop = FALSE]
  list(Q = t(P) %*% solve(Qi) %*% P, Om = t(P) %*% (Qi %*% Sp %*% t(Qi)) %*% P)
}
# one ordering: without controls, BDY on Sigma; with the control cluster 1, the
# residualization on cluster 1 and BDY on the Schur complement (2026b, eq. 7)
vd_order = function(A, Sigma, cl, ord, ctl) {
  k = nrow(Sigma)
  if (!ctl) {
    b = bdy(Sigma, cl, ord)
    return(vdq(A, Sigma, b$Q %*% b$Om, diag(b$Om)))
  }
  ic = which(cl == 1)
  ib = which(cl != 1)
  S11 = Sigma[ic, ic, drop = FALSE]
  G = Sigma[ib, ic, drop = FALSE] %*% solve(S11)
  Sbb1 = Sigma[ib, ib] - G %*% Sigma[ic, ib, drop = FALSE]
  b = bdy(Sbb1, cl[ib], ord)
  Q = matrix(0, k, k)
  Om = matrix(0, k, k)
  Q[ic, ic] = diag(length(ic))
  Q[ib, ic] = G
  Q[ib, ib] = b$Q
  Om[ic, ic] = S11
  Om[ib, ib] = b$Om
  stopifnot(max(abs(Q %*% Om %*% t(Q) - Sigma)) < 1e-8 * max(abs(Sigma)))
  vdq(A, Sigma, Q %*% Om, diag(Om))
}
allperms = function(v) {
  if (length(v) <= 1) return(list(v))
  out = list()
  for (i in seq_along(v)) for (p in allperms(v[-i])) out[[length(out) + 1]] = c(v[i], p)
  out
}
# net connectedness of each cluster in the cluster table of W (percent)
tnet = function(W, cl) {
  C = max(cl)
  T = matrix(NA, C, C)
  for (c in 1:C) for (d in 1:C) T[c, d] = sum(W[cl == c, cl == d, drop = FALSE])/sum(cl == c)
  (colSums(T) - diag(T)) - (rowSums(T) - diag(T))
}
# clustered decomposition averaged over the orderings (of the clusters other than the
# control cluster), with the positions and nets of every ordering
clavg = function(Phi, Sigma, H, cl, ctl, average = TRUE) {
  A = ma(Phi, H)
  C = max(cl)
  perm = if (ctl) 2:C else 1:C
  ords = if (average) allperms(perm) else list(perm)
  Ws = lapply(ords, function(o) vd_order(A, Sigma, cl, o, ctl))
  pos = t(sapply(ords, function(o) { p = rep(1, C); p[o] = seq_along(o) + ctl; p }))
  net = t(sapply(Ws, function(W) tnet(100*W, cl)))
  list(W = Reduce(`+`, Ws)/length(Ws), pos = matrix(pos, ncol = C), net = matrix(net, ncol = C))
}
# ordering diagnostics: averaged net, distribution, mean by position, first-versus-last
diag_tables = function(W, pos, net, cl, ctl) {
  C = max(cl)
  perm = if (ctl) 2:C else 1:C
  Cp = length(perm)
  N = cbind(tnet(W, cl), colMeans(net), apply(net, 2, sd), apply(net, 2, min), apply(net, 2, max))
  P = matrix(NA, Cp, Cp)
  D = matrix(NA, Cp, 5)
  for (ci in seq_along(perm)) {
    c = perm[ci]
    for (p in 1:Cp) P[ci, p] = mean(net[pos[, c] == p + ctl, c])
    others = setdiff(perm, c)
    key = apply(pos[, others, drop = FALSE], 1, function(r) paste(others[order(r)], collapse = "-"))
    f = which(pos[, c] == 1 + ctl)
    l = which(pos[, c] == C)
    d = net[f, c] - net[l[match(key[f], key[l])], c]
    D[ci, ] = c(mean(d), if (length(d) > 1) sd(d) else NA, min(d), max(d), length(d))
  }
  list(N = N, P = P, D = D)
}
dflsd = function(pos, net, cl, ctl) {
  C = max(cl)
  dg = diag_tables(matrix(0, length(cl), length(cl)), pos, net, cl, ctl)
  dfl = rep(NA, C)
  dfl[if (ctl) 2:C else 1:C] = dg$D[, 1]
  c(dfl, apply(net, 2, sd))
}
# simultaneous clustered orthogonalization; DC an arbitrary positive scaling
msqrt = function(M) {
  e = eigen((M + t(M))/2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(e$values), nrow(M)) %*% t(e$vectors)
}
sco_x = function(Sigma, cl, DC = NULL) {
  k = nrow(Sigma)
  d = sqrt(diag(Sigma))
  R1 = Sigma/outer(d, d)
  RC = R1 * outer(cl, cl, "==")
  R1h = msqrt(R1)
  R1ih = solve(R1h)
  Z = R1ih %*% msqrt(R1h %*% RC %*% R1h) %*% R1ih
  stopifnot(max(abs(Z %*% R1 %*% Z - RC)) < 1e-10)
  if (is.null(DC)) DC = diag(d)
  list(X = DC %*% Z %*% diag(1/d), R1 = R1, RC = RC)
}
scovd = function(Phi, Sigma, H, cl, DC = NULL) {
  A = ma(Phi, H)
  X = sco_x(Sigma, cl, DC)$X
  Q = solve(X)
  Om = X %*% Sigma %*% t(X)
  vdq(A, Sigma, Q %*% Om, diag(Om))
}

# measures of the 0.9.0 and 0.9.3 fixtures from one decomposition W (percent)
clrow = function(W, cl) {
  C = max(cl)
  k = nrow(W)
  own = com = ctg = to = rep(NA, C)
  T = matrix(NA, C, C)
  for (c in 1:C) {
    ic = which(cl == c)
    jc = which(cl != c)
    nc = length(ic)
    B = W[ic, ic, drop = FALSE]
    own[c] = sum(diag(B))/nc
    com[c] = (sum(B) - sum(diag(B)))/nc
    ctg[c] = if (length(jc)) sum(W[ic, jc])/nc else 0
    to[c] = if (length(jc)) sum(W[jc, ic])/nc else 0
    for (d in 1:C) T[c, d] = sum(W[ic, which(cl == d)])/nc
  }
  pr = c()
  for (c in 1:C) for (d in 1:C) if (d != c) pr = c(pr, T[c, d])
  same = outer(cl, cl, "==")
  wth = sum(W[same & !diag(k)])/k
  crs = sum(W[!same])/k
  tto = colSums(T) - diag(T)
  list(r090 = c(own, com, ctg, to, to - ctg, mean(own), mean(com), mean(ctg), pr),
       r093 = c(wth, crs, tto, tto - ctg), own = own, com = com, ctg = ctg, to = to, T = T,
       wth = wth, crs = crs, tto = tto)
}
names090 = function(C) {
  nm = c(paste0("cl_own_", 1:C), paste0("cl_comove_", 1:C), paste0("cl_contag_", 1:C),
         paste0("cl_to_", 1:C), paste0("cl_net_", 1:C), "cl_own", "cl_comove", "cl_contag")
  for (c in 1:C) for (d in 1:C) if (d != c) nm = c(nm, paste0("cl_", c, "_", d))
  nm
}
names093 = function(C) c("cl_within", "cl_cross", paste0("cl_tto_", 1:C), paste0("cl_tnet_", 1:C))

# all fixture files of a case. type: "avg" (average over orderings), "given", "sco";
# ctl: cluster 1 is the control cluster
clustfx = function(nm, Bt, Qt, H, cl, type = "avg", ctl = FALSE) {
  k = dim(Qt)[1]
  nt = dim(Qt)[3]
  C = max(cl)
  W = array(NA, c(k, k, nt))
  ts = c90 = c93 = c95 = NULL
  for (t in 1:nt) {
    B = if (length(dim(Bt)) == 3) Bt[, , t] else Bt
    Q = matrix(Qt[, , t], k, k)
    if (type == "sco") {
      W[, , t] = 100*scovd(B, Q, H, cl)
    } else {
      r = clavg(B, Q, H, cl, ctl, average = (type == "avg"))
      W[, , t] = 100*r$W
      if (type == "avg") {
        c95 = rbind(c95, dflsd(r$pos, r$net, cl, ctl))
        if (nt == 1) {
          wr(cbind(r$pos, r$net), paste0("R_", nm, "_clord.csv"))
          dg = diag_tables(100*r$W, r$pos, r$net, cl, ctl)
          Cp = nrow(dg$P)
          pad = function(M) cbind(M, matrix(NA, nrow(M), max(5, Cp) - ncol(M)))
          wr(rbind(pad(dg$N), pad(dg$P), pad(dg$D)), paste0("R_", nm, "_cldiag.csv"))
        }
      }
    }
    m = clrow(W[, , t], cl)
    ts = rbind(ts, permeas(W[, , t]/100))
    c90 = rbind(c90, m$r090)
    c93 = rbind(c93, m$r093)
  }
  colnames(ts) = measnames(k)
  colnames(c90) = names090(C)
  colnames(c93) = names093(C)
  per = seq_len(nt)
  write.csv(as.data.frame(cbind(period = per, ts)), file.path(out, paste0("R_", nm, "_ts.csv")), row.names = FALSE)
  avg = apply(W, 1:2, mean)
  wr(avg, paste0("R_", nm, "_fevd.csv"))
  write.csv(as.data.frame(cbind(period = per, c90)), file.path(out, paste0("R_", nm, "_clts.csv")), row.names = FALSE)
  write.csv(as.data.frame(cbind(period = per, c93)), file.path(out, paste0("R_", nm, "_clts093.csv")), row.names = FALSE)
  if (!is.null(c95) && nt > 1) {
    colnames(c95) = c(paste0("cl_dfl_", 1:C), paste0("cl_nsd_", 1:C))
    write.csv(as.data.frame(cbind(period = per, c95)), file.path(out, paste0("R_", nm, "_clts095.csv")),
              row.names = FALSE, na = "")
  }
  m = clrow(avg, cl)
  wr(rbind(cbind(m$own, m$com, m$ctg, m$to, m$to - m$ctg), c(mean(m$own), mean(m$com), mean(m$ctg), NA, NA)),
     paste0("R_", nm, "_cl.csv"))
  wr(m$T, paste0("R_", nm, "_ctg.csv"))
  wr(rbind(cbind(m$T, m$ctg), c(m$tto, NA), c(m$tto - m$ctg, NA)), paste0("R_", nm, "_cltab.csv"))
  wr(t(c(m$wth, m$crs)), paste0("R_", nm, "_clsys.csv"))
  cat(sprintf("%-26s %-5s periods=%4d mean TCI=%.6f net (table): %s\n", nm, type, nt, mean(ts[, "tci"]),
              paste(sprintf("%.3f", m$tto - m$ctg), collapse = " ")))
}
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

# ---------------------------------------------------------------------------
# checks of the corners and properties of simultaneous orthogonalization
x = getx("dy2012")
v = VAR(x, configuration = list(nlag = 4))
S = matrix(v$Q[, , 1], 4, 4)
A = ma(v$B, 10)
gr = TimeConnectedness(Phi = array(v$B, c(dim(v$B), 1)), Sigma = v$Q, nfore = 10, generalized = TRUE)$CT[, , 1]
cat(sprintf("SCO, one cluster, against TimeConnectedness(): max abs difference %.2e\n",
            max(abs(scovd(v$B, S, 10, rep(1, 4)) - gr))))
d = sqrt(diag(S))
Qf = diag(d) %*% msqrt(S/outer(d, d))
cat(sprintf("SCO, one variable per cluster, against Francis, Hansen and Tong Q = D1 R1^1/2: max abs difference %.2e\n",
            max(abs(scovd(v$B, S, 10, 1:4) - vdq(A, S, Qf, rep(1, 4))))))
cat(sprintf("SCO, D_C = D1 against D_C = diag(1, 2, 3, 4): max abs difference %.2e\n",
            max(abs(scovd(v$B, S, 10, c(1, 1, 2, 2)) - scovd(v$B, S, 10, c(1, 1, 2, 2), diag(1:4))))))
# the criterion sum_i corr(u_i, eps_i) at its maximum (2026c, Theorem 1)
crit = function(X, S) sum(diag(X %*% S)/sqrt(diag(S) * diag(X %*% S %*% t(X))))
x9 = getx("dy2009")
reg = c(US = 1, UK = 2, FRA = 2, GER = 2, HKG = 3, JPN = 3, AUS = 3, IDN = 3, KOR = 3, MYS = 3,
        PHL = 3, SGP = 3, TAI = 3, THA = 3, ARG = 1, BRA = 1, CHL = 1, MEX = 1, TUR = 2)
cl9 = unname(reg[colnames(x9)])
v9 = VAR(x9, configuration = list(nlag = 2))
S9 = matrix(v9$Q[, , 1], ncol(x9), ncol(x9))
sx = sco_x(S9, cl9)
cat(sprintf("criterion: SCO %.10f, trace formula (2026c, eq. 25) %.10f\n", crit(sx$X, S9),
            sum(sqrt(pmax(eigen(msqrt(sx$R1) %*% sx$RC %*% msqrt(sx$R1), symmetric = TRUE)$values, 0)))))
# a different admissible identification (rotated within the constraint) scores lower
U = qr.Q(qr(matrix(rnorm(ncol(x9)^2), ncol(x9))))
G = diag(sqrt(diag(S9))) %*% msqrt(sx$RC)
cat(sprintf("criterion of a rotated admissible identification: %.6f (lower)\n",
            crit(G %*% U %*% solve(G) %*% sx$X, S9)))
# the control-first decomposition by Schur complement = BDY with the control cluster first
clc = ifelse(colnames(x9) == "US", 1, cl9 + 1)
w1 = vd_order(ma(v9$B, 10), S9, clc, c(3, 2, 4), TRUE)
b = bdy(S9, clc, c(1, 3, 2, 4))
w2 = vdq(ma(v9$B, 10), S9, b$Q %*% b$Om, diag(b$Om))
cat(sprintf("control first: Schur complement against sequential BDY: max abs difference %.2e\n", max(abs(w1 - w2))))

# ---------------------------------------------------------------------------
# A. dy2009 weekly equity returns, static VAR(2), H = 10: US as a control cluster,
#    then the other Americas, Europe and Asia-Pacific (6 orderings); the same three
#    regions with the US in the Americas, without controls (6 orderings); and
#    simultaneous orthogonalization of the three regions
clustfx("dy2009_ctl", v9$B, staticQ(v9, ncol(x9)), 10, clc, "avg", TRUE)
clustfx("dy2009_ctl_given", v9$B, staticQ(v9, ncol(x9)), 10, clc, "given", TRUE)
clustfx("dy2009_reg_diag", v9$B, staticQ(v9, ncol(x9)), 10, cl9, "avg")
clustfx("dy2009_sco", v9$B, staticQ(v9, ncol(x9)), 10, cl9, "sco")
# B. seven clusters of dy2009 (5,040 orderings): US, Latin America, UK, euro area and
#    Turkey, Japan, East Asia, rest of Asia-Pacific
reg7 = c(US = 1, ARG = 2, BRA = 2, CHL = 2, MEX = 2, UK = 3, FRA = 4, GER = 4, TUR = 4, JPN = 5,
         HKG = 6, KOR = 6, TAI = 6, SGP = 6, AUS = 7, IDN = 7, MYS = 7, PHL = 7, THA = 7)
cl7 = unname(reg7[colnames(x9)])
clustfx("dy2009_seven", v9$B, staticQ(v9, ncol(x9)), 10, cl7, "avg")
# C. cgs2022 swap curves, static VAR(1), H = 10: six clusters (720 orderings), and
#    simultaneous orthogonalization
x = getx("cgs2022")
v = VAR(x, configuration = list(nlag = 1))
clustfx("cgs2022_diag", v$B, staticQ(v, ncol(x)), 10, rep(1:6, each = 3), "avg")
clustfx("cgs2022_sco", v$B, staticQ(v, ncol(x)), 10, rep(1:6, each = 3), "sco")
# D. rolling VAR(2) of dy2012, 200-observation windows, first 400 rows: SP500 as a
#    control cluster and the other three variables as clusters of their own (6
#    orderings); four singleton clusters (24 orderings); simultaneous orthogonalization
#    of two clusters
x = getx("dy2012")[1:400, ]
rv = rollvar(x, 2, 200)
clustfx("dy2012_ctl_roll", rv$B, rv$Q, 10, 1:4, "avg", TRUE)
clustfx("dy2012_four_roll", rv$B, rv$Q, 10, 1:4, "avg")
clustfx("dy2012_sco_roll", rv$B, rv$Q, 10, c(1, 1, 2, 2), "sco")
# E. TVP-VAR(1) of acg2020 (kappa2 0.96, Bayes prior from 200 observations, H = 12):
#    simultaneous orthogonalization of the European currencies and the yen, and the
#    ordering diagnostic of the same two clusters with GBP as a control cluster
x = getx("acg2020")
prior = BayesPrior(x = x, size = 199, nlag = 1)
tv = TVPVAR(x, configuration = list(l = c(0.99, 0.96), nlag = 1, prior = prior))
clustfx("acg2020_sco", tv$B_t, tv$Q_t, 12, c(1, 1, 2, 1), "sco")
clustfx("acg2020_ctl", tv$B_t, tv$Q_t, 12, c(2, 1, 3, 2), "avg", TRUE)
