# Reference fixtures for version 0.9.4: bootstrap bands for the spillover asymmetry
# measure (SAM = TCI of the positive system minus TCI of the negative system; Barunik,
# Kocenda and Vacha 2016) and its TO, FROM and NET, from the joint bootstrap of the
# total, positive and negative systems (estat bootstrap, option sam and stat(sam)).
#
# The bootstrap is the independent implementation of make_fixtures_091.R and
# make_fixtures_092.R (Greenwood-Nimmo, Kocenda and Nguyen 2024), extended to three
# VARs: every draw resamples the same rows of the residuals of the three systems,
# rebuilds the three series recursively from their own initial values, and is
# discarded if any of the three VARs is not stable. Bias correction applies to each
# measure, and so to the difference. A window whose VAR is not stable takes, system by
# system, the estimates of the most recent window with a stable VAR.
# The measures come from TimeConnectedness() of ConnectednessApproach; a second, own
# generalized decomposition is checked against it.
# The resampled rows are drawn here, written to R_boot_idx_sam_*.csv and read by Stata
# (test-only option testindices()).
#
# Run from the test/r folder:   Rscript make_fixtures_094.R

lib = Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ConnectednessApproach))
cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
suppressPackageStartupMessages(library(zoo))

out = file.path("..", "fixtures")

getx = function(d) {
  x = get(data(list = d))
  if (!is(x, "zoo")) x = zoo::zoo(as.matrix(x[, -1]), as.Date(x[, 1]))
  as.matrix(x)
}
q7 = function(v, p) as.numeric(stats::quantile(v, p, type = 7))
wr = function(m, f) write.table(m, file.path(out, f), sep = ",", row.names = FALSE, col.names = FALSE, na = "")

# OLS VAR with intercept: B (k x kp), intercept, recentred residuals, covariance E'E/n
varc = function(x, p) {
  k = ncol(x)
  n = nrow(x) - p
  Y = x[(p + 1):nrow(x), , drop = FALSE]
  Z = do.call(cbind, lapply(1:p, function(j) x[(p + 1 - j):(nrow(x) - j), , drop = FALSE]))
  fit = lm.fit(cbind(1, Z), Y)
  cf = fit$coefficients
  E = sweep(fit$residuals, 2, colMeans(fit$residuals))
  list(B = t(cf[-1, , drop = FALSE]), c = cf[1, ], E = E, Q = crossprod(E)/n)
}
maxroot = function(B) {
  k = nrow(B)
  p = ncol(B)/k
  F = if (p == 1) B else rbind(B, cbind(diag(k*(p - 1)), matrix(0, k*(p - 1), k)))
  max(Mod(eigen(F, only.values = TRUE)$values))
}
# own generalized decomposition (Pesaran and Shin 1998), rows normalized, percent
gfevd = function(B, Q, H) {
  k = nrow(B)
  p = ncol(B)/k
  A = list(diag(k))
  for (h in 2:H) {
    Ah = matrix(0, k, k)
    for (j in 1:min(p, h - 1)) Ah = Ah + B[, ((j - 1)*k + 1):(j*k), drop = FALSE] %*% A[[h - j]]
    A[[h]] = Ah
  }
  num = Reduce(`+`, lapply(A, function(a) (a %*% Q)^2))
  den = Reduce(`+`, lapply(A, function(a) diag(a %*% Q %*% t(a))))
  th = sweep(num, 2, diag(Q), "/") / den
  100*th/rowSums(th)
}
own = 0
# all measures of one system: TCI, TO, FROM, NET (cn_bootfull() up to NET)
meas1 = function(B, Q, H) {
  k = nrow(B)
  dca = TimeConnectedness(Phi = array(B, c(dim(B), 1)), Sigma = array(Q, c(k, k, 1)), nfore = H,
                          generalized = TRUE)
  ct = 100*dca$CT[, , 1]
  own <<- max(own, max(abs(ct - gfevd(B, Q, H))))
  c(dca$TCI[1], dca$TO[1, ], dca$FROM[1, ], dca$NET[1, ])
}
# measures of the three systems (vs: list of the total, positive and negative VARs):
# TCI, TO, FROM, NET of the total system; SAM, its TO, FROM and NET; TCI positive,
# TCI negative
meas3 = function(vs, H) {
  m = lapply(vs, function(v) meas1(v$B, v$Q, H))
  c(m[[1]], m[[2]] - m[[3]], m[[2]][1], m[[3]][1])
}
draws3 = function(x0s, vs, p, H, R, idx, state) {
  k = ncol(x0s[[1]])
  n = nrow(vs[[1]]$E)
  D = NULL
  while (is.null(D) || nrow(D) < R) {
    e = idx[, state$col]
    state$col = state$col + 1
    vb = list()
    ok = TRUE
    for (s in 1:3) {
      v = vs[[s]]
      Xb = rbind(x0s[[s]], matrix(NA, n, k))
      for (t in (p + 1):(p + n)) {
        val = v$c + v$E[e[t - p], ]
        for (j in 1:p) val = val + v$B[, ((j - 1)*k + 1):(j*k), drop = FALSE] %*% Xb[t - j, ]
        Xb[t, ] = val
      }
      vb[[s]] = varc(Xb, p)
      if (maxroot(vb[[s]]$B) >= 1) {
        ok = FALSE
        break
      }
    }
    if (!ok) {
      state$nd = state$nd + 1
      next
    }
    D = rbind(D, meas3(vb, H))
  }
  D
}
# bootstrap of window w (rows w to w + ws - 1) of the series Xs (list of three
# matrices: total, positive, negative)
bootwin3 = function(Xs, w, ws, p, H, R, idx, level) {
  state = new.env()
  state$col = 1
  state$nd = 0
  vs = list()
  x0s = list()
  carry = FALSE
  for (s in 1:3) {
    X = Xs[[s]]
    v = varc(X[w:(w + ws - 1), , drop = FALSE], p)
    u = w
    while (maxroot(v$B) >= 1 && u > 1) {
      u = u - 1
      v = varc(X[u:(u + ws - 1), , drop = FALSE], p)
    }
    if (u < w) carry = TRUE
    vs[[s]] = v
    x0s[[s]] = X[w:(w + p - 1), , drop = FALSE]
  }
  pt = meas3(vs, H)
  D1 = draws3(x0s, vs, p, H, R, idx, state)
  bias = colMeans(D1) - pt
  D2 = sweep(draws3(x0s, vs, p, H, R, idx, state), 2, bias)
  a = (1 - level/100)/2
  list(tab = cbind(estimate = pt, bias = bias, corrected = pt - bias, se = apply(D2, 2, sd),
                   lower = apply(D2, 2, q7, a) - bias, upper = apply(D2, 2, q7, 1 - a) - bias,
                   prob = 100*colMeans(sweep(D2, 2, bias) > 0)),
       D = D2, nd = state$nd, carry = carry)
}
k = 4
m3 = 1 + 3*k
isam = (m3 + 1):(2*m3)                      # SAM, its TO, FROM, NET
ipos = 2*m3 + 1
ineg = 2*m3 + 2

g = getx("g2020")
cat("g2020:", nrow(g), "rows;", colnames(g), "; first row", round(g[1, ], 6), "\n")

# ---------------------------------------------------------------------------
# A. static VAR(2) of the first 600 observations of g2020 (EUR GBP CHF JPY), option
#    asymmetric (positive and negative parts of the returns), horizon 10, reps 40, 90%
#    intervals. Files: the total system (TCI, TO, FROM, NET: estimate, bias, corrected,
#    se, lower, upper), the SAM table (and prob) and the TCI of each part
x = g[1:600, ]
p = 2
n = nrow(x) - p
set.seed(2094)
idx = matrix(sample.int(n, n*150, replace = TRUE), n, 150)
wr(idx, "R_boot_idx_sam_static.csv")
Xs = list(x, pmax(x, 0), pmin(x, 0))
b = bootwin3(Xs, 1, nrow(x), p, 10, 40, idx, 90)
wr(b$tab[1:m3, 1:6], "R_sam_static_tot.csv")
wr(b$tab[isam, ], "R_sam_static_sam.csv")
wr(b$tab[c(ipos, ineg), 1:6], "R_sam_static_parts.csv")
cat(sprintf("static       SAM %.6f corrected %.6f [%.4f, %.4f] P(>0) %.1f; TCI+ %.4f TCI- %.4f; discarded %d\n",
            b$tab[isam[1], 1], b$tab[isam[1], 3], b$tab[isam[1], 5], b$tab[isam[1], 6], b$tab[isam[1], 7],
            b$tab[ipos, 1], b$tab[ineg, 1], b$nd))

# A2. the same with 95% intervals, no bias correction is not tested here (as 0.9.1);
#     orthogonalized decompositions are covered by the Stata identities

# ---------------------------------------------------------------------------
# B. rolling VAR(1) with positive() and negative(): squared returns (total), squared
#    positive and squared negative returns (stand-ins for realized semivariances),
#    first 240 observations, windows of 200 rows (Stata window(201)), horizon 10, reps
#    25, 95% bands. File: one row per window with the bands of the TCI of the total
#    system, SAM (with prob), SAM TO, FROM, NET of each variable (with prob), TCI+ and
#    TCI-, then the discarded draws and the carry-forward flag
x = g[1:240, ]
p = 1
ws = 200
nw = nrow(x) - ws + 1
n = ws - p
set.seed(2095)
idx = matrix(sample.int(n, n*80, replace = TRUE), n, 80)
wr(idx, "R_boot_idx_sam_roll.csv")
Xs = list(x^2, pmax(x, 0)^2, pmin(x, 0)^2)
TD = array(NA, c(25, 2*m3 + 2, nw))
rows = NULL
for (w in 1:nw) {
  b = bootwin3(Xs, w, ws, p, 10, 25, idx, 95)
  TD[, , w] = b$D
  r = c(b$tab[1, 1:6])
  for (i in isam) r = c(r, b$tab[i, ])
  r = c(r, b$tab[ipos, 1:6], b$tab[ineg, 1:6], b$nd, !b$carry)
  rows = rbind(rows, r)
}
wr(rows, "R_sam_roll.csv")
cat(sprintf("roll         %d windows; SAM corrected first/last %.4f / %.4f; windows excluding zero %d; discarded %d; carried %d\n",
            nw, rows[1, 9], rows[nw, 9], sum(rows[, 11] > 0 | rows[, 12] < 0), sum(rows[, ncol(rows) - 1]),
            sum(rows[, ncol(rows)] == 0)))

# event probabilities for the event at window 10 (Stata event(209)) for measure m
evp = function(m, re, leads, ref, dec) {
  S = TD[, m, ]
  sb = mean(colMeans(S[, (re - ref):(re - 1), drop = FALSE]))
  sapply(leads, function(j) {
    d = S[, re + j]
    100*mean(if (dec) d < sb else d > sb)
  })
}
detp = function(m) {
  S = TD[, m, ]
  c(NA, sapply(2:nw, function(w) 100*mean(S[, w] > mean(S[, w - 1]))))
}
# variables: 1 EUR, 2 GBP, 3 CHF, 4 JPY
ev = rbind(evp(isam[1], 10, c(0, 1, 5), 1, FALSE),             # stat(sam)
           evp(isam[1 + 2*k + 3], 10, c(0, 1, 5), 3, TRUE),    # stat(sam net CHF) reference(3) decrease
           evp(isam[1 + 1], 10, c(0, 1, 5), 1, FALSE))         # stat(sam to EUR)
wr(ev, "R_sam_roll_event.csv")
dt = cbind(detp(isam[1]), detp(isam[1 + k + 2]))              # stat(sam) detect(0); stat(sam from GBP) detect(0)
wr(dt, "R_sam_roll_detect.csv")
cat(sprintf("events       %s\n", paste(apply(round(ev, 1), 1, paste, collapse = " "), collapse = " | ")))
cat(sprintf("own generalized decomposition against TimeConnectedness(): max abs difference %.2e\n", own))
