# Reference fixtures for version 0.9.2: bootstrap bands for the pairwise connectedness
# table and net pairwise directional connectedness (NPDC), and event and detection
# probabilities for any measure (option stat() of estat bootstrap).
# The bootstrap is the independent implementation of make_fixtures_091.R (Greenwood-
# Nimmo, Kocenda and Nguyen 2024), extended to keep, in every draw, all the measures of
# the decomposition: TCI, TO, FROM, NET, the pairwise table CT (percent; element (i, j)
# is the spillover from j to i) and the NPDC matrix of ConnectednessApproach
# (NPDC = CT - t(CT); element (i, j) is the net spillover from j to i).
# The resampled rows are those written by make_fixtures_091.R (R_boot_idx_*.csv), read
# here and by Stata (test-only option testindices()).
#
# Run from the test/r folder, after make_fixtures_091.R:   Rscript make_fixtures_092.R

if (!exists("ConnectednessApproach")) {
  lib = Sys.getenv("R_LIBS_USER")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(ConnectednessApproach))
  cat("ConnectednessApproach version:", as.character(packageVersion("ConnectednessApproach")), "\n")
}
suppressPackageStartupMessages(library(zoo))

out = file.path("..", "fixtures")

getx = function(d) {
  x = get(data(list = d))
  if (!is(x, "zoo")) x = zoo::zoo(as.matrix(x[, -1]), as.Date(x[, 1]))
  as.matrix(x)
}
q7 = function(v, p) as.numeric(stats::quantile(v, p, type = 7))
readidx = function(f) unname(as.matrix(read.csv(file.path(out, f), header = FALSE)))

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
# all measures from ConnectednessApproach: TCI, TO, FROM, NET, CT and NPDC by rows
meas = function(B, Q, H, gen) {
  k = nrow(B)
  dca = TimeConnectedness(Phi = array(B, c(dim(B), 1)), Sigma = array(Q, c(k, k, 1)), nfore = H,
                          generalized = gen)
  ct = 100*dca$CT[, , 1]
  np = dca$NPDC[, , 1]
  stopifnot(max(abs(np - (ct - t(ct)))) < 1e-10)
  c(dca$TCI[1], dca$TO[1, ], dca$FROM[1, ], dca$NET[1, ], as.vector(t(ct)), as.vector(t(np)))
}
draws = function(x0, v, p, H, gen, R, idx, state) {
  k = ncol(x0)
  n = nrow(v$E)
  D = NULL
  while (is.null(D) || nrow(D) < R) {
    e = idx[, state$col]
    state$col = state$col + 1
    Xb = rbind(x0, matrix(NA, n, k))
    for (t in (p + 1):(p + n)) {
      val = v$c + v$E[e[t - p], ]
      for (j in 1:p) val = val + v$B[, ((j - 1)*k + 1):(j*k), drop = FALSE] %*% Xb[t - j, ]
      Xb[t, ] = val
    }
    vb = varc(Xb, p)
    if (maxroot(vb$B) >= 1) {
      state$nd = state$nd + 1
      next
    }
    D = rbind(D, meas(vb$B, vb$Q, H, gen))
  }
  D
}
bootwin = function(X, w, ws, p, H, gen, R, idx, level) {
  state = new.env()
  state$col = 1
  state$nd = 0
  x = X[w:(w + ws - 1), , drop = FALSE]
  v = varc(x, p)
  s = w
  while (maxroot(v$B) >= 1 && s > 1) {
    s = s - 1
    v = varc(X[s:(s + ws - 1), , drop = FALSE], p)
  }
  pt = meas(v$B, v$Q, H, gen)
  x0 = x[1:p, , drop = FALSE]
  D1 = draws(x0, v, p, H, gen, R, idx, state)
  bias = colMeans(D1) - pt
  D2 = sweep(draws(x0, v, p, H, gen, R, idx, state), 2, bias)
  a = (1 - level/100)/2
  list(tab = cbind(estimate = pt, bias = bias, corrected = pt - bias, se = apply(D2, 2, sd),
                   lower = apply(D2, 2, q7, a) - bias, upper = apply(D2, 2, q7, 1 - a) - bias),
       D = D2)
}
k = 4
ipw = function(i, j) 1 + 3*k + (i - 1)*k + j          # CT[i, j]: from j to i
inp = function(i, j) 1 + 3*k + k*k + (i - 1)*k + j    # NPDC[i, j]: net from j to i
wr = function(m, f) write.table(m, file.path(out, f), sep = ",", row.names = FALSE, col.names = FALSE, na = "")

# ---------------------------------------------------------------------------
# A. static VAR(2) of the first 800 observations of dy2012, horizon 10, generalized,
#    reps 40, 90% intervals: the pairwise table and the NPDC matrix (k*k rows each, by
#    rows of the matrix; columns estimate, bias, corrected, se, lower, upper)
x = getx("dy2012")[1:800, ]
p = 2
idx = readidx("R_boot_idx_static.csv")
b = bootwin(x, 1, nrow(x), p, 10, TRUE, 40, idx, 90)
wr(b$tab[ipw(1, 1):ipw(k, k), ], "R_boot_static_pw.csv")
wr(b$tab[inp(1, 1):inp(k, k), ], "R_boot_static_npdc.csv")
cat(sprintf("static       pairwise SP500 <- USDX %.6f (corrected %.6f); NPDC USDX -> SP500 %.6f [%.4f, %.4f]\n",
            b$tab[ipw(1, 4), 1], b$tab[ipw(1, 4), 3], b$tab[inp(1, 4), 3], b$tab[inp(1, 4), 5],
            b$tab[inp(1, 4), 6]))

# B. rolling VAR(2) of the first 230 observations, windows of 200 rows (Stata
#    window(202)), horizon 10, reps 25, 95% bands: long file, one row per window and
#    ordered pair (from a, to b), with the pairwise spillover CT[b, a] and NPDC[b, a]
x = getx("dy2012")[1:230, ]
ws = 200
nw = nrow(x) - ws + 1
idx = readidx("R_boot_idx_roll.csv")
TD = array(NA, c(25, 1 + 3*k + 2*k*k, nw))
lng = NULL
for (w in 1:nw) {
  b = bootwin(x, w, ws, p, 10, TRUE, 25, idx, 95)
  TD[, , w] = b$D
  for (a in 1:k) for (bb in 1:k) {
    if (a == bb) next
    lng = rbind(lng, c(w, a, bb, b$tab[ipw(bb, a), ], b$tab[inp(bb, a), ]))
  }
}
wr(lng, "R_boot_roll_pw.csv")

# event probabilities for the event at window 10 (Stata event(209)) for a measure m
evp = function(m, re, leads, ref, thr, dec) {
  S = TD[, m, ]
  sb = mean(colMeans(S[, (re - ref):(re - 1), drop = FALSE]))
  sapply(leads, function(j) {
    d = S[, re + j]
    if (thr > 0) {
      d = 100*(d - sb)/sb
      100*mean(if (dec) d <= -thr else d >= thr)
    } else 100*mean(if (dec) d < sb else d > sb)
  })
}
# detection probabilities: rise from the previous window (of at least thr percent)
detp = function(m, thr) {
  S = TD[, m, ]
  c(NA, sapply(2:nw, function(w) {
    sb = mean(S[, w - 1])
    d = S[, w]
    if (thr > 0) 100*mean(100*(d - sb)/sb >= thr) else 100*mean(d > sb)
  }))
}
# variables: 1 SP500, 2 R_10Y, 3 DJUBSCOM, 4 USDX
ev = rbind(evp(inp(2, 1), 10, c(0, 1, 5), 1, 0, FALSE),   # stat(npdc SP500 R_10Y)
           evp(ipw(1, 4), 10, c(0, 1, 5), 3, 2, FALSE),   # stat(pairwise USDX SP500) reference(3) threshold(2)
           evp(1 + 2*k + 3, 10, c(0, 1, 5), 1, 0, TRUE),  # stat(net DJUBSCOM) decrease
           evp(1 + 2, 10, c(0, 1, 5), 1, 1, FALSE))       # stat(to R_10Y) threshold(1)
wr(ev, "R_boot_roll_statevent.csv")
dt = cbind(detp(inp(2, 1), 0), detp(1 + k + 4, 2))        # stat(npdc SP500 R_10Y) detect(0); stat(from USDX) detect(2)
wr(dt, "R_boot_roll_statdetect.csv")
cat(sprintf("roll         %d windows x %d pairs; events: %s\n", nw, k*(k - 1),
            paste(apply(round(ev, 1), 1, paste, collapse = " "), collapse = " | ")))
cat(sprintf("detect       npdc flagged (>= 90): %d; from USDX flagged: %d\n",
            sum(dt[, 1] >= 90, na.rm = TRUE), sum(dt[, 2] >= 90, na.rm = TRUE)))
