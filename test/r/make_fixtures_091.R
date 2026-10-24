# Reference fixtures for version 0.9.1: bootstrap inference for VAR connectedness
# (Greenwood-Nimmo, Kocenda and Nguyen 2024). The bootstrap is implemented here
# independently of the Stata code, following the authors' Gauss programs
# (2.VAR_DY_bootstrap1.g, 2.VAR_DY_bootstrap2.g, 3.Prob_Events_Analysis_BAB_*.g):
#   - OLS VAR(p) with an intercept on each window; residuals recentred;
#   - recursive design: the first p observations of the window are kept, and
#     x_t = c + sum_j A_j x_{t-j} + e*_t with e*_t resampled with replacement;
#   - draws whose VAR has a companion root of modulus >= 1 are discarded and redrawn;
#   - a window whose own VAR is not stable uses the VAR of the most recent window
#     with a stable one (the authors carry those estimates forward);
#   - bias = mean of a first round of draws - estimate; a second round of draws minus
#     the bias gives the corrected draws (centred on the estimate); the reported
#     interval is their percentile interval recentred on the bias-corrected estimate,
#     estimate - bias (quantiles minus the bias);
#   - event probability = share of corrected draws at window re + j above the mean
#     corrected draw at window re - 1 (or of the m windows before re); with a
#     threshold, at least thr percent above (below with decrease).
# Random draws cannot be matched across programs, so the resampled rows are drawn
# here once, written to R_boot_idx_<case>.csv, and read by Stata through the
# test-only option testindices() of estat bootstrap: draw r of each window uses
# column r (discarded draws included; the second round continues from the first).
# The decomposition and the measures come from ConnectednessApproach
# (TimeConnectedness()).
#
# Run from the test/r folder:   Rscript make_fixtures_091.R

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
  as.matrix(x)
}
q7 = function(v, p) as.numeric(stats::quantile(v, p, type = 7))

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
# TCI, TO, FROM, NET from ConnectednessApproach
meas = function(B, Q, H, gen, all) {
  k = nrow(B)
  dca = TimeConnectedness(Phi = array(B, c(dim(B), 1)), Sigma = array(Q, c(k, k, 1)), nfore = H,
                          generalized = gen)
  if (!all) return(dca$TCI[1])
  c(dca$TCI[1], dca$TO[1, ], dca$FROM[1, ], dca$NET[1, ])
}
# R draws from the VAR v with initial values x0 (p x k); idx: index matrix;
# state$col: next column
draws = function(x0, v, p, H, gen, R, all, idx, state) {
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
    D = rbind(D, meas(vb$B, vb$Q, H, gen, all))
  }
  D
}
# window w of the data X (ws rows); a window whose VAR is not stable takes the VAR of
# the most recent window with a stable VAR (initial values stay those of window w)
bootwin = function(X, w, ws, p, H, gen, R, all, idx, level) {
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
  pt = meas(v$B, v$Q, H, gen, all)
  x0 = x[1:p, , drop = FALSE]
  D1 = draws(x0, v, p, H, gen, R, all, idx, state)
  bias = colMeans(D1) - pt
  D2 = sweep(draws(x0, v, p, H, gen, R, all, idx, state), 2, bias)
  a = (1 - level/100)/2
  list(tab = cbind(estimate = pt, bias = bias, corrected = pt - bias, se = apply(D2, 2, sd),
                   lower = apply(D2, 2, q7, a) - bias, upper = apply(D2, 2, q7, 1 - a) - bias),
       tci = D2[, 1], nd = state$nd)
}

# ---------------------------------------------------------------------------
# A. static VAR(2) of the first 800 observations of dy2012, horizon 10: generalized
#    (reps 40, 90% intervals) and orthogonalized (reps 30, 95%)
x = getx("dy2012")[1:800, ]
p = 2
n = nrow(x) - p
set.seed(2024)
idx = matrix(sample.int(n, n*100, replace = TRUE), n, 100)
write.table(idx, file.path(out, "R_boot_idx_static.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
for (cs in list(list(nm = "static_gen", gen = TRUE, R = 40, level = 90),
                list(nm = "static_orth", gen = FALSE, R = 30, level = 95))) {
  b = bootwin(x, 1, nrow(x), p, 10, cs$gen, cs$R, TRUE, idx, cs$level)
  write.table(b$tab, file.path(out, paste0("R_boot_", cs$nm, ".csv")), sep = ",", row.names = FALSE,
              col.names = FALSE)
  cat(sprintf("%-12s TCI %.6f bias %.6f corrected %.6f [%.4f, %.4f] discarded %d\n", cs$nm,
              b$tab[1, 1], b$tab[1, 2], b$tab[1, 3], b$tab[1, 5], b$tab[1, 6], b$nd))
}

# B. rolling VAR(2) of the first 230 observations of dy2012, windows of 200 rows
#    (Stata window(202)), horizon 10, reps 25, 95% bands; event probabilities for the
#    event at window 10 and the probabilities of a rise from the previous window
x = getx("dy2012")[1:230, ]
ws = 200
nw = nrow(x) - ws + 1
n = ws - p
set.seed(2025)
idx = matrix(sample.int(n, n*60, replace = TRUE), n, 60)
write.table(idx, file.path(out, "R_boot_idx_roll.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
ts = NULL
TD = matrix(NA, 25, nw)
for (w in 1:nw) {
  b = bootwin(x, w, ws, p, 10, TRUE, 25, FALSE, idx, 95)
  ts = rbind(ts, b$tab[1, ])
  TD[, w] = b$tci
}
colnames(ts) = c("tci", "tci_bias", "tci_bc", "tci_se", "tci_lo", "tci_hi")
prob = c(NA, sapply(2:nw, function(w) 100*mean(TD[, w] > mean(TD[, w - 1]))))
# window 1 has no previous window (prob missing), so the file starts at window 2
write.csv(data.frame(period = 1:nw, ts, prob = prob)[-1, ], file.path(out, "R_boot_roll_ts.csv"),
          row.names = FALSE)
evp = function(re, leads, m, thr, dec) {
  sb = mean(colMeans(TD[, (re - m):(re - 1), drop = FALSE]))
  sapply(leads, function(j) {
    d = TD[, re + j]
    if (thr > 0) {
      d = 100*(d - sb)/sb
      100*mean(if (dec) d <= -thr else d >= thr)
    } else 100*mean(if (dec) d < sb else d > sb)
  })
}
ev = rbind(evp(10, c(0, 1, 5), 1, 0, FALSE), evp(10, c(0, 1, 5), 3, 2, TRUE))
write.table(ev, file.path(out, "R_boot_roll_event.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
cat(sprintf("roll         windows %d, mean TCI %.6f, mean bias %.6f; event: %s\n", nw, mean(ts[, 1]),
            mean(ts[, 2]), paste(round(ev[1, ], 1), collapse = " ")))
