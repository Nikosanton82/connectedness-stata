# Monte Carlo coverage of the bootstrap intervals of the spillover asymmetry measure
# (version 0.9.4). Under a symmetric data generating process the positive and negative
# parts of the series are mirror images in distribution (-min(x, 0) = max(-x, 0) and
# -x has the distribution of x), and connectedness does not change when all series
# change sign, so the SAM is zero. The joint bootstrap of estat bootstrap, sam is
# replicated here (independent code, own generalized decomposition, no package
# needed) on M simulated samples, and the share of intervals that contain zero is
# reported for each measure of the SAM table, and, as a benchmark, for the TCI of the
# total system (a correctly specified VAR) against its population value.
#
# DGP: Gaussian VAR(1), k = 3, T = 500, no intercept, B = 0.3 I + 0.1 (J - I),
# residual correlations 0.4. Bootstrap: VAR(1), horizon 10, 199 draws per round, two
# rounds (bias correction), 90% and 95% intervals.
#
# Run:   Rscript sam_coverage_094.R [M] [cores] [T]

args = commandArgs(TRUE)
M = if (length(args) > 0) as.integer(args[1]) else 200
cores = if (length(args) > 1) as.integer(args[2]) else 1
R = 199
H = 10
p = 1
k = 3
TT = if (length(args) > 2) as.integer(args[3]) else 500

q7 = function(v, pr) as.numeric(stats::quantile(v, pr, type = 7))
varc = function(x, p) {
  n = nrow(x) - p
  Y = x[(p + 1):nrow(x), , drop = FALSE]
  Z = do.call(cbind, lapply(1:p, function(j) x[(p + 1 - j):(nrow(x) - j), , drop = FALSE]))
  fit = lm.fit(cbind(1, Z), Y)
  cf = fit$coefficients
  E = sweep(fit$residuals, 2, colMeans(fit$residuals))
  list(B = t(cf[-1, , drop = FALSE]), c = cf[1, ], E = E, Q = crossprod(E)/n)
}
maxroot = function(B) {
  kk = nrow(B)
  pp = ncol(B)/kk
  F = if (pp == 1) B else rbind(B, cbind(diag(kk*(pp - 1)), matrix(0, kk*(pp - 1), kk)))
  max(Mod(eigen(F, only.values = TRUE)$values))
}
gfevd = function(B, Q, H) {
  kk = nrow(B)
  pp = ncol(B)/kk
  A = list(diag(kk))
  for (h in 2:H) {
    Ah = matrix(0, kk, kk)
    for (j in 1:min(pp, h - 1)) Ah = Ah + B[, ((j - 1)*kk + 1):(j*kk), drop = FALSE] %*% A[[h - j]]
    A[[h]] = Ah
  }
  num = Reduce(`+`, lapply(A, function(a) (a %*% Q)^2))
  den = Reduce(`+`, lapply(A, function(a) diag(a %*% Q %*% t(a))))
  th = sweep(num, 2, diag(Q), "/") / den
  100*th/rowSums(th)
}
meas1 = function(v) {
  ct = gfevd(v$B, v$Q, H)
  to = colSums(ct) - diag(ct)
  fr = rowSums(ct) - diag(ct)
  c(sum(fr)/k, to, fr, to - fr)
}
sam = function(vs) c(meas1(vs[[2]]) - meas1(vs[[3]]), meas1(vs[[1]])[1])
simrec = function(x0, v, e) {
  n = length(e)
  Xb = rbind(x0, matrix(NA, n, k))
  Bt = t(v$B)
  for (t in (p + 1):(p + n)) Xb[t, ] = v$c + v$E[e[t - p], ] + Xb[t - 1, ] %*% Bt
  Xb
}
draws = function(x0s, vs, R, nd) {
  n = nrow(vs[[1]]$E)
  D = NULL
  while (is.null(D) || nrow(D) < R) {
    e = sample.int(n, n, replace = TRUE)
    vb = list()
    ok = TRUE
    for (s in 1:3) {
      vb[[s]] = varc(simrec(x0s[[s]], vs[[s]], e), p)
      if (maxroot(vb[[s]]$B) >= 1) { ok = FALSE; break }
    }
    if (!ok) next
    D = rbind(D, sam(vb))
  }
  D
}
B0 = 0.3*diag(k) + 0.1*(1 - diag(k))
S0 = 0.4 + 0.6*diag(k)
L0 = t(chol(S0))
tci0 = meas1(list(B = B0, Q = S0))[1]
one = function(m) {
  set.seed(9400 + m)
  u = matrix(rnorm((TT + 100)*k), ncol = k) %*% t(L0)
  x = matrix(0, TT + 100, k)
  for (t in 2:(TT + 100)) x[t, ] = B0 %*% x[t - 1, ] + u[t, ]
  x = x[101:(TT + 100), ]
  Xs = list(x, pmax(x, 0), pmin(x, 0))
  vs = lapply(Xs, varc, p = p)
  x0s = lapply(Xs, function(X) X[1:p, , drop = FALSE])
  pt = sam(vs)
  D1 = draws(x0s, vs, R)
  bias = colMeans(D1) - pt
  D2 = sweep(draws(x0s, vs, R), 2, bias)
  res = NULL
  for (lv in c(90, 95)) {
    a = (1 - lv/100)/2
    lo = apply(D2, 2, q7, a) - bias
    hi = apply(D2, 2, q7, 1 - a) - bias
    res = c(res, as.numeric(lo <= c(rep(0, length(lo) - 1), tci0) & hi >= c(rep(0, length(lo) - 1), tci0)))
  }
  c(res, pt, pt - bias)
}
t0 = Sys.time()
out = if (cores > 1) parallel::mclapply(1:M, one, mc.cores = cores) else lapply(1:M, one)
out = do.call(rbind, out)
m3 = 2 + 3*k
lab = c("SAM", paste0("TO_", 1:k), paste0("FROM_", 1:k), paste0("NET_", 1:k), "TCI")
cov90 = colMeans(out[, 1:m3])
cov95 = colMeans(out[, (m3 + 1):(2*m3)])
est = out[, (2*m3 + 1):(3*m3)]
bc = out[, (3*m3 + 1):(4*m3)]
cat(sprintf("SAM bootstrap coverage under a symmetric DGP: M = %d samples, T = %d, k = %d, R = %d per round\n",
            M, TT, k, R))
cat(sprintf("%-8s %9s %9s %11s %11s\n", "measure", "cover90", "cover95", "mean est", "mean bc"))
cat(sprintf("(TCI: the total system, against its population value %.4f)\n", tci0))
for (i in 1:m3) cat(sprintf("%-8s %9.3f %9.3f %11.4f %11.4f\n", lab[i], cov90[i], cov95[i], mean(est[, i]), mean(bc[, i])))
cat(sprintf("Monte Carlo standard error of a coverage of 0.90: %.3f; of 0.95: %.3f\n",
            sqrt(0.9*0.1/M), sqrt(0.95*0.05/M)))
cat("time:", format(Sys.time() - t0), "\n")
