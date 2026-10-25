# Reference fixtures for version 0.10.0: BEKK(1,1) spot variance spillovers (Fengler
# and Herwartz 2018) and volatility impulse responses (Hafner and Herwartz 2006).
#
# The BEKK models are estimated with the R package BEKKs (Fülle, Lange, Hafner and
# Herwartz; bekk_fit(), BHHH with analytic scores) on the demeaned returns, from the
# "simple" starting values; standard errors from the outer product of its analytic
# scores. The spillovers are computed here from the fitted parameters, independently
# of the Stata code and of the formulas Stata uses:
#   - Cov(vech(e e') | H) from the moments of e = H^(1/2) xi: Cov(e_i e_j, e_k e_l) =
#     H_ik H_jl + H_il H_jk + sum_a (kappa_a - 3) B_ia B_ja B_ka B_la, B = H^(1/2)
#     (instead of the Kronecker form of Proposition 1);
#   - Theta_m from the closed form of Theta(L) = (I - A_sum L)^-1 (I - B L):
#     Theta_0 = I, Theta_m = A_sum^m - A_sum^(m-1) B (instead of the recursion);
#   - Psi_m = Theta_m Sigma_(t+M-m|t)^(1/2), lambda_ij = sum_m Psi_ij^2 / row sum.
# The VIRFs come from BEKKs' virf(). Stata reads the parameters (bekkparams()) and is
# compared exactly. Standard errors come from the outer product of exact scores (own
# code, see score_exact()); BEKKs' scores carry a derivative of the first covariance
# H_1 = e'e/T, which does not depend on the parameters. BEKKs' iterations can stop
# before the maximum (its positive diagonal of C); Stata's own estimates are compared
# with the maximum found here (mle_exact()).
#
# Run from the test/r folder:   Rscript make_fixtures_0100.R

lib = Sys.getenv("R_LIBS_USER")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(ConnectednessApproach)
  library(BEKKs)
  library(zoo)
})
cat("ConnectednessApproach", as.character(packageVersion("ConnectednessApproach")),
    "; BEKKs", as.character(packageVersion("BEKKs")), "\n")
out = file.path("..", "fixtures")
wr = function(m, f) write.table(m, file.path(out, f), sep = ",", row.names = FALSE, col.names = FALSE, na = ".")
wcsv = function(df, f) write.csv(as.data.frame(df), file.path(out, f), row.names = FALSE)
getx = function(d) {
  x = get(data(list = d))
  if (!is(x, "zoo")) x = zoo::zoo(as.matrix(x[, -1]), as.Date(x[, 1]))
  as.matrix(x)
}
thin = function(m) {
  nt = nrow(m)
  m[(m[, "period"] %% 50 == 1) | (m[, "period"] == nt), , drop = FALSE]
}
pairs = function(k) {
  pr = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) pr = c(pr, paste0(i, "_", j))
  pr
}
pairmat = function(A, lower) {
  k = nrow(A)
  o = c()
  for (i in 1:(k - 1)) for (j in (i + 1):k) o = c(o, if (lower) A[j, i] else A[i, j])
  o
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

# parameters (BEKKs order) -> C, A, G
unpack = function(th, k, type) {
  C = matrix(0, k, k)
  C[lower.tri(C, diag = TRUE)] = th[1:(k*(k + 1)/2)]
  m = k*(k + 1)/2
  if (type == "bekk") {
    A = matrix(th[(m + 1):(m + k*k)], k, k)
    G = matrix(th[(m + k*k + 1):(m + 2*k*k)], k, k)
  } else if (type == "dbekk") {
    A = diag(th[(m + 1):(m + k)])
    G = diag(th[(m + k + 1):(m + 2*k)])
  } else {
    A = sqrt(th[m + 1])*diag(k)
    G = sqrt(th[m + 2])*diag(k)
  }
  list(C = C, A = A, G = G)
}
start = function(e, type) {
  k = ncol(e)
  U = crossprod(e)/nrow(e)
  C0 = diag(0.05*diag(U))
  for (i in 1:k) for (j in 1:k) if (i != j) C0[i, j] = U[i, j]/sqrt(U[i, i]*U[j, j])*sqrt(C0[i, i]*C0[j, j])
  L = t(chol(C0))
  th = L[lower.tri(L, diag = TRUE)]
  if (type == "bekk") c(th, c(0.3*diag(k)), c(0.92*diag(k)))
  else if (type == "dbekk") c(th, rep(0.3, k), rep(0.92, k))
  else c(th, 0.05, 0.9)
}
msqrt = function(X) {
  e = eigen((X + t(X))/2, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(X)) %*% t(e$vectors)
}
vechidx = function(k) which(lower.tri(diag(k), diag = TRUE), arr.ind = TRUE)  # (row, col), vech order
dupl = function(k) {
  ix = vechidx(k)
  D = matrix(0, k*k, nrow(ix))
  for (c in 1:nrow(ix)) { i = ix[c, 1]; j = ix[c, 2]; D[(j - 1)*k + i, c] = 1; D[(i - 1)*k + j, c] = 1 }
  D
}
# Cov(vech(e e')) for e = H^(1/2) xi, independent xi with fourth moments kap
sigvech = function(H, kap) {
  k = nrow(H)
  B = msqrt(H)
  ix = vechidx(k)
  n = nrow(ix)
  S = matrix(0, n, n)
  for (p in 1:n) for (q in 1:n) {
    i = ix[p, 1]; j = ix[p, 2]; kk = ix[q, 1]; l = ix[q, 2]
    S[p, q] = H[i, kk]*H[j, l] + H[i, l]*H[j, kk] + sum((kap - 3)*B[i, ]*B[j, ]*B[kk, ]*B[l, ])
  }
  S
}
# spot variance spillovers at every date (percent), sub-indices
fhspill = function(e, Hs, C, A, G, M, kap, chol) {
  Tn = nrow(e)
  k = ncol(e)
  n = k*(k + 1)/2
  D = dupl(k)
  Dp = solve(t(D) %*% D) %*% t(D)
  Aa = Dp %*% t(kronecker(A, A)) %*% D
  Bb = Dp %*% t(kronecker(G, G)) %*% D
  As = Aa + Bb
  mp = function(X, m) { R = diag(n); if (m > 0) for (i in 1:m) R = R %*% X; R }
  Th = lapply(0:(M - 1), function(m) if (m == 0) diag(n) else mp(As, m) - mp(As, m - 1) %*% Bb)
  isv = vechidx(k)[, 1] == vechidx(k)[, 2]
  CC = C %*% t(C)
  F = matrix(NA, Tn, n*n)
  S5 = matrix(NA, Tn, 5)
  for (t in 1:Tn) {
    Hf = CC + t(A) %*% (e[t, ] %o% e[t, ]) %*% A + t(G) %*% Hs[[t]] %*% G
    R = vector("list", M)
    for (h in 1:M) {
      if (h > 1) Hf = CC + t(A) %*% Hf %*% A + t(G) %*% Hf %*% G
      Sg = sigvech(Hf, kap)
      R[[h]] = if (chol) t(chol((Sg + t(Sg))/2)) else msqrt(Sg)
    }
    num = Reduce(`+`, lapply(0:(M - 1), function(m) (Th[[m + 1]] %*% R[[M - m]])^2))
    lam = num/rowSums(num)
    F[t, ] = 100*as.vector(t(lam))
    off = function(r, c) { X = lam[r, c, drop = FALSE]; if (all(r == c)) X = X - diag(diag(X), nrow(X)); sum(X) }
    S5[t, 1:4] = 100*c(off(isv, isv), off(!isv, !isv), sum(lam[isv, !isv]), sum(lam[!isv, isv]))/n
    S5[t, 5] = S5[t, 3] - S5[t, 4]
  }
  list(F = F, S5 = S5)
}
covnames = function(k) { o = c(); for (i in 1:k) for (j in i:k) o = c(o, paste0("cov_", i, "_", j)); o }

score_fun = list(bekk = BEKKs:::score_bekk, dbekk = BEKKs:::score_dbekk, sbekk = BEKKs:::score_sbekk)
ll_fun = list(bekk = BEKKs:::loglike_bekk, dbekk = BEKKs:::loglike_dbekk, sbekk = BEKKs:::loglike_sbekk)

# Per-observation scores of the Gaussian log likelihood, own code. H_1 = e'e/T does not
# depend on the parameters, so dvec(H_1) = 0. BEKKs' score functions start instead
# from dvec(H_1)/dtheta_C = d vec(CC')/dtheta_C (h1c = TRUE reproduces them); its
# scores are therefore not the exact gradient of its own likelihood, its BHHH steps
# and its standard errors differ slightly, and its iterations can stop short of the
# maximum. The exact scores are used for the standard errors written here and for the
# exact-gradient BHHH below.
# OPG standard errors; missing for a parameter whose scores vanish (a diagonal
# element of C at zero), the others from the rest of the matrix
opgse = function(S) {
  M = crossprod(S)
  kp = which(diag(M) > 1e-10*max(diag(M)))
  se = rep(NA, ncol(S))
  se[kp] = sqrt(diag(solve(M[kp, kp, drop = FALSE])))
  se
}
kmat = function(k) {
  K = matrix(0, k*k, k*k)
  for (i in 1:k) for (j in 1:k) K[(i - 1)*k + j, (j - 1)*k + i] = 1
  K
}
score_exact = function(th, e, type, h1c = FALSE) {
  TT = nrow(e); k = ncol(e); k2 = k*k; m = k*(k + 1)/2; np = length(th)
  p = unpack(th, k, type)
  C = p$C; A = p$A; G = p$G
  CC = C %*% t(C); GG = t(kronecker(G, G)); K = kmat(k); Ik = diag(k)
  ix = vechidx(k)
  SC = matrix(0, k2, m)
  for (c in 1:m) SC[(ix[c, 2] - 1)*k + ix[c, 1], c] = 1
  JC = (diag(k2) + K) %*% kronecker(C, Ik) %*% SC
  dg = (0:(k - 1))*(k + 1) + 1
  H = crossprod(e)/TT
  dH = matrix(0, k2, np)
  if (h1c) dH[, 1:m] = JC
  S = matrix(0, TT, np)
  for (t in 1:TT) {
    if (t > 1) {
      Hp = H
      X = e[t - 1, ] %o% e[t - 1, ]
      H = CC + t(A) %*% X %*% A + t(G) %*% Hp %*% G
      if (type == "sbekk") dH = cbind(JC, c(X), c(Hp)) + th[m + 2]*dH
      else {
        Ja = kronecker(Ik, t(A) %*% X) + kronecker(t(A) %*% X, Ik) %*% K
        Jg = kronecker(Ik, t(G) %*% Hp) + kronecker(t(G) %*% Hp, Ik) %*% K
        if (type == "dbekk") { Ja = Ja[, dg]; Jg = Jg[, dg] }
        dH = cbind(JC, Ja, Jg) + GG %*% dH
      }
    }
    Hi = solve((H + t(H))/2)
    u = Hi %*% e[t, ]
    S[t, ] = -0.5*(c(Hi) - c(u %*% t(u))) %*% dH
  }
  S
}
# sign normalization (the likelihood does not change): positive diagonal of C,
# A[1,1] > 0, G[1,1] > 0
normth = function(th, k, type) {
  p = unpack(th, k, type)
  C = p$C
  for (i in 1:k) if (C[i, i] < 0) C[, i] = -C[, i]
  m = k*(k + 1)/2
  o = th
  o[1:m] = C[lower.tri(C, diag = TRUE)]
  if (type == "bekk") {
    if (th[m + 1] < 0) o[(m + 1):(m + k*k)] = -th[(m + 1):(m + k*k)]
    if (th[m + k*k + 1] < 0) o[(m + k*k + 1):(m + 2*k*k)] = -th[(m + k*k + 1):(m + 2*k*k)]
  } else if (type == "dbekk") {
    if (th[m + 1] < 0) o[(m + 1):(m + k)] = -th[(m + 1):(m + k)]
    if (th[m + k + 1] < 0) o[(m + k + 1):(m + 2*k)] = -th[(m + k + 1):(m + 2*k)]
  }
  o
}
# log likelihood without sign restrictions: BEKKs' likelihood at the normalized
# parameters (BEKKs requires a positive diagonal of C, A[1,1] > 0 and G[1,1] > 0)
llfree = function(th, e, type) ll_fun[[type]](normth(th, ncol(e), type), e)
# The maximum of the likelihood, independently of the Stata code: (1) BHHH with the
# step grid and stopping rule of bhh_bekk() and the exact scores, no sign
# restrictions, at most 200 steps; (2) BFGS (optim()); (3) Newton steps with the Hessian from central
# differences of the exact gradient (absolute eigenvalues, step halving) until
# g'|H|^-1 g < 1e-10 with a negative definite Hessian
mle_exact = function(th, e, type, max_iter = 200, crit = 1e-14) {
  steps = c(9.9, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5, 0.25, 0.1, 0.01, 0.005, 0.001, 0.0005, 0.0001,
            0.00005, 0.00001, 0)
  llf = function(t) llfree(t, e, type)
  grf = function(t) colSums(score_exact(t, e, type))
  lik = llf(th)
  it1 = 0; conv = FALSE
  while (it1 < max_iter && !conv) {
    S = score_exact(th, e, type)
    M = crossprod(S)
    d = tryCatch(as.vector(solve(M, colSums(S))), error = function(z) {
      eg = eigen(M, symmetric = TRUE)
      iv = ifelse(eg$values > max(eg$values)*1e-14, 1/eg$values, 0)
      as.vector(eg$vectors %*% (iv*(t(eg$vectors) %*% colSums(S))))
    })
    lc = sapply(steps, function(s) if (s == 0) lik else llf(th + 0.1*s*d))
    bi = which.max(lc)
    it1 = it1 + 1
    if ((lc[bi] - lik)^2/abs(lik) < crit) conv = TRUE
    th = th + 0.1*steps[bi]*d
    lik = lc[bi]
  }
  ll1 = lik
  S = score_exact(th, e, type)
  sc = opgse(S)
  sc[is.na(sc)] = 1e-4
  o = optim(th, function(t) { v = llf(t); if (v <= -1e25) 1e25 else -v }, function(t) -grf(t),
            method = "BFGS", control = list(maxit = 2000, reltol = 1e-16, parscale = sc))
  if (-o$value > lik) { th = o$par; lik = -o$value }
  it2 = o$counts[2]
  conv = FALSE; it3 = 0; hmax = NA
  while (it3 < 30) {
    g = grf(th)
    np = length(th)
    Hs = matrix(0, np, np)
    for (j in 1:np) {
      h = 1e-5*max(abs(th[j]), 1e-2)
      tp = th; tp[j] = tp[j] + h; tm = th; tm[j] = tm[j] - h
      Hs[, j] = (grf(tp) - grf(tm))/(2*h)
    }
    Hs = (Hs + t(Hs))/2
    eg = eigen(Hs, symmetric = TRUE)
    hmax = max(eg$values)
    d = as.vector(eg$vectors %*% ((t(eg$vectors) %*% g)/abs(eg$values)))
    if (sum(g*d) < 1e-10 && hmax < 0) { conv = TRUE; break }
    s = 1
    while (llf(th + s*d) < lik && s > 1e-10) s = s/2
    if (llf(th + s*d) < lik) break
    th = th + s*d
    lik = llf(th)
    it3 = it3 + 1
  }
  list(theta = normth(th, ncol(e), type), ll = lik, ll1 = ll1, iter = c(it1, it2, it3), conv = conv, hmax = hmax,
       grad = sum(g*d))
}

bekkcase = function(nm, d, type, M, ekurt = FALSE, chol = FALSE, virfs = list(), mle = TRUE) {
  x = getx(d)
  e = sweep(x, 2, colMeans(x))
  k = ncol(e)
  n = k*(k + 1)/2
  spec = bekk_spec(model = list(type = type, asymmetric = FALSE), init_values = start(e, type))
  fit = bekk_fit(spec, e, max_iter = 1000, crit = 1e-14)
  th = as.vector(fit$theta)
  S = score_fun[[type]](matrix(th, ncol = 1), e)
  seb = sqrt(diag(solve(crossprod(S))))
  S1 = score_exact(th, e, type, h1c = TRUE)
  S2 = score_exact(th, e, type)
  se = opgse(S2)
  wr(rbind(th, se, c(fit$log_likelihood, fit$iter, rep(NA, length(th) - 2)), seb), paste0("R_bekk_", nm, "_theta.csv"))
  # the maximum of the likelihood (cases with mle = TRUE)
  if (mle) {
    mx = mle_exact(start(e, type), e, type)
    S3 = score_exact(mx$theta, e, type)
    se3 = opgse(S3)
    wr(rbind(mx$theta, se3, c(mx$ll, sum(mx$iter), mx$conv, mx$hmax, rep(NA, length(th) - 4))),
       paste0("R_bekk_", nm, "_mle.csv"))
  }
  cat(sprintf("%-14s own scores with BEKKs' H_1 convention = BEKKs' scores: max abs difference %.2e\n", nm,
              max(abs(S1 - S))))
  cat(sprintf("%-14s standard errors, BEKKs' scores against exact scores: max relative difference %.2e\n", nm,
              max(abs(seb - se)/se, na.rm = TRUE)))
  if (mle) cat(sprintf("%-14s maximum: ll %.6f (BEKKs %+.6f; after the BHHH stage %.6f), iterations %s, converged %s, largest Hessian eigenvalue %.3g, g'|H|^-1 g %.1e; max |theta - BEKKs| %.2e\n",
              nm, mx$ll, mx$ll - fit$log_likelihood, mx$ll1, paste(mx$iter, collapse = "/"), mx$conv, mx$hmax, mx$grad,
              max(abs(mx$theta - th))))
  Hm = as.matrix(fit$H_t)
  Hs = lapply(1:nrow(e), function(t) matrix(Hm[t, ], k, k))
  hu = t(sapply(Hs, function(H) { o = c(); for (i in 1:k) for (j in i:k) o = c(o, H[i, j]); o }))
  colnames(hu) = covnames(k)
  wcsv(thin(cbind(period = 1:nrow(e), hu)), paste0("R_bekk_", nm, "_H.csv"))
  p = unpack(th, k, type)
  if (ekurt) {
    xi = t(sapply(1:nrow(e), function(t) solve(msqrt(Hs[[t]]), e[t, ])))
    kap = colMeans(xi^4)
  } else kap = rep(3, k)
  wr(t(kap), paste0("R_bekk_", nm, "_kurt.csv"))
  fh = fhspill(e, Hs, p$C, p$A, p$G, M, kap, chol)
  ts = t(apply(fh$F, 1, function(r) permeas(matrix(r, n, n, byrow = TRUE)/100)))
  colnames(ts) = measnames(n)
  wcsv(thin(cbind(period = 1:nrow(e), ts)), paste0("R_", "bekk_", nm, "_ts.csv"))
  wr(matrix(colMeans(fh$F), n, n, byrow = TRUE), paste0("R_bekk_", nm, "_fevd.csv"))
  sub = fh$S5
  colnames(sub) = c("bekk_ovar", "bekk_ocov", "bekk_ccov", "bekk_cvar", "bekk_ncross")
  wcsv(thin(cbind(period = 1:nrow(e), sub)), paste0("R_bekk_", nm, "_sub.csv"))
  wr(t(colMeans(sub)), paste0("R_bekk_", nm, "_subavg.csv"))
  stopifnot(max(abs(rowSums(sub[, 1:4]) - ts[, "tci"])) < 1e-8)
  for (v in names(virfs)) {
    a = virfs[[v]]
    vr = virf(fit, time = a$time, q = a$q, index_series = a$index, n.ahead = a$n, time_shock = a$res)
    wr(as.matrix(vr$VIRF), paste0("R_bekk_", nm, "_virf_", v, ".csv"))
  }
  cat(sprintf("%-14s %-5s ll %.6f iter %d valid %s; mean TCI %.4f; sub-indices %s; kurtosis %s\n", nm, type,
              fit$log_likelihood, fit$iter, fit$BEKK_valid, mean(ts[, "tci"]),
              paste(sprintf("%.3f", colMeans(sub)), collapse = " "), paste(sprintf("%.3f", kap), collapse = " ")))
  invisible(list(e = e, th = th, Hs = Hs, p = p))
}

# Proposition 1 against a simulation: Cov(vech(e e')) at a fixed H, Gaussian and
# Student-t (unit variance, 8 degrees of freedom: kappa = 3 + 6/(8 - 4) = 4.5) shocks
set.seed(10)
Hx = matrix(c(1, 0.4, 0.2, 0.4, 2, 0.5, 0.2, 0.5, 1.5), 3, 3)
B = msqrt(Hx)
for (nu in c(Inf, 8)) {
  xi = if (is.infinite(nu)) matrix(rnorm(3e6), ncol = 3) else matrix(rt(3e6, nu)*sqrt((nu - 2)/nu), ncol = 3)
  ee = xi %*% B
  V = t(apply(ee, 1, function(r) { X = r %o% r; X[lower.tri(X, diag = TRUE)] }))
  kap = if (is.infinite(nu)) rep(3, 3) else rep(3 + 6/(nu - 4), 3)
  cat(sprintf("Proposition 1, %s shocks: max relative difference, simulation (1e6 draws) against formula %.4f\n",
              if (is.infinite(nu)) "Gaussian" else "t(8)", max(abs(cov(V) - sigvech(Hx, kap)))/max(abs(sigvech(Hx, kap)))))
}

# A. g2020 (EUR GBP CHF JPY), full BEKK, horizon 5 as Fengler and Herwartz, Gaussian
#    kurtosis, symmetric root; VIRFs at t = 2000 (EUR at its 5% quantile; the
#    residuals of the date). BEKKs' estimates only: for these returns (EUR and CHF
#    almost collinear under the EUR/CHF floor of 2011 to 2015) the likelihood of the
#    full BEKK rises towards the stationarity boundary and has no interior maximum
bekkcase("g2020_full", "g2020", "bekk", 5, mle = FALSE,
         virfs = list(q = list(time = 2000, q = 0.05, index = 1, n = 10, res = FALSE),
                      res = list(time = 2000, q = 0.05, index = 1, n = 10, res = TRUE)))
# B. g2020, diagonal BEKK, horizon 10, empirical kurtosis; VIRF of CHF at its 95% quantile
bekkcase("g2020_diag", "g2020", "dbekk", 10, ekurt = TRUE,
         virfs = list(q = list(time = 3000, q = 0.95, index = 3, n = 12, res = FALSE)))
# C. cgp2024 (GRE SOL WIN GEO), scalar BEKK, horizon 5, Cholesky identification
bekkcase("cgp2024_scalar", "cgp2024", "sbekk", 5, chol = TRUE,
         virfs = list(q = list(time = 1500, q = 0.05, index = 2, n = 10, res = FALSE)))
# D. cgp2024, full BEKK, horizon 10
bekkcase("cgp2024_full", "cgp2024", "bekk", 10)
