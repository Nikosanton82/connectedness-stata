"""
Reference port of the ConnectednessApproach R package (Gabauer, version 1.0.4, GPL-3),
used only to generate test fixtures for the Stata package `connectedness`.

Each function mirrors the R function of the same name line by line, including its
conventions (horizon h = 0..H-1, rolling window = window - nlag observations,
BayesPrior training sample = window - 2*nlag regression observations,
NPDC[i,j] = CT[i,j] - CT[j,i]). The definitive fixtures come from R itself
(test/r/make_fixtures.R); this port lets the Stata code be checked where R is absent.
"""
import numpy as np


# ----------------------------------------------------------------------------- helpers
def embed(x, d):
    """R's stats::embed: row t holds [x_t, x_{t-1}, ..., x_{t-d+1}]."""
    T, k = x.shape
    return np.hstack([x[d - 1 - j:T - j, :] for j in range(d)])


def ginv(a):
    """MASS::ginv: Moore-Penrose inverse with tolerance sqrt(eps) * max singular value."""
    u, s, vt = np.linalg.svd(a)
    tol = np.sqrt(np.finfo(float).eps) * s.max()
    sinv = np.where(s > tol, 1.0 / s, 0.0)
    return (vt.T * sinv) @ u.T


# ----------------------------------------------------------------------------- VAR (OLS)
def VAR(x, nlag=1):
    """Equation-by-equation OLS with intercept; returns B (k x k*nlag) and Q (k x k)."""
    k = x.shape[1]
    z = embed(x, nlag + 1)
    Z = z[:, k:]
    n = Z.shape[0]
    Zc = Z - Z.mean(0)
    XX = Zc.T @ Zc                      # cov(Z) * (n-1); the scale cancels in beta
    B, Res = [], []
    for i in range(k):
        y = z[:, i]
        yX = Zc.T @ (y - y.mean())
        beta = np.linalg.solve(XX, yX)
        resid = (y - y.mean()) - Z @ beta
        resid = resid - resid.mean()
        B.append(beta)
        Res.append(resid)
    Res = np.column_stack(Res)
    return np.array(B), Res.T @ Res / n


# ----------------------------------------------------------------------------- priors
def BayesPrior(x, size, nlag):
    p = x.shape[1]
    size = size - nlag
    yt = x[nlag:size + nlag, :].T
    m = p + nlag * p * p
    Zt = []
    for i in range(nlag, size + nlag):            # R: (nlag+1):(size+nlag), 1-based
        ztemp = [np.eye(p)]
        for j in range(1, nlag + 1):
            xlag = x[i - j, :]
            xtemp = np.zeros((p, p * p))
            for jj in range(p):
                xtemp[jj, jj * p:(jj + 1) * p] = xlag
            ztemp.append(xtemp)
        Zt.append(np.hstack(ztemp))
    vbar = np.zeros((m, m))
    xhy = np.zeros(m)
    for i in range(size):
        vbar += Zt[i].T @ Zt[i]
        xhy += Zt[i].T @ yt[:, i]
    vbar = ginv(vbar)
    aols = vbar @ xhy
    sse2 = np.zeros((p, p))
    for i in range(size):
        e = yt[:, i] - Zt[i] @ aols
        sse2 += np.outer(e, e)
    return dict(bprior=aols[p:], Vprior=vbar[p:, p:], Q=sse2 / size)


def MinnesotaPrior(gamma, k, nlag):
    m = nlag * k * k
    v = np.array([gamma / np.ceil((j + 1) / k) ** 2 for j in range(m // k) for i in range(k)])
    return dict(bprior=np.zeros(m), Vprior=np.diag(v), Q=None)


def UninformativePrior(k, nlag):
    m = nlag * k * k
    return dict(bprior=np.zeros(m), Vprior=4.0 * np.eye(m), Q=None)


# ----------------------------------------------------------------------------- TVP-VAR
def TVPVAR(x, kappa1=0.99, kappa2=0.99, nlag=1, prior=None):
    T, r = x.shape
    if prior is None:
        prior = UninformativePrior(r, nlag)
    Q = prior['Q'] if prior['Q'] is not None else np.cov(x, rowvar=False)
    x = x - x.mean(0)
    m = nlag * r * r
    kk = nlag * r
    yy = x[nlag:, :]
    xx = embed(x, nlag + 1)[:, r:]
    # Regressor rows: block for time irep (1-based) = rows (irep-1)*r .. irep*r-1;
    # the first nlag blocks are zero.
    Flag = np.zeros((T * r, m))
    for i in range(T - nlag):
        for j in range(nlag):
            xtemp = xx[i, j * r:(j + 1) * r]
            for q in range(r):
                Flag[(i + nlag) * r + q, j * r * r + q * r: j * r * r + (q + 1) * r] = xtemp
    beta_update = np.zeros((m, T))
    Sb = np.zeros((m, m, T))
    Q_t = np.zeros((r, r, T))
    beta_t = np.zeros((kk, kk, T))
    B = None
    for irep in range(1, T + 1):
        it = irep - 1
        if irep == 1:
            Q_t[:, :, it] = Q
        else:
            if irep <= nlag + 1:
                Gf = 0.1 * np.outer(x[it, :], x[it, :])
            else:
                e = yy[irep - nlag - 1, :] - xx[irep - nlag - 1, :] @ B[:r, :kk].T
                Gf = np.outer(e, e)
            Q_t[:, :, it] = kappa2 * Q_t[:, :, it - 1] + (1 - kappa2) * Gf
        if irep <= nlag + 1:
            beta_pred = prior['bprior'].copy()
            beta_update[:, it] = beta_pred
            Rb = prior['Vprior'].copy()
        else:
            beta_pred = beta_update[:, it - 1].copy()
            Rb = Sb[:, :, it - 1] / kappa1
        if irep >= nlag + 1:
            F = Flag[it * r:(it + 1) * r, :]
            Rx = Rb @ F.T
            KV = Q_t[:, :, it] + F @ Rx
            KG = Rx @ ginv(KV)
            beta_update[:, it] = beta_pred + KG @ (x[it, :] - F @ beta_pred)
            Sb[:, :, it] = Rb - KG @ (F @ Rb)
        bb = beta_update[:, it]
        biga = np.zeros((r, kk))
        s = 0
        for ii in range(nlag):
            for iii in range(r):
                biga[iii, ii * r:(ii + 1) * r] = bb[s:s + r]
                s += r
        if nlag > 1:
            B = np.vstack([biga, np.hstack([np.eye(r * (nlag - 1)), np.zeros((r * (nlag - 1), r))])])
        else:
            B = biga
        if np.max(np.abs(np.linalg.eigvals(B))) <= 1 or irep == 1:
            beta_t[:, :, it] = B
        else:
            beta_t[:, :, it] = beta_t[:, :, it - 1]
            beta_update[:, it] = 0.99 * beta_update[:, it - 1]
    return beta_t[:r, :, :], Q_t


# ----------------------------------------------------------------------------- IRF / FEVD
def PSI(Phi, nfore):
    """Moving-average matrices PSI_0..PSI_nfore (R function IRF with orth=FALSE)."""
    k = Phi.shape[0]
    p = Phi.shape[1] // k
    out = [np.eye(k)]
    for il in range(1, nfore + 1):
        tmp = np.zeros((k, k))
        for j in range(1, min(il, p) + 1):
            tmp += Phi[:, (j - 1) * k:j * k] @ out[il - j]
        out.append(tmp)
    return out


def FEVD(Phi, Sigma, nfore=10, generalized=True):
    psi = PSI(Phi, nfore)[:nfore]          # h = 0..nfore-1
    if generalized:
        sd = np.sqrt(np.diag(Sigma))
        num = sum(((P @ Sigma) / sd) ** 2 for P in psi)
        return num / num.sum(1, keepdims=True)
    L = np.linalg.cholesky(Sigma)          # lower; R uses t(chol(Sigma))
    denom = np.diag(sum(P @ Sigma @ P.T for P in psi))
    num = sum((P @ L) ** 2 for P in psi)
    return num / denom[:, None]


def ConnectednessTable(ct):
    """ct: k x k FEVD in percent (already averaged when used for the table)."""
    k = ct.shape[0]
    own = np.diag(np.diag(ct))
    TO = (ct - own).sum(0)
    FROM = ct.sum(1) - np.diag(ct)
    NET = TO - FROM
    TCI = FROM.mean()
    cTCI = TCI * k / (k - 1)
    NPDC = ct - ct.T
    NPT = (NPDC < 0).sum(1)
    with np.errstate(invalid='ignore', divide='ignore'):
        INFLUENCE = 100 * np.abs(NPDC / (ct.T + ct))
    PCI = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            PCI[i, j] = 200 * (ct[i, j] + ct[j, i]) / (ct[i, i] + ct[i, j] + ct[j, i] + ct[j, j])
    return dict(TO=TO, FROM=FROM, NET=NET, TCI=TCI, cTCI=cTCI, NPDC=NPDC, NPT=NPT,
                INFLUENCE=INFLUENCE, PCI=PCI, INCOWN=ct.sum(0))


def TimeConnectedness(B_t, Q_t, nfore=10, generalized=True, corrected=False):
    t = Q_t.shape[2]
    CT = np.stack([FEVD(B_t[:, :, i], Q_t[:, :, i], nfore, generalized) for i in range(t)], axis=2)
    res = [ConnectednessTable(100 * CT[:, :, i]) for i in range(t)]
    out = dict(CT=CT,
               TCI=np.array([r['cTCI'] if corrected else r['TCI'] for r in res]),
               TO=np.array([r['TO'] for r in res]), FROM=np.array([r['FROM'] for r in res]),
               NET=np.array([r['NET'] for r in res]), NPT=np.array([r['NPT'] for r in res]),
               NPDC=np.stack([r['NPDC'] for r in res], 2), PCI=np.stack([r['PCI'] for r in res], 2),
               INFLUENCE=np.stack([r['INFLUENCE'] for r in res], 2))
    out['TABLE'] = ConnectednessTable(100 * CT.mean(2))
    out['CTavg'] = 100 * CT.mean(2)
    return out


def ConnectednessApproach(x, nlag=1, nfore=10, window=None, corrected=False, model='VAR',
                          generalized=True, kappa1=0.99, kappa2=0.99, prior='BayesPrior', gamma=0.01):
    T, k = x.shape
    if window is None:
        ws, t0 = T, 1
    else:
        ws = window - nlag
        t0 = T - ws + 1
    if model == 'TVP-VAR':
        if prior == 'MinnesotaPrior':
            pr = MinnesotaPrior(gamma, k, nlag)
        elif prior == 'UninformativePrior':
            pr = UninformativePrior(k, nlag)
        else:
            pr = BayesPrior(x, ws, nlag)
        B_t, Q_t = TVPVAR(x, kappa1, kappa2, nlag, pr)
    else:
        B_t = np.zeros((k, k * nlag, t0))
        Q_t = np.zeros((k, k, t0))
        for i in range(t0):
            B_t[:, :, i], Q_t[:, :, i] = VAR(x[i:i + ws, :], nlag)
    out = TimeConnectedness(B_t, Q_t, nfore, generalized, corrected)
    out['B_t'], out['Q_t'] = B_t, Q_t
    return out


# ----------------------------------------------------------------------------- group measures
# Time-domain versions of the R functions AggregatedConnectedness, ConditionalConnectedness,
# InclusiveConnectedness, ExclusiveConnectedness, InternalConnectedness and
# ExternalConnectedness. CT is the k x k x t FEVD array in shares (rows sum to one),
# groups are lists of 0-based indices. Each returns the transformed array (shares) used
# for the per-period measures and the averaged table.
def group_transform(CT, kind, groups):
    k, _, t = CT.shape
    if kind == 'aggregated':
        m = len(groups)
        out = np.zeros((m, m, t))
        for i, gi in enumerate(groups):
            for j, gj in enumerate(groups):
                out[i, j, :] = CT[np.ix_(gi, gj)].sum(axis=(0, 1)) / len(gi)
        return out
    if kind == 'conditional':
        g = groups[0]
        sub = CT[np.ix_(g, g)]
        return sub / sub.sum(axis=1, keepdims=True)
    out = CT.copy()
    if kind == 'inclusive':
        keep = np.zeros((k, k), bool)
        for i in groups[0]:
            keep[:, i] = True
            keep[i, :] = True
        out[~keep] = 0
    elif kind == 'exclusive':
        for i in groups[0]:
            out[:, i, :] = 0
            out[i, :, :] = 0
    elif kind == 'internal':
        for a, ga in enumerate(groups):
            for b, gb in enumerate(groups):
                if a != b:
                    out[np.ix_(ga, gb)] = 0
    elif kind == 'external':
        for g in groups:
            out[np.ix_(g, g)] = 0
    return out


def GroupConnectedness(CT, kind, groups):
    """Per-period TCI/TO/FROM/NET/NPT/NPDC on the transformed array, the averaged table,
    and the group TCIs: TCI_ext (aggregated) or gTCI (internal, external)."""
    k = CT.shape[0]
    X = group_transform(CT, kind, groups)
    n, _, t = X.shape
    res = [ConnectednessTable(100 * X[:, :, i]) for i in range(t)]
    out = dict(X=X, TABLE=ConnectednessTable(100 * X.mean(2)), CTavg=100 * X.mean(2),
               TCI=np.array([r['TCI'] for r in res]), cTCI=np.array([r['cTCI'] for r in res]),
               TO=np.array([r['TO'] for r in res]), FROM=np.array([r['FROM'] for r in res]),
               NET=np.array([r['NET'] for r in res]), NPT=np.array([r['NPT'] for r in res]),
               NPDC=np.stack([r['NPDC'] for r in res], 2))
    if kind == 'aggregated':
        w = np.array([len(g) for g in groups])
        out['TCIext'] = (out['TO'] * (k - w) / k).sum(1)
        out['cTCIext'] = (out['TO'] * (k - w) / (k - 1)).sum(1)
    if kind in ('internal', 'external'):
        out['gTCI'] = np.column_stack([out['TO'][:, g].sum(1) / k for g in groups])
        out['cgTCI'] = np.column_stack([out['TO'][:, g].sum(1) / (k - 1) for g in groups])
    return out
