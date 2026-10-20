"""
Reference port (NumPy) of the ConnectednessApproach 1.0.4 extensions added in
connectedness 0.3.0 and later: frequency connectedness (FrequencyConnectedness, FEVD
type "frequency", frequencyConnectedness::getPartition), joint connectedness
(JointConnectedness), extended joint connectedness (ExtendedJointConnectedness), and
the quantile and LAD VARs (QVAR with quantreg's Frisch-Newton method; LADVAR).
It mirrors the Mata code in ado/cn_core.mata step by step and is used to check that
code against the R fixtures written by test/r/make_fixtures_ext.R.
"""
import numpy as np


def psi(Phi, H):
    """Moving-average matrices Psi_0..Psi_H (H+1 of them); Phi is k x k*p."""
    k = Phi.shape[0]
    p = Phi.shape[1] // k
    out = [np.eye(k)]
    for h in range(1, H + 1):
        P = np.zeros((k, k))
        for j in range(1, min(h, p) + 1):
            P += Phi[:, (j - 1) * k:j * k] @ out[h - j]
        out.append(P)
    return out


def get_indices(n, up, down):
    """frequencyConnectedness:::getIndeces, 0-based positions of the DFT."""
    m = n // 2
    space = np.arange(m + 1) / m * np.pi
    o = ((space >= down) & (space < up)).astype(int)
    if n % 2 == 0:
        o = np.concatenate([o, o[1:-1][::-1]])
    else:
        o = np.concatenate([o, o[1:][::-1]])
    return np.where(o == 1)[0]


def get_partition(partition, H):
    part = [get_indices(H + 1, partition[i], partition[i + 1]) for i in range(len(partition) - 1)]
    if any(len(q) == 0 for q in part):
        raise ValueError("empty frequency band")
    return part


def freq_fevd(Phi, Sigma, H, gen, rng):
    """List of n = H+1 k x k FEVD matrices, one per DFT frequency (R FEVD(type='frequency'))."""
    k = Sigma.shape[0]
    n = H + 1
    P = psi(Phi, H)                                   # H+1 matrices
    W = np.array([q.reshape(-1) for q in P])          # n x k^2, row h = vec (row-major) of Psi_h
    j = np.arange(n)
    Om = np.exp(-2j * np.pi * np.outer(j, j) / n)     # DFT matrix (R fft convention)
    F = Om @ W
    Ps = [F[r].reshape(k, k) for r in range(n)]
    den = np.zeros(k)
    for r in rng:
        den += np.real(np.diag(Ps[r] @ Sigma @ Ps[r].conj().T))
    den = den / H
    out = []
    if gen:
        sd2 = np.diag(Sigma)
        f = [np.abs(Ps[r] @ Sigma) ** 2 / n / np.outer(den, sd2) for r in range(n)]
        tot = sum(f[r] for r in rng).sum(axis=1)
        return [q / tot[:, None] for q in f]
    L = np.linalg.cholesky(Sigma)
    return [np.abs(Ps[r] @ L) ** 2 / n / den[:, None] for r in range(n)]


def measures(ct):
    k = ct.shape[0]
    to = ct.sum(0) - np.diag(ct)
    fr = ct.sum(1) - np.diag(ct)
    net = to - fr
    tci = fr.mean()
    npt = np.array([sum(ct[i, j] < ct[j, i] for j in range(k) if j != i) for i in range(k)])
    npdc, pci, inf = [], [], []
    for i in range(k - 1):
        for j in range(i + 1, k):
            npdc.append(ct[j, i] - ct[i, j])
            pci.append(200 * (ct[i, j] + ct[j, i]) / (ct[i, i] + ct[i, j] + ct[j, i] + ct[j, j]))
            inf.append(100 * abs(ct[j, i] - ct[i, j]) / (ct[i, j] + ct[j, i]))
    return dict(tci=tci, ctci=tci * k / (k - 1), to=to, fr=fr, net=net, npt=npt,
                npdc=np.array(npdc), pci=np.array(pci), inf=np.array(inf))


def freq_band_measures(ctb, ctt):
    """Band measures; PCI and influence use the total table in the denominators (R)."""
    m = measures(ctb)
    k = ctb.shape[0]
    pci, inf = [], []
    for i in range(k - 1):
        for j in range(i + 1, k):
            pci.append(200 * (ctb[i, j] + ctb[j, i]) / (ctt[i, i] + ctt[i, j] + ctt[j, i] + ctt[j, j]))
            inf.append(100 * abs(ctb[j, i] - ctb[i, j]) / (ctt[i, j] + ctt[j, i]))
    m['pci'] = np.array(pci)
    m['inf'] = np.array(inf)
    return m


def frequency_tables(Phi, Sigma, H, gen, partition, scenario="ABS"):
    """Band tables in percent (list) and their total."""
    k = Sigma.shape[0]
    part = get_partition(partition, H)
    rng = np.unique(np.concatenate(part))
    f = freq_fevd(Phi, Sigma, H, gen, rng)
    bands = []
    for q in part:
        t = sum(f[r] for r in q)
        if scenario == "WTH":
            t = t / t.sum() * k
        bands.append(100 * t)
    return bands, sum(bands)


def gfevd(Phi, Sigma, H):
    k = Sigma.shape[0]
    P = psi(Phi, H - 1)
    sd = np.sqrt(np.diag(Sigma))
    num = sum(((q @ Sigma) / sd) ** 2 for q in P)
    return num / num.sum(1, keepdims=True)


def joint_from(Phi, Sigma, H):
    k = Sigma.shape[0]
    P = psi(Phi, H - 1)
    Xi = sum(q @ Sigma @ q.T for q in P)
    s = np.zeros(k)
    for i in range(k):
        o = [j for j in range(k) if j != i]
        Si = np.linalg.inv(Sigma[np.ix_(o, o)])
        for q in P:
            a = q[i] @ Sigma[:, o]
            s[i] += a @ Si @ a
    return 100 * s / np.diag(Xi)


def joint(Phi, Sigma, H):
    g = 100 * gfevd(Phi, Sigma, H)
    off = g - np.diag(np.diag(g))
    gsoi = off.sum(1).mean()
    fr = joint_from(Phi, Sigma, H)
    jsoi = fr.mean()
    lam = gsoi / jsoi
    to = off.sum(0) / lam
    return dict(tci=jsoi, to=to, fr=fr, net=to - fr, ct=g)


def extjoint(Phi, Sigma, H):
    g = 100 * gfevd(Phi, Sigma, H)
    off = g - np.diag(np.diag(g))
    s = joint_from(Phi, Sigma, H)
    lam = s / off.sum(1)
    j = off * lam[:, None]
    np.fill_diagonal(j, 100 - j.sum(1))
    return j


# ----------------------------------------------------------------------------
# Quantile VAR and LAD VAR (0.4.0)
def rqfnb(X, yy, tau, beta=0.99995, eps=1e-6, maxit=500):
    n, p = X.shape
    A = X.T
    c = -yy
    b = (1 - tau) * X.sum(0)
    x = np.full(n, 1 - tau); u = np.ones(n); d = np.ones(n)
    def stepy(d, rhs):
        ada = (A * d) @ A.T
        L = np.linalg.cholesky(ada)
        return np.linalg.solve(L.T, np.linalg.solve(L, rhs)), L
    y, L = stepy(d, A @ c)
    s = c - A.T @ y
    small = np.abs(s) < eps
    z = np.maximum(s, 0) + small * eps
    w = np.maximum(-s, 0) + small * eps
    s = u - x
    gap = z @ x + w @ s
    nit = 0
    while gap > eps and nit < maxit:
        nit += 1
        d = 1 / (z / x + w / s)
        ds = z - w
        dz = d * ds
        dy = b - A @ x + A @ dz
        rhs = dy.copy()
        dy, L = stepy(d, dy)
        ds = A.T @ dy - ds
        dx = d * ds; ds = -dx
        dz = -z * (dx / x + 1); dw = -w * (ds / s + 1)
        def lim(v, dv):
            m = dv < 0
            return np.min(-v[m] / dv[m]) if m.any() else 1e20
        deltap = min(lim(x, dx), lim(s, ds)); deltad = min(lim(z, dz), lim(w, dw))
        deltap = min(beta * deltap, 1); deltad = min(beta * deltad, 1)
        if min(deltap, deltad) < 1:
            mu = x @ z + s @ w
            g = (mu + deltap * dx @ z + deltad * dz @ x + deltap * deltad * dz @ dx
                 + deltap * ds @ w + deltad * dw @ s + deltap * deltad * ds @ dw)
            mu = mu * (g / mu) ** 3 / (2 * n)
            dr = d * (mu * (1 / s - 1 / x) + dx * dz / x - ds * dw / s)
            dy = rhs + A @ dr
            dy = np.linalg.solve(L.T, np.linalg.solve(L, dy))
            u = A.T @ dy
            dxdz = dx * dz; dsdw = ds * dw
            dx = d * (u - z + w) - dr
            ds = -dx
            dz = -z + (mu - z * dx - dxdz) / x
            dw = -w + (mu - w * ds - dsdw) / s
            deltap = min(lim(x, dx), lim(s, ds)); deltad = min(lim(z, dz), lim(w, dw))
            deltap = min(beta * deltap, 1); deltad = min(beta * deltad, 1)
        x = x + deltap * dx; s = s + deltap * ds
        y = y + deltad * dy; z = z + deltad * dz; w = w + deltad * dw
        gap = z @ x + w @ s
    return -y


def qvar(x, nlag, tau, lad=False):
    """QVAR()/LADVAR(): B (k x k*nlag) and the uncentered residual covariance."""
    from ca_reference import embed
    k = x.shape[1]
    z = embed(x, nlag + 1)
    Y = z[:, :k]
    X = np.column_stack([np.ones(len(z)), z[:, k:]])
    B, R = [], []
    for i in range(k):
        b = rqfnb(X, Y[:, i], 0.5 if lad else tau[i])
        if lad:
            r = np.abs(Y[:, i] - X @ b)
            ix = np.argsort(r, kind="stable")[:X.shape[1]]
            b = np.linalg.solve(X[ix], Y[ix, i])
        B.append(b[1:])
        R.append(Y[:, i] - X @ b)
    R = np.array(R).T
    return np.array(B), R.T @ R / len(R)
