"""
Reference port (NumPy) of the ConnectednessApproach 1.0.4 extensions added in
connectedness 0.3.0: frequency connectedness (FrequencyConnectedness, FEVD type
"frequency", frequencyConnectedness::getPartition), joint connectedness
(JointConnectedness) and extended joint connectedness (ExtendedJointConnectedness).
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
