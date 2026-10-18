"""
Generate full-precision test fixtures for the Stata certification scripts.

Usage:  python3 make_fixtures.py   (run from test/python)
Output: test/fixtures/<case>_fevd.csv  averaged FEVD table in percent (k x k)
        test/fixtures/<case>_ts.csv    per-period measures (column `period` = row of the results
                                       frame; long series thinned to every 10th period)
The reference port reproduces the published R tables for dy2009, dy2012, gg2018
and cg2021 to every printed digit (see README in this folder).
"""
import os
import numpy as np
import pandas as pd
from ca_reference import ConnectednessApproach, ConnectednessTable, GroupConnectedness

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'data')
OUT = os.path.join(HERE, '..', 'fixtures')

CASES = {
    # name: (dataset, options)
    'dy2009_orth_static':  ('dy2009', dict(nlag=2, nfore=10, generalized=False)),
    'dy2009_orth_roll200': ('dy2009', dict(nlag=2, nfore=10, generalized=False, window=200)),
    'dy2012_gen_static':   ('dy2012', dict(nlag=4, nfore=10)),
    'dy2012_gen_roll200':  ('dy2012', dict(nlag=4, nfore=10, window=200)),
    'acg2020_tvp_bayes':   ('acg2020', dict(nlag=1, nfore=12, window=200, model='TVP-VAR',
                                            kappa1=0.99, kappa2=0.96)),
    'acg2020_tvp_minn':    ('acg2020', dict(nlag=1, nfore=12, model='TVP-VAR', prior='MinnesotaPrior',
                                            gamma=0.1, kappa1=0.99, kappa2=0.99)),
    'acg2020_tvp_unin':    ('acg2020', dict(nlag=1, nfore=12, model='TVP-VAR', prior='UninformativePrior')),
    'acg2020_tvp_lag2':    ('acg2020', dict(nlag=2, nfore=12, window=200, model='TVP-VAR')),
    'acg2020_var_roll100': ('acg2020', dict(nlag=1, nfore=12, window=100)),
    'gg2018_tvp_bayes':    ('gg2018', dict(nlag=1, nfore=10, window=200, model='TVP-VAR')),
    'cg2021_tvp_corr':     ('cg2021', dict(nlag=1, nfore=10, window=200, model='TVP-VAR', corrected=True)),
    'cgs2022_tvp_corr':    ('cgs2022', dict(nlag=1, nfore=10, window=250, model='TVP-VAR', corrected=True)),
    'cegg2022_tvp_minn':   ('cegg2022', dict(nlag=1, nfore=20, model='TVP-VAR', corrected=True,
                                             prior='MinnesotaPrior', gamma=0.1)),
}

# Group measures: (base case, type, groups as 0-based index lists)
GROUP_CASES = {
    'gg2018_internal':     ('gg2018_tvp_bayes', 'internal', [[0, 1, 2, 3], [4, 5, 6, 7]]),
    'gg2018_external':     ('gg2018_tvp_bayes', 'external', [[0, 1, 2, 3], [4, 5, 6, 7]]),
    'cgs2022_conditional': ('cgs2022_tvp_corr', 'conditional', [list(range(0, 18, 3))]),
    'cgs2022_aggregated':  ('cgs2022_tvp_corr', 'aggregated', [list(range(3*i, 3*i+3)) for i in range(6)]),
    'cegg2022_internal':   ('cegg2022_tvp_minn', 'internal', [[0, 1, 2], [3, 4, 5, 6], [7, 8, 9, 10]]),
    'cegg2022_external':   ('cegg2022_tvp_minn', 'external', [[0, 1, 2], [3, 4, 5, 6], [7, 8, 9, 10]]),
    'cegg2022_inclusive':  ('cegg2022_tvp_minn', 'inclusive', [[0, 1, 2]]),
    'cegg2022_exclusive':  ('cegg2022_tvp_minn', 'exclusive', [[0, 1, 2]]),
}

# Asymmetric connectedness (Adekoya et al. 2022): percentage changes of the levels in
# aaacgo2022, then all, positive and negative changes
ASYM_OPTS = dict(nlag=1, nfore=10, window=200, model='TVP-VAR', prior='MinnesotaPrior', gamma=0.1)


def load(name):
    d = pd.read_stata(os.path.join(DATA, f'cn_{name}.dta'))
    return d.drop(columns=[c for c in ('t', 'date', 'month') if c in d]).to_numpy(float)


def per_period(ct):
    """Measures in the column order used by the Stata results frame."""
    k = ct.shape[0]
    r = ConnectednessTable(ct)
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    npdc = [ct[j, i] - ct[i, j] for i, j in pairs]            # net from i to j
    pci = [r['PCI'][i, j] for i, j in pairs]
    inf = [r['INFLUENCE'][i, j] for i, j in pairs]
    return np.concatenate([[r['TCI'], r['cTCI']], r['TO'], r['FROM'], r['NET'], r['NPT'],
                           npdc, pci, inf])


def group_period(g, t, kind):
    n = g['TO'].shape[1]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    npdc = [g['NPDC'][j, i, t] for i, j in pairs]   # NPDC[i,j] = CT[i,j]-CT[j,i]; net i->j
    row = np.concatenate([[g['TCI'][t], g['cTCI'][t]], g['TO'][t], g['FROM'][t], g['NET'][t],
                          g['NPT'][t], npdc])
    if kind == 'aggregated':
        row = np.concatenate([row, [g['TCIext'][t], g['cTCIext'][t]]])
    if kind in ('internal', 'external'):
        row = np.concatenate([row, g['gTCI'][t], g['cgTCI'][t]])
    return row


def group_colnames(n, kind, m):
    pairs = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]
    c = (['tci', 'ctci'] + [f'to_{i}' for i in range(1, n + 1)] + [f'from_{i}' for i in range(1, n + 1)]
         + [f'net_{i}' for i in range(1, n + 1)] + [f'npt_{i}' for i in range(1, n + 1)]
         + [f'npdc_{i}_{j}' for i, j in pairs])
    if kind == 'aggregated':
        c += ['tciext', 'ctciext']
    if kind in ('internal', 'external'):
        c += [f'gtci_{g}' for g in range(1, m + 1)] + [f'cgtci_{g}' for g in range(1, m + 1)]
    return c


def write_case(case, CTavg, ts, cols):
    pd.DataFrame(CTavg).to_csv(os.path.join(OUT, f'{case}_fevd.csv'),
                               header=False, index=False, float_format='%.17g')
    df = pd.DataFrame(ts, columns=cols)
    df.insert(0, 'period', np.arange(1, len(df) + 1))
    step = 50 if len(df) > 2000 else 10     # keep fixtures small: every step-th period + last
    if len(df) > 300:
        df = df[(df['period'] % step == 1) | (df['period'] == len(df))]
    df.to_csv(os.path.join(OUT, f'{case}_ts.csv'), index=False, float_format='%.17g')


def colnames(k):
    pairs = [(i, j) for i in range(1, k + 1) for j in range(i + 1, k + 1)]
    return (['tci', 'ctci'] + [f'to_{i}' for i in range(1, k + 1)] + [f'from_{i}' for i in range(1, k + 1)]
            + [f'net_{i}' for i in range(1, k + 1)] + [f'npt_{i}' for i in range(1, k + 1)]
            + [f'npdc_{i}_{j}' for i, j in pairs] + [f'pci_{i}_{j}' for i, j in pairs]
            + [f'inf_{i}_{j}' for i, j in pairs])


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    results = {}
    for case, (ds, opts) in CASES.items():
        x = load(ds)
        k = x.shape[1]
        res = ConnectednessApproach(x, **opts)
        results[case] = res
        CT = res['CT']
        ts = np.array([per_period(100 * CT[:, :, t]) for t in range(CT.shape[2])])
        write_case(case, res['CTavg'], ts, colnames(k))
        print(f'{case:22s} k={k:2d} periods={CT.shape[2]:5d} TCI={res["TABLE"]["TCI"]:.4f}')
    for case, (base, kind, groups) in GROUP_CASES.items():
        g = GroupConnectedness(results[base]['CT'], kind, groups)
        n = g['TO'].shape[1]
        ts = np.array([group_period(g, t, kind) for t in range(g['TO'].shape[0])])
        write_case(case, g['CTavg'], ts, group_colnames(n, kind, len(groups)))
        print(f'{case:22s} n={n:2d} TCI={g["TABLE"]["TCI"]:.4f} cTCI={g["TABLE"]["cTCI"]:.4f}')
    lv = load('aaacgo2022')
    y = 100 * (lv[1:] - lv[:-1]) / lv[:-1]
    for part, z in (('all', y), ('pos', np.where(y < 0, 0, y)), ('neg', np.where(y > 0, 0, y))):
        res = ConnectednessApproach(z, **ASYM_OPTS)
        CT = res['CT']
        ts = np.array([per_period(100 * CT[:, :, t]) for t in range(CT.shape[2])])
        write_case(f'aaacgo2022_asym_{part}', res['CTavg'], ts, colnames(z.shape[1]))
        print(f'aaacgo2022_asym_{part:4s}  k={z.shape[1]:2d} TCI={res["TABLE"]["TCI"]:.4f}')
