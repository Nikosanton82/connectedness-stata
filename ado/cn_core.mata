*! version 0.3.0  20oct2026
*! Mata engine for -connectedness-
*! Port of the ConnectednessApproach R package (D. Gabauer, version 1.0.4, GPL-3)
*! License: GPL-3

version 16.1

mata:
mata set matastrict on

// =============================================================================
// Priors for the TVP-VAR
// =============================================================================

struct cn_prior {
    real colvector b        // prior mean of the stacked coefficients (p*k*k)
    real matrix    V        // prior variance (p*k*k x p*k*k)
    real matrix    Q        // initial error covariance (k x k)
}

// R's stats::embed(): row t holds (x_t, x_{t-1}, ..., x_{t-d+1})
real matrix cn_embed(real matrix x, real scalar d)
{
    real scalar T, j
    real matrix out

    T   = rows(x)
    out = J(T-d+1, 0, .)
    for (j=0; j<d; j++) {
        out = out, x[|d-j, 1 \ T-j, cols(x)|]
    }
    return(out)
}

// Moore-Penrose inverse with MASS::ginv()'s tolerance, sqrt(eps)*max(s)
real matrix cn_ginv(real matrix A)
{
    real matrix    U, Vt
    real vector    s, sinv
    real scalar    tol

    svd(A, U, s, Vt)
    if (rows(s)==1) s = s'
    tol  = sqrt(epsilon(1)) * max(s)
    sinv = (s :> tol) :/ (s + (s :<= tol))
    return(Vt' * diag(sinv) * U')
}

// Kronecker block of regressors for one observation:
// F = (I_k # x_{t-1}, I_k # x_{t-2}, ..., I_k # x_{t-p})     (k x p*k*k)
real matrix cn_regblock(real matrix x, real scalar t, real scalar p)
{
    real scalar j, k
    real matrix F

    k = cols(x)
    F = J(k, 0, .)
    for (j=1; j<=p; j++) {
        F = F, (I(k) # x[t-j, .])
    }
    return(F)
}

// BayesPrior(): OLS VAR with intercept on the first size0-p observations
struct cn_prior scalar cn_bayesprior(real matrix x, real scalar size0, real scalar p)
{
    struct cn_prior scalar pr
    real scalar    k, m, size, i
    real matrix    vbar, sse, Z
    real colvector xhy, aols, e

    k    = cols(x)
    size = size0 - p
    m    = k + p*k*k
    vbar = J(m, m, 0)
    xhy  = J(m, 1, 0)
    for (i=p+1; i<=size+p; i++) {
        Z    = I(k), cn_regblock(x, i, p)
        vbar = vbar + cross(Z, Z)
        xhy  = xhy  + cross(Z, x[i, .]')
    }
    vbar = cn_ginv(vbar)
    aols = vbar * xhy
    sse  = J(k, k, 0)
    for (i=p+1; i<=size+p; i++) {
        Z   = I(k), cn_regblock(x, i, p)
        e   = x[i, .]' - Z*aols
        sse = sse + e*e'
    }
    pr.b = aols[|k+1 \ m|]
    pr.V = vbar[|k+1, k+1 \ m, m|]
    pr.Q = sse / size
    return(pr)
}

// MinnesotaPrior(): zero mean, variance gamma/lag^2
struct cn_prior scalar cn_minnesotaprior(real matrix x, real scalar gamma, real scalar p)
{
    struct cn_prior scalar pr
    real scalar    k, m, i, j
    real colvector v

    k = cols(x)
    m = p*k*k
    v = J(m, 1, .)
    for (j=1; j<=p*k; j++) {
        for (i=1; i<=k; i++) {
            v[(j-1)*k+i] = gamma / (ceil(j/k)^2)
        }
    }
    pr.b = J(m, 1, 0)
    pr.V = diag(v)
    pr.Q = variance(x)
    return(pr)
}

// UninformativePrior(): zero mean, variance 4
struct cn_prior scalar cn_uninfprior(real matrix x, real scalar p)
{
    struct cn_prior scalar pr
    real scalar m

    m    = p*cols(x)*cols(x)
    pr.b = J(m, 1, 0)
    pr.V = 4*I(m)
    pr.Q = variance(x)
    return(pr)
}

// =============================================================================
// Estimators
// =============================================================================

// VAR(p) by equation-by-equation OLS with intercept (R function VAR()).
// B: k x k*p, ordered (A_1, ..., A_p); Q: residual covariance divided by n
void cn_var(real matrix x, real scalar p, real matrix B, real matrix Q)
{
    real scalar k, n
    real matrix z, Y, Z, Zc, Yc, Beta, Res

    k    = cols(x)
    z    = cn_embed(x, p+1)
    Y    = z[|1, 1 \ rows(z), k|]
    Z    = z[|1, k+1 \ rows(z), cols(z)|]
    n    = rows(Z)
    Zc   = Z :- mean(Z)
    Yc   = Y :- mean(Y)
    Beta = lusolve(cross(Zc, Zc), cross(Zc, Yc))
    Res  = Yc - Z*Beta
    Res  = Res :- mean(Res)
    B    = Beta'
    Q    = cross(Res, Res) / n
}

// TVP-VAR with forgetting factors (R function TVPVAR(); Koop and Korobilis 2013).
// Bt: T x k*k*p, row t = rowshape(B_t, 1) with B_t = (A_1t, ..., A_pt)
// Qt: T x k*k,   row t = rowshape(Sigma_t, 1)
void cn_tvpvar(real matrix x0, real scalar kappa1, real scalar kappa2,
               real scalar p, struct cn_prior scalar pr,
               real matrix Bt, real matrix Qt)
{
    real scalar    T, r, m, kk, it, ii
    real matrix    x, xx, yy, Qcur, Qprev, Rb, Sb, Sbprev, F, Rx, KV, KG
    real matrix    B, Bkeep, biga, Gf
    real colvector bpred, bupd, bupdprev, e

    T  = rows(x0)
    r  = cols(x0)
    m  = p*r*r
    kk = p*r
    x  = x0 :- mean(x0)
    xx = cn_embed(x, p+1)
    xx = xx[|1, r+1 \ rows(xx), cols(xx)|]
    yy = x[|p+1, 1 \ T, r|]

    Bt = J(T, r*kk, .)
    Qt = J(T, r*r, .)
    B        = J(kk, kk, 0)
    Bkeep    = J(kk, kk, 0)
    Qprev    = pr.Q
    Sbprev   = J(m, m, 0)
    bupdprev = pr.b

    for (it=1; it<=T; it++) {
        // error covariance
        if (it==1) {
            Qcur = pr.Q
        }
        else {
            if (it <= p+1) {
                Gf = 0.1 * (x[it, .]' * x[it, .])
            }
            else {
                e  = yy[it-p, .]' - B[|1, 1 \ r, kk|] * xx[it-p, .]'
                Gf = e * e'
            }
            Qcur = kappa2*Qprev + (1-kappa2)*Gf
        }
        // coefficients: prediction step
        if (it <= p+1) {
            bpred = pr.b
            bupd  = bpred
            Rb    = pr.V
        }
        else {
            bpred = bupdprev
            bupd  = bpred
            Rb    = Sbprev / kappa1
        }
        Sb = J(m, m, 0)
        // update step
        if (it >= p+1) {
            F    = cn_regblock(x, it, p)
            Rx   = Rb * F'
            KV   = Qcur + F*Rx
            KG   = Rx * cn_ginv(KV)
            bupd = bpred + KG*(x[it, .]' - F*bpred)
            Sb   = Rb - KG*(F*Rb)
        }
        // coefficient matrix and companion form
        biga = J(r, 0, .)
        for (ii=1; ii<=p; ii++) {
            biga = biga, rowshape(bupd[|(ii-1)*r*r+1 \ ii*r*r|]', r)
        }
        if (p > 1) {
            B = biga \ (I(r*(p-1)), J(r*(p-1), r, 0))
        }
        else {
            B = biga
        }
        // keep only stable draws; otherwise shrink towards the previous state
        if (it==1) {
            Bkeep = B
        }
        else if (max(abs(eigenvalues(B))) <= 1) {
            Bkeep = B
        }
        else {
            bupd = 0.99*bupdprev
        }
        Bt[it, .] = rowshape(Bkeep[|1, 1 \ r, kk|], 1)
        Qt[it, .] = rowshape(Qcur, 1)
        Qprev    = Qcur
        Sbprev   = Sb
        bupdprev = bupd
    }
}

// =============================================================================
// Forecast error variance decompositions and connectedness measures
// =============================================================================

// FEVD at horizon H, summing h = 0, ..., H-1 (R function FEVD(), type "time").
// Row i, column j: share of variable i's forecast error variance due to shock j.
real matrix cn_fevd(real matrix B, real matrix S, real scalar H, real scalar gen)
{
    real scalar    k, p, h, j
    real matrix    PS, P, num, L
    real colvector den
    real rowvector sd

    k  = rows(B)
    p  = cols(B)/k
    PS = J(k, k*H, 0)
    PS[|1, 1 \ k, k|] = I(k)
    for (h=1; h<H; h++) {
        P = J(k, k, 0)
        for (j=1; j<=min((h, p)); j++) {
            P = P + B[|1, (j-1)*k+1 \ k, j*k|] * PS[|1, (h-j)*k+1 \ k, (h-j+1)*k|]
        }
        PS[|1, h*k+1 \ k, (h+1)*k|] = P
    }
    num = J(k, k, 0)
    if (gen) {
        sd = sqrt(diagonal(S))'
        for (h=0; h<H; h++) {
            P   = PS[|1, h*k+1 \ k, (h+1)*k|] * S
            num = num + (P :/ sd):^2
        }
        return(num :/ rowsum(num))
    }
    L = cholesky(S)
    if (hasmissing(L)) return(J(k, k, .))
    den = J(k, 1, 0)
    for (h=0; h<H; h++) {
        P   = PS[|1, h*k+1 \ k, (h+1)*k|]
        num = num + (P*L):^2
        den = den + diagonal(P*S*P')
    }
    return(num :/ den)
}

// Per-period measures from a k x k FEVD in percent. Order:
// TCI, cTCI, TO(k), FROM(k), NET(k), NPT(k), NPDC(pairs), PCI(pairs), INF(pairs)
// with pairs (1,2), (1,3), ..., (k-1,k); NPDC_ij = net transmission from i to j.
real rowvector cn_measures(real matrix ct)
{
    real scalar    k, i, j, c, np, tci
    real rowvector to, from, net, npt, npdc, pci, inf

    k    = rows(ct)
    np   = k*(k-1)/2
    to   = colsum(ct) - diagonal(ct)'
    from = rowsum(ct)' - diagonal(ct)'
    net  = to - from
    tci  = sum(from)/k
    npt  = J(1, k, 0)
    npdc = J(1, np, .)
    pci  = J(1, np, .)
    inf  = J(1, np, .)
    for (i=1; i<=k; i++) {
        for (j=1; j<=k; j++) {
            if (j != i) npt[i] = npt[i] + (ct[i, j] < ct[j, i])
        }
    }
    c = 0
    for (i=1; i<k; i++) {
        for (j=i+1; j<=k; j++) {
            c++
            npdc[c] = ct[j, i] - ct[i, j]
            pci[c]  = 200*(ct[i, j] + ct[j, i]) / (ct[i, i] + ct[i, j] + ct[j, i] + ct[j, j])
            inf[c]  = 100*abs(ct[j, i] - ct[i, j]) / (ct[i, j] + ct[j, i])
        }
    }
    return((tci, tci*k/(k-1), to, from, net, npt, npdc, pci, inf))
}

// Connectedness table: rows (variables, TO, Inc.Own, NET, NPT), columns (variables, FROM)
real matrix cn_table(real matrix ct, real scalar corrected)
{
    real scalar    k
    real rowvector m
    real matrix    tab

    k   = rows(ct)
    m   = cn_measures(ct)
    tab = J(k+4, k+1, .z)
    tab[|1, 1 \ k, k|]       = ct
    tab[|1, k+1 \ k, k+1|]   = m[|2+k+1 \ 2+2*k|]'                       // FROM
    tab[|k+1, 1 \ k+1, k|]   = m[|3 \ 2+k|]                             // TO
    tab[k+1, k+1]            = sum(m[|3 \ 2+k|])
    tab[|k+2, 1 \ k+2, k|]   = colsum(ct)                               // incl. own
    tab[|k+3, 1 \ k+3, k|]   = m[|2+2*k+1 \ 2+3*k|]                     // NET
    tab[k+3, k+1]            = (corrected ? m[2] : m[1])                // (c)TCI
    tab[|k+4, 1 \ k+4, k|]   = m[|2+3*k+1 \ 2+4*k|]                     // NPT
    return(tab)
}

// Names of the per-period measures, matching cn_measures()
string rowvector cn_measnames(real scalar k)
{
    real scalar       i, j
    string rowvector  nm, pr

    nm = ("tci", "ctci")
    for (i=1; i<=k; i++) nm = nm, "to_"   + strofreal(i)
    for (i=1; i<=k; i++) nm = nm, "from_" + strofreal(i)
    for (i=1; i<=k; i++) nm = nm, "net_"  + strofreal(i)
    for (i=1; i<=k; i++) nm = nm, "npt_"  + strofreal(i)
    pr = J(1, 0, "")
    for (i=1; i<k; i++) {
        for (j=i+1; j<=k; j++) pr = pr, strofreal(i) + "_" + strofreal(j)
    }
    nm = nm, ("npdc_" :+ pr), ("pci_" :+ pr), ("inf_" :+ pr)
    return(nm)
}

// Variable labels for the per-period measures
string rowvector cn_measlabels(string rowvector v)
{
    real scalar       i, j, k
    string rowvector  lb, pn, pp

    k  = cols(v)
    lb = ("Total connectedness index (TCI)", "Corrected TCI")
    for (i=1; i<=k; i++) lb = lb, "TO others: "   + v[i]
    for (i=1; i<=k; i++) lb = lb, "FROM others: " + v[i]
    for (i=1; i<=k; i++) lb = lb, "NET: "         + v[i]
    for (i=1; i<=k; i++) lb = lb, "Net pairwise transmitter count: " + v[i]
    pn = J(1, 0, "")
    pp = J(1, 0, "")
    for (i=1; i<k; i++) {
        for (j=i+1; j<=k; j++) {
            pn = pn, v[i] + " - " + v[j]
            pp = pp, v[i] + " , " + v[j]
        }
    }
    lb = lb, ("NPDC: " :+ pn), ("PCI: " :+ pp), ("Influence: " :+ pp)
    return(lb)
}

// =============================================================================
// Frequency connectedness (Barunik and Krehlik 2018): R functions
// FrequencyConnectedness() and FEVD(type="frequency"), with the frequency bands of
// frequencyConnectedness::getPartition()
// =============================================================================

// Moving-average matrices Psi_0, ..., Psi_H side by side (k x k*(H+1))
real matrix cn_psi(real matrix B, real scalar H)
{
    real scalar k, p, h, j
    real matrix PS, P

    k  = rows(B)
    p  = cols(B)/k
    PS = J(k, k*(H+1), 0)
    PS[|1, 1 \ k, k|] = I(k)
    for (h=1; h<=H; h++) {
        P = J(k, k, 0)
        for (j=1; j<=min((h, p)); j++) {
            P = P + B[|1, (j-1)*k+1 \ k, j*k|] * PS[|1, (h-j)*k+1 \ k, (h-j+1)*k|]
        }
        PS[|1, h*k+1 \ k, (h+1)*k|] = P
    }
    return(PS)
}

// getIndeces(n, up, down): DFT positions (1 to n) whose frequency lies in [down, up)
real colvector cn_getindices(real scalar n, real scalar up, real scalar down)
{
    real scalar    m, L
    real colvector sp, o

    m  = floor(n/2)
    sp = (0::m) / m * pi()
    o  = (sp :>= down) :& (sp :< up)
    L  = rows(o)
    if (mod(n, 2) == 0) o = o \ o[((L-1)::2)]
    else                o = o \ o[(L::2)]
    return(selectindex(o))
}

// n x nb indicator of the DFT positions in each band; part holds the band
// boundaries in radians in decreasing order (R's partition argument)
real matrix cn_partition(real rowvector part, real scalar H)
{
    real scalar    b, nb
    real matrix    P
    real colvector ix

    nb = cols(part) - 1
    P  = J(H+1, nb, 0)
    for (b=1; b<=nb; b++) {
        ix = cn_getindices(H+1, part[b], part[b+1])
        if (rows(ix) == 0) {
            errprintf("frequency band %g contains no frequency at horizon %g;\n", b, H)
            errprintf("increase horizon() or change the bands\n")
            exit(error(198))
        }
        P[ix, b] = J(rows(ix), 1, 1)
    }
    return(P)
}

// Band connectedness tables (percent), side by side: k x k*nb.
// P: indicator from cn_partition(); scen: "abs" or "wth"
real matrix cn_freqtables(real matrix B, real matrix S, real scalar H, real scalar gen,
                          real matrix P, string scalar scen)
{
    real scalar    k, n, nb, r, b, h
    real matrix    PS, W, Om, L, f, FQ, out, blk
    complex matrix F, Pw
    real colvector rng, den, tot
    real rowvector s2

    k   = rows(S)
    n   = H + 1
    nb  = cols(P)
    PS  = cn_psi(B, H)
    W   = J(n, k*k, .)
    for (h=0; h<=H; h++) W[h+1, .] = rowshape(PS[|1, h*k+1 \ k, (h+1)*k|], 1)
    Om  = (0::H) * (0..H)
    F   = exp(C(0, -2*pi()/n) :* Om) * W               // row r: DFT of Psi at frequency r
    rng = selectindex(rowsum(P) :> 0)
    den = J(k, 1, 0)
    for (r=1; r<=rows(rng); r++) {
        Pw  = rowshape(F[rng[r], .], k)
        den = den + Re(diagonal(Pw * S * Pw'))
    }
    den = den / H
    FQ  = J(n, k*k, .)                                  // row r: rowshape of the FEVD at r
    if (gen) {
        s2  = diagonal(S)'
        tot = J(k, 1, 0)
        for (r=1; r<=n; r++) {
            Pw = rowshape(F[r, .], k)
            f  = (abs(Pw * S):^2 / n) :/ (den * s2)
            FQ[r, .] = rowshape(f, 1)
            if (rowsum(P[r, .]) > 0) tot = tot + rowsum(f)
        }
        for (r=1; r<=n; r++) FQ[r, .] = rowshape(rowshape(FQ[r, .], k) :/ tot, 1)
    }
    else {
        L = cholesky(S)
        if (hasmissing(L)) return(J(k, k*nb, .))
        for (r=1; r<=n; r++) {
            Pw = rowshape(F[r, .], k)
            FQ[r, .] = rowshape((abs(Pw * L):^2 / n) :/ den, 1)
        }
    }
    out = J(k, k*nb, .)
    for (b=1; b<=nb; b++) {
        blk = rowshape(colsum(FQ[selectindex(P[., b]), .]), k)
        if (scen == "wth") blk = blk / sum(blk) * k
        out[|1, (b-1)*k+1 \ k, b*k|] = 100 * blk
    }
    return(out)
}

// Band measures as cn_measures(), except that the pairwise connectedness and
// influence indices use the total table in their denominators (as in R)
real rowvector cn_bandmeasures(real matrix cb, real matrix ct)
{
    real scalar    k, i, j, c, np
    real rowvector m, pci, inf

    k   = rows(cb)
    np  = k*(k-1)/2
    m   = cn_measures(cb)
    pci = J(1, np, .)
    inf = J(1, np, .)
    c   = 0
    for (i=1; i<k; i++) {
        for (j=i+1; j<=k; j++) {
            c++
            pci[c] = 200*(cb[i, j] + cb[j, i]) / (ct[i, i] + ct[i, j] + ct[j, i] + ct[j, j])
            inf[c] = 100*abs(cb[j, i] - cb[i, j]) / (ct[i, j] + ct[j, i])
        }
    }
    m[|2+4*k+np+1 \ 2+4*k+2*np|]   = pci
    m[|2+4*k+2*np+1 \ 2+4*k+3*np|] = inf
    return(m)
}

// =============================================================================
// Joint connectedness (Lastrapes and Wiesen 2021; R JointConnectedness()) and
// extended joint connectedness (Balcilar, Gabauer and Umar 2021; R
// ExtendedJointConnectedness())
// =============================================================================

// Joint spillovers FROM all others to each variable, in percent (row vector)
real rowvector cn_jointfrom(real matrix B, real matrix S, real scalar H)
{
    real scalar    k, i, h
    real matrix    PS, P, Xi, So, Gi
    real rowvector s, a
    real colvector o

    k  = rows(S)
    PS = cn_psi(B, H-1)
    Xi = J(k, k, 0)
    s  = J(1, k, 0)
    for (h=0; h<H; h++) {
        P  = PS[|1, h*k+1 \ k, (h+1)*k|]
        Xi = Xi + P*S*P'
    }
    for (i=1; i<=k; i++) {
        o  = selectindex((1::k) :!= i)
        So = S[., o]
        Gi = cn_ginv(S[o, o])
        for (h=0; h<H; h++) {
            a    = PS[|i, h*k+1 \ i, (h+1)*k|] * So
            s[i] = s[i] + a * Gi * a'
        }
    }
    return(100 * s :/ diagonal(Xi)')
}

// Joint connectedness of one period. ct: generalized FEVD (percent);
// returns (TCI, TO, FROM, NET); tb: this period's connectedness table as in R
real rowvector cn_joint(real matrix B, real matrix S, real scalar H,
                        real matrix ct, real matrix tb)
{
    real scalar    k, gsoi, jsoi
    real matrix    off
    real rowvector from, to

    k    = rows(S)
    ct   = 100 * cn_fevd(B, S, H, 1)
    off  = ct - diag(diagonal(ct))
    gsoi = mean(rowsum(off))
    from = cn_jointfrom(B, S, H)
    jsoi = mean(from')
    to   = colsum(off) / (gsoi/jsoi)
    tb   = J(k+3, k+1, 0)
    tb[|1, 1 \ k, k|]     = ct
    tb[|1, k+1 \ k, k+1|] = from'
    tb[|k+1, 1 \ k+1, k|] = to
    tb[k+1, k+1]          = sum(to)
    tb[|k+2, 1 \ k+2, k|] = colsum(ct)
    tb[|k+3, 1 \ k+3, k|] = to - from
    tb[k+3, k+1]          = jsoi
    return((jsoi, to, from, to - from))
}

// Extended joint connectedness table of one period (percent)
real matrix cn_extjoint(real matrix B, real matrix S, real scalar H)
{
    real matrix    g, off, jt
    real rowvector s

    g   = 100 * cn_fevd(B, S, H, 1)
    off = g - diag(diagonal(g))
    s   = cn_jointfrom(B, S, H)
    jt  = off :* (s' :/ rowsum(off))
    jt  = jt + diag(100 :- rowsum(jt))
    return(jt)
}

// Names and labels of the joint connectedness measures
string rowvector cn_jointnames(real scalar k)
{
    real scalar      i
    string rowvector nm

    nm = "tci"
    for (i=1; i<=k; i++) nm = nm, "to_"   + strofreal(i)
    for (i=1; i<=k; i++) nm = nm, "from_" + strofreal(i)
    for (i=1; i<=k; i++) nm = nm, "net_"  + strofreal(i)
    return(nm)
}

string rowvector cn_jointlabels(string rowvector v)
{
    real scalar      i, k
    string rowvector lb

    k  = cols(v)
    lb = "Joint total connectedness index (jSOI)"
    for (i=1; i<=k; i++) lb = lb, "Joint TO others: "   + v[i]
    for (i=1; i<=k; i++) lb = lb, "Joint FROM others: " + v[i]
    for (i=1; i<=k; i++) lb = lb, "Joint NET: "         + v[i]
    return(lb)
}

// =============================================================================
// Driver called by connectedness.ado
// =============================================================================

void cn_main(string scalar vars, string scalar touse, string scalar tvar,
             string scalar model, real scalar p, real scalar H, real scalar window,
             real scalar gen, real scalar corrected, string scalar prior,
             real scalar kappa1, real scalar kappa2, real scalar gamma,
             string scalar tabname, string scalar fevdname, string scalar scalname,
             string scalar names, string scalar dvar, string scalar approach,
             real rowvector part, string scalar scen, string scalar blabels)
{
    external real matrix   cn__M, cn__CT
    external real colvector cn__tv, cn__dv
    external string rowvector cn__mn, cn__ml, cn__fn, cn__fl
    struct cn_prior scalar pr
    real scalar    T, k, ws, t0, nT, t, nb, b, kk
    real matrix    x, Bt, Qt, B, Q, ct, avg, P, TB, tb, tabsum, blk
    real colvector tv
    real rowvector ma, m, cv, sc
    string rowvector nm, bl, fn, fl, mn0, ml0

    x  = st_data(., tokens(vars), touse)
    tv = st_data(., tvar, touse)
    T  = rows(x)
    k  = cols(x)
    kk = k*k
    nm = tokens(names)
    bl = tokens(blabels)

    if (window==0) {
        ws = T
        t0 = 1
    }
    else {
        ws = window - p
        t0 = T - ws + 1
    }

    if (model=="tvpvar") {
        if (prior=="minnesota")          pr = cn_minnesotaprior(x, gamma, p)
        else if (prior=="uninformative") pr = cn_uninfprior(x, p)
        else                              pr = cn_bayesprior(x, ws, p)
        cn_tvpvar(x, kappa1, kappa2, p, pr, Bt, Qt)
        nT = T
    }
    else {
        Bt = J(t0, k*k*p, .)
        Qt = J(t0, k*k, .)
        for (t=1; t<=t0; t++) {
            cn_var(x[|t, 1 \ t+ws-1, k|], p, B, Q)
            Bt[t, .] = rowshape(B, 1)
            Qt[t, .] = rowshape(Q, 1)
        }
        nT = t0
    }

    // names of the stored measures and of the decompositions
    nb = 0
    fn = J(1, 0, "")
    fl = J(1, 0, "")
    for (t=1; t<=k; t++) {
        for (b=1; b<=k; b++) {
            fn = fn, "fevd_" + strofreal(t) + "_" + strofreal(b)
            fl = fl, "FEVD share of " + nm[t] + " due to " + nm[b]
        }
    }
    if (approach=="joint") {
        cn__mn = cn_jointnames(k)
        cn__ml = cn_jointlabels(nm)
        fl     = "Generalized " :+ fl
    }
    else {
        cn__mn = cn_measnames(k)
        cn__ml = cn_measlabels(nm)
    }
    cn__fn = fn
    cn__fl = fl
    if (approach=="frequency") {
        P   = cn_partition(part, H)
        nb  = cols(P)
        mn0 = cn__mn
        ml0 = cn__ml
        for (b=1; b<=nb; b++) {
            cn__mn = cn__mn, mn0 :+ ("_f" + strofreal(b))
            cn__ml = cn__ml, ml0 :+ (" [" + bl[b] + "]")
            cn__fn = cn__fn, fn :+ ("_f" + strofreal(b))
            cn__fl = cn__fl, fl :+ (" [" + bl[b] + "]")
        }
    }

    cn__CT = J(nT, kk*(nb+1), .)
    cn__M  = J(nT, cols(cn__mn), .)
    tabsum = J(k+3, k+1, 0)
    for (t=1; t<=nT; t++) {
        B = rowshape(Bt[t, .], k)
        Q = rowshape(Qt[t, .], k)
        if (approach=="frequency") {
            TB = cn_freqtables(B, Q, H, gen, P, scen)
            ct = J(k, k, 0)
            for (b=1; b<=nb; b++) ct = ct + TB[|1, (b-1)*k+1 \ k, b*k|]
            m  = cn_measures(ct)
            cv = rowshape(ct, 1)
            for (b=1; b<=nb; b++) {
                blk = TB[|1, (b-1)*k+1 \ k, b*k|]
                m   = m, cn_bandmeasures(blk, ct)
                cv  = cv, rowshape(blk, 1)
            }
        }
        else if (approach=="joint") {
            m      = cn_joint(B, Q, H, ct, tb)
            tabsum = tabsum + tb
            cv     = rowshape(ct, 1)
        }
        else {
            if (approach=="extjoint") ct = cn_extjoint(B, Q, H)
            else                      ct = 100 * cn_fevd(B, Q, H, gen)
            m  = cn_measures(ct)
            cv = rowshape(ct, 1)
        }
        if (hasmissing(cv) | hasmissing(m)) {
            errprintf("connectedness measures could not be computed for period %g ", t)
            errprintf("(error covariance not positive definite)\n")
            exit(error(506))
        }
        cn__CT[t, .] = cv
        cn__M[t, .]  = m
    }
    cn__tv = tv[|T-nT+1 \ T|]
    if (dvar != "") {
        cn__dv = st_data(., dvar, touse)
        cn__dv = cn__dv[|T-nT+1 \ T|]
    }

    // time averages: tables and summary scalars
    avg = rowshape(mean(cn__CT[|1, 1 \ nT, kk|]), k)
    if (approach=="joint") {
        tb = tabsum / nT
        tb[k+2, k+1] = .z
        st_matrix(tabname, tb)
        st_matrixrowstripe(tabname, (J(k+3, 1, ""), (nm, ("TO", "Inc_Own", "NET"))'))
        ma = (tb[k+3, k+1], .)
    }
    else {
        cn_posttable(tabname, avg, corrected, nm)
        ma = cn_measures(avg)
    }
    st_matrixcolstripe(tabname, (J(k+1, 1, ""), (nm, "FROM")'))
    st_matrix(fevdname, avg)
    st_matrixrowstripe(fevdname, (J(k, 1, ""), nm'))
    st_matrixcolstripe(fevdname, (J(k, 1, ""), nm'))
    sc = (ma[1], ma[2], T, nT, tv[1], tv[T], cn__tv[1])
    for (b=1; b<=nb; b++) {
        avg = rowshape(mean(cn__CT[|1, b*kk+1 \ nT, (b+1)*kk|]), k)
        cn_posttable(tabname + "f" + strofreal(b), avg, corrected, nm)
        st_matrix(fevdname + "f" + strofreal(b), avg)
        st_matrixrowstripe(fevdname + "f" + strofreal(b), (J(k, 1, ""), nm'))
        st_matrixcolstripe(fevdname + "f" + strofreal(b), (J(k, 1, ""), nm'))
        ma = cn_measures(avg)
        sc = sc, ma[1], ma[2]
    }
    st_matrix(scalname, sc)
}

// Connectedness table of a time-averaged decomposition, with row and column names
void cn_posttable(string scalar tabname, real matrix avg, real scalar corrected,
                  string rowvector nm)
{
    real scalar k

    k = rows(avg)
    st_matrix(tabname, cn_table(avg, corrected))
    st_matrixrowstripe(tabname, (J(k+4, 1, ""), (nm, ("TO", "Inc_Own", "NET", "NPT"))'))
    st_matrixcolstripe(tabname, (J(k+1, 1, ""), (nm, "FROM")'))
}

// Write the per-period results into the current (results) frame
void cn_postframe(string scalar tvar, string scalar tfmt, real scalar savefevd,
                  string scalar dvar, string scalar dfmt)
{
    external real matrix    cn__M, cn__CT
    external real colvector cn__tv, cn__dv
    external string rowvector cn__mn, cn__ml, cn__fn, cn__fl
    real scalar      i
    real rowvector   idx

    st_addobs(rows(cn__M))
    idx = st_addvar("double", tvar)
    st_store(., idx, cn__tv)
    st_varformat(tvar, tfmt)
    if (dvar != "") {
        idx = st_addvar("double", dvar)
        st_store(., idx, cn__dv)
        st_varformat(dvar, dfmt)
    }
    idx = st_addvar("double", cn__mn)
    st_store(., idx, cn__M)
    for (i=1; i<=cols(cn__mn); i++) st_varlabel(cn__mn[i], substr(cn__ml[i], 1, 80))
    if (savefevd) {
        idx = st_addvar("double", cn__fn)
        st_store(., idx, cn__CT)
        for (i=1; i<=cols(cn__fn); i++) st_varlabel(cn__fn[i], substr(cn__fl[i], 1, 80))
    }
}


// =============================================================================
// Group measures (estat group): aggregated, conditional, inclusive, exclusive,
// internal and external connectedness (time domain), as in the R package
// =============================================================================

// Transform one period's FEVD (shares, k x k). gid: group id per variable (0 = none)
real matrix cn_gtransform(real matrix ct, string scalar type, real colvector gid,
                          real scalar m)
{
    real scalar    a, b
    real matrix    X, S, G1, G2, ing
    real colvector ia, ib

    if (type=="aggregated") {
        X = J(m, m, 0)
        for (a=1; a<=m; a++) {
            ia = selectindex(gid:==a)
            for (b=1; b<=m; b++) {
                ib = selectindex(gid:==b)
                X[a, b] = sum(ct[ia, ib]) / rows(ia)
            }
        }
        return(X)
    }
    if (type=="conditional") {
        ia = selectindex(gid:==1)
        S  = ct[ia, ia]
        return(S :/ rowsum(S))
    }
    // k x k masks: G1[i,j] = gid[i], G2[i,j] = gid[j] (Mata's colon operators do not
    // combine a column vector with a row vector, so the masks are built explicitly)
    G1  = gid * J(1, rows(gid), 1)
    G2  = G1'
    ing = (G1 :> 0) :| (G2 :> 0)
    if (type=="inclusive")  return(ct :* ing)
    if (type=="exclusive")  return(ct :* (1 :- ing))
    if (type=="internal")   return(ct :* (1 :- ((G1:>0) :& (G2:>0) :& (G1 :!= G2))))
    if (type=="external")   return(ct :* (1 :- ((G1:>0) :& (G1 :== G2))))
    _error(3498, "unknown group type")
    return(J(0, 0, .))
}

// Measures of a transformed table (percent): TCI, cTCI, TO, FROM, NET, NPT, NPDC
real rowvector cn_gmeasures(real matrix X)
{
    real scalar    n, np
    real rowvector m

    n  = rows(X)
    np = n*(n-1)/2
    m  = cn_measures(X)
    return(m[|1 \ 2+4*n+np|])
}

string rowvector cn_fevdnames(real scalar k)
{
    real scalar      i, j
    string rowvector fn

    fn = J(1, 0, "")
    for (i=1; i<=k; i++) {
        for (j=1; j<=k; j++) fn = fn, "fevd_" + strofreal(i) + "_" + strofreal(j)
    }
    return(fn)
}

void cn_group(string scalar type, real colvector gid, real scalar m,
              real scalar corrected, string scalar tabname, string scalar scalname,
              string scalar gname)
{
    external real matrix   cn__G
    real scalar    k, T, n, t, g, np, nx
    real matrix    F, X, avg, gt
    real rowvector meas, w, to
    real colvector ig

    k  = rows(gid)
    F  = st_data(., cn_fevdnames(k)) / 100
    T  = rows(F)
    n  = (type=="aggregated" ? m : (type=="conditional" ? sum(gid:==1) : k))
    np = n*(n-1)/2
    nx = (type=="aggregated" ? 2 : (type=="internal" | type=="external" ? 2*m : 0))
    cn__G = J(T, 2 + 4*n + np + nx, .)
    avg   = J(n, n, 0)
    w     = J(1, m, .)
    for (g=1; g<=m; g++) w[g] = sum(gid:==g)
    for (t=1; t<=T; t++) {
        X    = 100 * cn_gtransform(rowshape(F[t, .], k), type, gid, m)
        avg  = avg + X
        meas = cn_gmeasures(X)
        to   = meas[|3 \ 2+n|]
        if (type=="aggregated") {
            meas = meas, sum(to :* (k :- w))/k, sum(to :* (k :- w))/(k-1)
        }
        else if (type=="internal" | type=="external") {
            gt = J(1, 2*m, .)
            for (g=1; g<=m; g++) {
                ig = selectindex(gid:==g)
                gt[g]   = sum(to[ig]) / k
                gt[m+g] = sum(to[ig]) / (k-1)
            }
            meas = meas, gt
        }
        cn__G[t, .] = meas
    }
    avg  = avg / T
    meas = cn_gmeasures(avg)
    st_matrix(tabname, cn_table(avg, corrected))
    st_matrix(scalname, (meas[1], meas[2], n, T))
    if (nx > 0) st_matrix(gname, mean(cn__G[|1, 2+4*n+np+1 \ T, cols(cn__G)|]))
}

// Post group measures into the current frame (created empty by the caller)
void cn_postgroup(string scalar type, string scalar labnames, real colvector tv,
                  string scalar tvar, string scalar tfmt, real colvector dv,
                  string scalar dvar, string scalar dfmt, string scalar gnames)
{
    external real matrix cn__G
    real scalar      n, np, i, m
    real rowvector   idx
    string rowvector v, nm, lb, gn

    v  = tokens(labnames)
    n  = cols(v)
    np = n*(n-1)/2
    nm = cn_measnames(n)
    lb = cn_measlabels(v)
    nm = nm[|1 \ 2+4*n+np|]
    lb = lb[|1 \ 2+4*n+np|]
    if (type=="aggregated") {
        nm = nm, ("tciext", "ctciext")
        lb = lb, ("TCI_ext (ConnectednessApproach)", "Corrected TCI_ext (ConnectednessApproach)")
    }
    else if (type=="internal" | type=="external") {
        gn = tokens(gnames)
        m  = cols(gn)
        for (i=1; i<=m; i++) {
            nm = nm, "gtci_" + strofreal(i)
            lb = lb, "Group TCI: " + gn[i]
        }
        for (i=1; i<=m; i++) {
            nm = nm, "cgtci_" + strofreal(i)
            lb = lb, "Corrected group TCI: " + gn[i]
        }
    }
    st_addobs(rows(cn__G))
    idx = st_addvar("double", tvar)
    st_store(., idx, tv)
    st_varformat(tvar, tfmt)
    if (dvar != "") {
        idx = st_addvar("double", dvar)
        st_store(., idx, dv)
        st_varformat(dvar, dfmt)
    }
    idx = st_addvar("double", nm)
    st_store(., idx, cn__G)
    for (i=1; i<=cols(nm); i++) st_varlabel(nm[i], substr(lb[i], 1, 80))
}

// =============================================================================
// Network plot data (cnnetwork)
// =============================================================================

// Reads npdc_i_j or pci_i_j and net_i from the current frame (all rows, or row r)
// and returns node and edge matrices for drawing.
void cn_netdata(string scalar method, real scalar k, real scalar r,
                real scalar threshold, string scalar nodename, string scalar edgename,
                string scalar suf)
{
    real scalar    i, j, c, np, wmax, v, th
    real matrix    P, W, N, E
    real rowvector pv, net
    string rowvector pn, nn

    np = k*(k-1)/2
    pn = J(1, np, "")
    c  = 0
    for (i=1; i<k; i++) {
        for (j=i+1; j<=k; j++) {
            c++
            pn[c] = method + "_" + strofreal(i) + "_" + strofreal(j) + suf
        }
    }
    nn = J(1, k, "")
    for (i=1; i<=k; i++) nn[i] = "net_" + strofreal(i) + suf
    P   = st_data(., pn)
    net = st_data(., nn)
    if (r < .) {
        pv  = P[r, .]
        net = net[r, .]
    }
    else {
        pv  = mean(P)
        net = mean(net)
    }
    W = J(k, k, 0)
    c = 0
    for (i=1; i<k; i++) {
        for (j=i+1; j<=k; j++) {
            c++
            v = pv[c]
            if (method=="pci") {
                W[i, j] = v
                W[j, i] = v
            }
            else if (v > 0) W[i, j] = v
            else            W[j, i] = -v
        }
    }
    wmax = max(W)
    if (wmax > 0) W = W / wmax
    // nodes: x, y, label clock position, (unused), size weight, transmitter (1) or receiver (0)
    N = J(k, 6, .)
    for (i=1; i<=k; i++) {
        th = 2*pi()*(i-1)/k
        N[i, 1] = cos(th)
        N[i, 2] = sin(th)
        // label position (clock) outside the circle
        N[i, 3] = (cos(th) > 0.38 ? 3 : (cos(th) < -0.38 ? 9 : (sin(th) > 0 ? 12 : 6)))
        N[i, 4] = .
        N[i, 5] = (method=="pci" ? 1 : 1 + abs(net[i])/max((max(abs(net)), 1e-12)))
        N[i, 6] = (method=="pci" ? 1 : (net[i] >= 0))
    }
    // edges: x1, y1, x2, y2, weight, width class 1-5 (shortened to stop at nodes)
    E = J(0, 6, .)
    for (i=1; i<=k; i++) {
        for (j=1; j<=k; j++) {
            if (i==j | W[i, j] <= threshold | W[i, j] == 0) continue
            if (method=="pci" & j < i) continue
            E = E \ (N[i, 1] + 0.12*(N[j, 1]-N[i, 1]), N[i, 2] + 0.12*(N[j, 2]-N[i, 2]),
                     N[j, 1] - 0.12*(N[j, 1]-N[i, 1]), N[j, 2] - 0.12*(N[j, 2]-N[i, 2]),
                     W[i, j], ceil(5*W[i, j]))
        }
    }
    if (rows(E) == 0) E = J(1, 6, .)
    st_matrix(nodename, N)
    st_matrix(edgename, E)
}

end
