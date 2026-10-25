* cert_connectedness.do
* Certification script for -connectedness-: compares Stata results with reference
* fixtures (Python port in test/python; R itself when test/r has been run).
* Run through run_all.do from the test folder.

version 16.1
set more off
set linesize 80
global cn_fails 0
global cn_tests 0

capture mata: mata drop cn_tsdiff()
mata:
real scalar cn_tsdiff(string scalar resframe, string scalar fxframe)
{
    string scalar    cur
    string rowvector v
    real colvector   per
    real matrix      A, B

    cur = st_framecurrent()
    st_framecurrent(fxframe)
    v   = st_varname((2..st_nvar()))
    per = st_data(., "period")
    A   = st_data(., v)
    st_framecurrent(resframe)
    B   = st_data(per, v)
    st_framecurrent(cur)
    return(mreldif(B, A))
}
end

capture program drop cn_check
program define cn_check
    syntax , case(string) tol(real) [prefix(string) MATrix(string) FRame(string) BANDs(integer 0) ///
        SUFfixes(string)]
    forvalues b = 1/`bands' {
        local suffixes `suffixes' f`b'
    }
    if "`matrix'" == "" local matrix e(fevd)
    if "`frame'" == ""  local frame `e(frame)'
    local f "fixtures/`prefix'`case'"
    capture confirm file "`f'_ts.csv"
    if _rc {
        di as txt "  `prefix'`case': no fixture, skipped"
        exit
    }
    tempname F E
    capture frame drop __fx
    frame create __fx
    local d1 = 0
    local ttxt "table "
    capture confirm file "`f'_fevd.csv"
    if !_rc {
        frame __fx {
            quietly import delimited using "`f'_fevd.csv", clear varnames(nonames) asdouble
            mkmat _all, matrix(`F')
        }
        matrix `E' = `matrix'
        local d1 = mreldif(`E', `F')
        foreach b of local suffixes {
            frame __fx {
                quietly import delimited using "`f'_fevd_`b'.csv", clear varnames(nonames) asdouble
                mkmat _all, matrix(`F')
            }
            matrix `E' = e(fevd_`b')
            local d1 = max(`d1', mreldif(`E', `F'))
        }
    }
    else local ttxt "(no table)"
    frame __fx: quietly import delimited using "`f'_ts.csv", clear varnames(1) asdouble
    mata: st_local("d2", strofreal(cn_tsdiff("`frame'", "__fx"), "%12.0g"))
    frame drop __fx
    local ok = (`d1' < `tol') & (`d2' < `tol')
    global cn_tests = $cn_tests + 1
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-26s "`prefix'`case'" " `ttxt'" as res %9.2e `d1' ///
        as txt "  series " as res %9.2e `d2' as txt "  " cond(`ok', "pass", "FAIL")
end

capture program drop cn_expect
program define cn_expect
    gettoken rc 0 : 0
    capture noisily `0'
    global cn_tests = $cn_tests + 1
    if _rc != `rc' {
        global cn_fails = $cn_fails + 1
        local got = _rc
        di as err `"  FAIL: expected r(`rc'), got r(`got'): `0'"'
    }
    else di as txt "  expected error r(`rc'): pass"
end

timer clear
local opt frame(cn_test) replace notable

di _n as txt "{hline 78}" _n "1. Static and rolling VAR (Diebold and Yilmaz 2009, 2012)" _n "{hline 78}"
use ../data/cn_dy2009, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
connectedness `v', lags(2) horizon(10) orthogonalized `opt'
foreach p in "" "R_" {
    cn_check, case(dy2009_orth_static) tol(1e-8) prefix(`p')
}
timer on 1
connectedness `v', lags(2) horizon(10) orthogonalized window(200) `opt'
timer off 1
foreach p in "" "R_" {
    cn_check, case(dy2009_orth_roll200) tol(1e-8) prefix(`p')
}

use ../data/cn_dy2012, clear
tsset t
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) `opt'
foreach p in "" "R_" {
    cn_check, case(dy2012_gen_static) tol(1e-8) prefix(`p')
}
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(10) window(200) `opt'
foreach p in "" "R_" {
    cn_check, case(dy2012_gen_roll200) tol(1e-8) prefix(`p')
}

di _n as txt "{hline 78}" _n "2. TVP-VAR (Antonakakis, Chatziantoniou and Gabauer 2020)" _n "{hline 78}"
use ../data/cn_acg2020, clear
tsset month
quietly ds t date month, not
local v `r(varlist)'
connectedness `v', model(tvpvar) horizon(12) window(200) kappa1(0.99) kappa2(0.96) `opt'
foreach p in "" "R_" {
    cn_check, case(acg2020_tvp_bayes) tol(1e-6) prefix(`p')
}
connectedness `v', model(tvpvar) horizon(12) prior(minnesota) gamma(0.1) `opt'
foreach p in "" "R_" {
    cn_check, case(acg2020_tvp_minn) tol(1e-6) prefix(`p')
}
connectedness `v', model(tvpvar) horizon(12) prior(uninformative) `opt'
foreach p in "" "R_" {
    cn_check, case(acg2020_tvp_unin) tol(1e-6) prefix(`p')
}
connectedness `v', model(tvpvar) lags(2) horizon(12) window(200) `opt'
foreach p in "" "R_" {
    cn_check, case(acg2020_tvp_lag2) tol(1e-6) prefix(`p')
}
connectedness `v', horizon(12) window(100) `opt'
foreach p in "" "R_" {
    cn_check, case(acg2020_var_roll100) tol(1e-8) prefix(`p')
}

di _n as txt "{hline 78}" _n "3. TVP-VAR replications (Gabauer and Gupta 2018; Chatziantoniou and Gabauer 2021)" _n "{hline 78}"
use ../data/cn_gg2018, clear
tsset month
quietly ds t date month, not
connectedness `r(varlist)', model(tvpvar) window(200) `opt'
foreach p in "" "R_" {
    cn_check, case(gg2018_tvp_bayes) tol(1e-6) prefix(`p')
}
use ../data/cn_cg2021, clear
tsset t
quietly ds t date, not
timer on 2
connectedness `r(varlist)', model(tvpvar) window(200) corrected `opt'
timer off 2
foreach p in "" "R_" {
    cn_check, case(cg2021_tvp_corr) tol(1e-6) prefix(`p')
}

di _n as txt "{hline 78}" _n "4. Stored results, postestimation and graphs" _n "{hline 78}"
use ../data/cn_dy2012, clear
tsset t
connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) window(200) datevar(date) savefevd ///
    frame(cn_test) replace
global cn_tests = $cn_tests + 1
capture noisily {
    assert "`e(cmd)'" == "connectedness"
    assert e(k) == 4 & e(lags) == 4 & e(horizon) == 10 & e(window) == 200
    assert e(T0) == _N - 200 + 4 + 1
    assert rowsof(e(table)) == 8 & colsof(e(table)) == 5
    assert abs(e(cTCI) - e(TCI)*4/3) < 1e-10
    frame cn_test {
        assert _N == e(T0)
        confirm variable t date tci ctci to_4 from_4 net_4 npt_4 npdc_3_4 pci_1_2 inf_2_4 fevd_4_4
        tempvar rs
        egen double `rs' = rowtotal(fevd_1_*)
        assert abs(`rs' - 100) < 1e-8
        assert abs(net_1 + net_2 + net_3 + net_4) < 1e-8
    }
    estat table, digits(3)
    connectedness
    cngraph tci, nodraw
    cngraph net, nodraw
    cngraph npdc, nodraw
    cngraph pci, select(SP500) nodraw
    estimates store roll200
    connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) model(tvpvar) window(200) ///
        frame(cn_tvp) replace notable
    cngraph tci, compare(roll200) nodraw
}
if _rc {
    global cn_fails = $cn_fails + 1
    di as err "  FAIL: stored results and postestimation (r(" _rc "))"
}
else di as txt "  stored results, estat, cngraph: pass"

di _n as txt "{hline 78}" _n "5. Group measures (Gabauer and Gupta 2018; Stenfors et al. 2022;" ///
    _n "   Chatziantoniou et al. 2023)" _n "{hline 78}"
use ../data/cn_gg2018, clear
tsset month
quietly ds t date month, not
connectedness `r(varlist)', model(tvpvar) window(200) `opt'
local G2 `""US: US_FPU US_MPU US_TPU US_CPU" "JP: JP_FPU JP_MPU JP_TPU JP_CPU""'
foreach ty in internal external {
    estat group, type(`ty') groups(`G2') frame(cn_grp) replace notable
    matrix GT = r(table)
    matrix GT = GT[1..8, 1..8]
    foreach p in "" "R_" {
        cn_check, case(gg2018_`ty') tol(1e-6) prefix(`p') matrix(GT) frame(cn_grp)
    }
}

use ../data/cn_cgs2022, clear
tsset t
quietly ds t date, not
timer on 3
connectedness `r(varlist)', model(tvpvar) window(250) corrected `opt'
timer off 3
foreach p in "" "R_" {
    cn_check, case(cgs2022_tvp_corr) tol(1e-6) prefix(`p')
}
estat group, type(conditional) groups(X2YCHF X2YDKK X2YEUR X2YGBP X2YNOK X2YSEK) ///
    frame(cn_grp) replace notable
matrix GT = r(table)
matrix GT = GT[1..6, 1..6]
foreach p in "" "R_" {
    cn_check, case(cgs2022_conditional) tol(1e-6) prefix(`p') matrix(GT) frame(cn_grp)
}
estat group, type(aggregated) frame(cn_grp) replace ///
    groups("CHF: X2YCHF X5YCHF X10YCHF" "DKK: X2YDKK X5YDKK X10YDKK" ///
    "EUR: X2YEUR X5YEUR X10YEUR" "GBP: X2YGBP X5YGBP X10YGBP" ///
    "NOK: X2YNOK X5YNOK X10YNOK" "SEK: X2YSEK X5YSEK X10YSEK")
matrix GT = r(table)
matrix GT = GT[1..6, 1..6]
foreach p in "" "R_" {
    cn_check, case(cgs2022_aggregated) tol(1e-6) prefix(`p') matrix(GT) frame(cn_grp)
}

use ../data/cn_cegg2022, clear
tsset t
quietly ds t date, not
connectedness `r(varlist)', model(tvpvar) horizon(20) prior(minnesota) gamma(0.1) corrected `opt'
foreach p in "" "R_" {
    cn_check, case(cegg2022_tvp_minn) tol(1e-6) prefix(`p')
}
local G3 `""Oil: Demand_Shock Supply_Shock Risk_Shock" "Exporters: Brazil Canada Norway Russia" "Importers: China Europe India Japan""'
foreach ty in internal external {
    estat group, type(`ty') groups(`G3') frame(cn_grp) replace notable
    matrix GT = r(table)
    matrix GT = GT[1..11, 1..11]
    foreach p in "" "R_" {
        cn_check, case(cegg2022_`ty') tol(1e-6) prefix(`p') matrix(GT) frame(cn_grp)
    }
}
foreach ty in inclusive exclusive {
    estat group, type(`ty') groups(Demand_Shock Supply_Shock Risk_Shock) frame(cn_grp) replace notable
    matrix GT = r(table)
    matrix GT = GT[1..11, 1..11]
    foreach p in "" "R_" {
        cn_check, case(cegg2022_`ty') tol(1e-6) prefix(`p') matrix(GT) frame(cn_grp)
    }
}

di _n as txt "{hline 78}" _n "6. Asymmetric connectedness (Adekoya et al. 2022)" _n "{hline 78}"
use ../data/cn_aaacgo2022, clear
tsset t
quietly ds t date, not
foreach v of varlist `r(varlist)' {
    quietly generate double r_`v' = 100*(`v' - L.`v')/L.`v'
}
connectedness r_*, model(tvpvar) window(200) prior(minnesota) gamma(0.1) asymmetric ///
    frame(cn_asy) replace notable
foreach p in "" "R_" {
    cn_check, case(aaacgo2022_asym_all) tol(1e-6) prefix(`p') matrix(e(fevd)) frame(cn_asy)
    cn_check, case(aaacgo2022_asym_pos) tol(1e-6) prefix(`p') matrix(e(fevd_pos)) frame(cn_asy_pos)
    cn_check, case(aaacgo2022_asym_neg) tol(1e-6) prefix(`p') matrix(e(fevd_neg)) frame(cn_asy_neg)
}

di _n as txt "{hline 78}" _n "7. New postestimation: estat group, cnnetwork, overlays, dates" _n "{hline 78}"
global cn_tests = $cn_tests + 1
capture noisily {
    connectedness, part(positive)
    cngraph tci, nodraw
    cngraph net, nodraw
    cngraph npdc, select(r_Brent) nodraw
    cnnetwork, nodraw
    cnnetwork, method(pci) threshold(0.2) nodraw
    use ../data/cn_gg2018, clear
    tsset month
    quietly ds t date month, not
    connectedness `r(varlist)', model(tvpvar) window(200) datevar(date) `opt'
    assert "`e(first)'" == "01feb1987"
    estat group, type(internal) groups(`G2') frame(cn_grp) replace
    assert rowsof(r(gtci)) == 1 & colsof(r(gtci)) == 4
    cngraph gtci, frame(cn_grp) nodraw
    cngraph net, frame(cn_grp) nodraw
    cnnetwork, frame(cn_grp) nodraw
    estat group, type(aggregated) groups(`G2') frame(cn_grp) replace notable
    cngraph tciext, frame(cn_grp) nodraw
    cnnetwork, at(tm(2008m9)) nodraw
}
if _rc {
    global cn_fails = $cn_fails + 1
    di as err "  FAIL: new postestimation features (r(" _rc "))"
}
else di as txt "  estat group, cnnetwork, cngraph overlays and dates: pass"

di _n as txt "{hline 78}" _n "8. Error handling" _n "{hline 78}"
use ../data/cn_dy2012, clear
tsset t
cn_expect 198 connectedness SP500 R_10Y, model(garch) `opt'
cn_expect 198 connectedness SP500 R_10Y, prior(minnesota) `opt'
cn_expect 198 connectedness SP500 R_10Y, model(tvpvar) kappa1(1) `opt'
cn_expect 184 connectedness SP500 R_10Y, orthogonalized generalized `opt'
cn_expect 2001 connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) window(20) `opt'
cn_expect 498 connectedness SP500 R_10Y in 1/100 if mod(t, 7) != 0, `opt'
cn_expect 110 connectedness SP500 R_10Y, frame(cn_test)
cn_expect 102 connectedness SP500, `opt'
quietly connectedness SP500 R_10Y DJUBSCOM USDX, `opt'
cn_expect 198 estat group, type(internal) groups("SP500 R_10Y DJUBSCOM USDX") frame(cn_grp) replace
cn_expect 198 estat group, type(internal) groups("SP500 R_10Y" "R_10Y USDX") frame(cn_grp) replace
cn_expect 198 estat group, type(internal) groups("SP500 NOTAVAR" "USDX") frame(cn_grp) replace
cn_expect 198 estat group, type(bogus) groups("SP500" "USDX") frame(cn_grp) replace
quietly connectedness SP500 R_10Y DJUBSCOM USDX, nofevd `opt'
cn_expect 198 estat group, type(external) groups("SP500 R_10Y" "USDX") frame(cn_grp) replace
cn_expect 198 connectedness, part(positive)

di _n as txt "{hline 78}" _n "9. Frequency, joint and extended joint connectedness (0.3.0)" _n "{hline 78}"
* fixtures from test/r/make_fixtures_ext.R (R package only)
use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
connectedness `v', lags(4) horizon(100) approach(frequency) bands(5 20) `opt'
cn_check, case(dy2012_freq_static) tol(1e-8) prefix(R_) bands(3)
connectedness `v', lags(4) horizon(10) approach(frequency) `opt'
cn_check, case(dy2012_freq_default) tol(1e-8) prefix(R_) bands(2)
connectedness `v', lags(4) horizon(100) approach(frequency) bands(5 20) scenario(wth) `opt'
cn_check, case(dy2012_freq_wth) tol(1e-8) prefix(R_) bands(3)
timer on 4
connectedness `v', lags(4) horizon(50) window(200) approach(frequency) bands(5) `opt'
timer off 4
cn_check, case(dy2012_freq_roll200) tol(1e-8) prefix(R_) bands(2)
connectedness `v', lags(4) horizon(10) window(200) approach(extjoint) `opt'
cn_check, case(dy2012_extjoint_roll200) tol(1e-8) prefix(R_)

use ../data/cn_dy2009, clear
tsset t
quietly ds t date, not
connectedness `r(varlist)', lags(2) horizon(100) orthogonalized approach(frequency) bands(5) `opt'
cn_check, case(dy2009_freq_orth) tol(1e-8) prefix(R_) bands(2)

use ../data/cn_cgg2022, clear
tsset t
quietly ds t date, not
timer on 5
connectedness `r(varlist)', model(tvpvar) horizon(100) window(200) approach(frequency) bands(5) `opt'
timer off 5
cn_check, case(cgg2022_freq_tvp) tol(1e-6) prefix(R_) bands(2)

use ../data/cn_acg2020, clear
tsset month
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) corrected ///
    approach(frequency) bands(3 12) `opt'
cn_check, case(acg2020_freq_tvp_corr) tol(1e-6) prefix(R_) bands(3)

use ../data/cn_lw2021, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
connectedness `v', lags(2) horizon(30) approach(joint) `opt'
cn_check, case(lw2021_joint_static) tol(1e-8) prefix(R_)
connectedness `v', lags(2) horizon(30) window(200) approach(joint) `opt'
cn_check, case(lw2021_joint_roll200) tol(1e-8) prefix(R_)
connectedness `v', lags(2) horizon(30) approach(extjoint) `opt'
cn_check, case(lw2021_extjoint_static) tol(1e-8) prefix(R_)

use ../data/cn_bgu2021, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
connectedness `v', model(tvpvar) horizon(20) prior(minnesota) gamma(0.1) approach(joint) `opt'
cn_check, case(bgu2021_joint_tvp) tol(1e-6) prefix(R_)
timer on 6
connectedness `v', model(tvpvar) horizon(20) prior(minnesota) gamma(0.1) approach(extjoint) `opt'
timer off 6
cn_check, case(bgu2021_extjoint_tvp) tol(1e-6) prefix(R_)

global cn_tests = $cn_tests + 1
capture noisily {
    use ../data/cn_dy2012, clear
    tsset t
    connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) horizon(50) window(200) datevar(date) ///
        approach(frequency) bands(5) `opt'
    assert e(nbands) == 2 & "`e(bands)'" == "1-5 5-Inf" & "`e(approach)'" == "frequency"
    assert reldif(e(TCI), e(TCI_f1) + e(TCI_f2)) < 1e-10
    connectedness, band(2)
    cngraph tci, nodraw
    cngraph net, band(1) nodraw
    cngraph npdc, band(2) select(SP500) nodraw
    cnnetwork, band(1) nodraw
    connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) approach(joint) `opt'
    assert e(cTCI) == . & rowsof(e(table)) == 7
    cngraph tci, nodraw
    cngraph net, nodraw
    connectedness SP500 R_10Y DJUBSCOM USDX, lags(4) approach(extjoint) window(200) `opt'
    cnnetwork, method(pci) nodraw
}
if _rc {
    global cn_fails = $cn_fails + 1
    di as err "  FAIL: frequency and joint postestimation (r(" _rc "))"
}
else di as txt "  frequency bands, joint tables, cngraph and cnnetwork band(): pass"
cn_expect 198 connectedness SP500 R_10Y, approach(joint) orthogonalized `opt'
cn_expect 198 connectedness SP500 R_10Y, approach(extjoint) corrected `opt'
cn_expect 198 connectedness SP500 R_10Y, approach(bogus) `opt'
cn_expect 198 connectedness SP500 R_10Y, bands(5) `opt'
cn_expect 198 connectedness SP500 R_10Y, approach(frequency) horizon(1) `opt'
cn_expect 198 connectedness SP500 R_10Y, approach(frequency) partition(1 2 0) `opt'
cn_expect 184 connectedness SP500 R_10Y, approach(frequency) bands(5) partition(3 1 0) `opt'
quietly connectedness SP500 R_10Y DJUBSCOM USDX, approach(frequency) `opt'
cn_expect 198 estat group, type(internal) groups("SP500 R_10Y" "USDX") frame(cn_grp) replace
cn_expect 198 connectedness, band(3)
quietly connectedness SP500 R_10Y DJUBSCOM USDX, approach(joint) `opt'
cn_expect 198 cngraph npdc, nodraw
cn_expect 198 cnnetwork, nodraw
cn_expect 198 connectedness, band(1)

di _n as txt "{hline 78}" _n "10. Quantile and LAD VAR connectedness (0.4.0)" _n "{hline 78}"
use ../data/cn_cgs2021, clear
tsset t
local v USD1Y EUR1Y GBP1Y JPY1Y
timer on 7
connectedness `v' in 1/700, model(qvar) horizon(20) window(200) corrected `opt'
timer off 7
cn_check, case(cgs2021_qvar_roll200) tol(1e-6) prefix(R_)
assert e(quantile) == 0.5
connectedness `v', model(qvar) quantile(0.05) horizon(20) corrected `opt'
cn_check, case(cgs2021_qvar_q05) tol(1e-6) prefix(R_)
connectedness `v', model(qvar) quantile(0.2 0.8 0.5 0.95) lags(2) horizon(20) `opt'
cn_check, case(cgs2021_qvar_qoq) tol(1e-6) prefix(R_)

use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
connectedness `v', model(lad) lags(4) `opt'
cn_check, case(dy2012_lad_static) tol(1e-6) prefix(R_)
connectedness `v' in 1/600, model(lad) lags(4) window(200) `opt'
cn_check, case(dy2012_lad_roll200) tol(1e-6) prefix(R_)

use ../data/cn_gcat2022, clear
tsset t
quietly ds t date, not
connectedness `r(varlist)' in 1/600, model(qvar) horizon(100) window(200) approach(frequency) ///
    bands(5) `opt'
cn_check, case(gcat2022_qvar_freq_roll200) tol(1e-6) prefix(R_) bands(2)

use ../data/cn_jcggh2022, clear
tsset t
quietly ds t date, not
connectedness `r(varlist)', model(qvar) quantile(0.9) approach(extjoint) `opt'
cn_check, case(jcggh2022_qvar_extjoint_q90) tol(1e-6) prefix(R_)

use ../data/cn_dy2012, clear
tsset t
cn_expect 198 connectedness SP500 R_10Y, model(qvar) quantile(0.2 0.5 0.8) `opt'
cn_expect 125 connectedness SP500 R_10Y, model(qvar) quantile(1) `opt'
cn_expect 198 connectedness SP500 R_10Y, quantile(0.5) `opt'
cn_expect 198 connectedness SP500 R_10Y, model(lad) prior(minnesota) `opt'
cn_expect 198 connectedness SP500 R_10Y, model(qvar) kappa2(0.96) `opt'

di _n as txt "{hline 78}" _n "11. R2 decomposed and model-free connectedness (0.5.0)" _n "{hline 78}"
use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
connectedness `v', approach(r2) lags(1) window(200) `opt'
cn_check, case(dy2012_r2_roll200) tol(1e-8) prefix(R_) suffixes(c l)
connectedness `v', approach(r2) lags(1) window(200) cormethod(spearman) `opt'
cn_check, case(dy2012_r2_spearman_roll200) tol(1e-8) prefix(R_) suffixes(c l)
timer on 8
connectedness `v', approach(r2) lags(1) window(200) cormethod(kendall) `opt'
timer off 8
cn_check, case(dy2012_r2_kendall_roll200) tol(1e-8) prefix(R_) suffixes(c l)
connectedness `v', approach(r2) lags(0) window(200) `opt'
cn_check, case(dy2012_r2_lag0_roll200) tol(1e-8) prefix(R_)
connectedness `v', approach(r2) lags(2) corrected `opt'
cn_check, case(dy2012_r2_lag2_static_corr) tol(1e-8) prefix(R_) suffixes(c l)
connectedness `v', approach(modelfree) `opt'
cn_check, case(dy2012_mf_static) tol(1e-8) prefix(R_)
connectedness `v', approach(modelfree) window(200) `opt'
cn_check, case(dy2012_mf_roll200) tol(1e-8) prefix(R_)
connectedness `v', approach(modelfree) window(200) cormethod(kendall) `opt'
cn_check, case(dy2012_mf_kendall_roll200) tol(1e-8) prefix(R_)
connectedness `v', approach(modelfree) cormethod(spearman) `opt'
cn_check, case(dy2012_mf_spearman_static) tol(1e-8) prefix(R_)

use ../data/cn_cgs2022, clear
tsset t
connectedness X2YCHF-X10YSEK, approach(r2) lags(1) `opt'
cn_check, case(cgs2022_r2_static) tol(1e-8) prefix(R_) suffixes(c l)

global cn_tests = $cn_tests + 1
capture noisily {
    use ../data/cn_dy2012, clear
    tsset t
    connectedness SP500 R_10Y DJUBSCOM USDX, approach(r2) lags(1) window(200) datevar(date) `opt'
    assert "`e(cormethod)'" == "pearson" & e(TCI_c) < . & e(TCI_l) < .
    assert reldif(e(TCI), e(TCI_c) + e(TCI_l)) < 1e-10
    connectedness, component(lagged)
    cngraph tci, nodraw
    cngraph net, component(c) nodraw
    cnnetwork, component(l) nodraw
    connectedness SP500 R_10Y DJUBSCOM USDX, approach(modelfree) window(200) `opt'
    estat group, type(inclusive) groups(SP500 R_10Y) frame(cn_grp) replace notable
}
if _rc {
    global cn_fails = $cn_fails + 1
    di as err "  FAIL: R2 and model-free postestimation (r(" _rc "))"
}
else di as txt "  R2 components, model-free group measures, cngraph and cnnetwork: pass"
cn_expect 198 connectedness SP500 R_10Y, approach(r2) model(tvpvar) `opt'
cn_expect 198 connectedness SP500 R_10Y, approach(r2) cormethod(bogus) `opt'
cn_expect 198 connectedness SP500 R_10Y, cormethod(kendall) `opt'
cn_expect 198 connectedness SP500 R_10Y, lags(0) `opt'
cn_expect 2001 connectedness SP500 R_10Y DJUBSCOM USDX, approach(r2) lags(2) window(12) `opt'
quietly connectedness SP500 R_10Y DJUBSCOM USDX, approach(r2) lags(0) `opt'
cn_expect 198 connectedness, component(lagged)

di _n as txt "{hline 78}" _n "12. LASSO, ridge and elastic-net VARs (0.6.0)" _n "{hline 78}"
* folds(interleaved): the folds used by test/r/make_fixtures_ext.R
use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
local f folds(interleaved)
connectedness `v', model(lasso) lags(2) `f' `opt'
cn_check, case(dy2012_lasso_static) tol(1e-8) prefix(R_)
connectedness `v', model(ridge) lags(2) `f' `opt'
cn_check, case(dy2012_ridge_static) tol(1e-8) prefix(R_)
timer on 9
connectedness `v', model(elasticnet) lags(2) `f' `opt'
timer off 9
cn_check, case(dy2012_elastic_static) tol(1e-8) prefix(R_)
connectedness `v', model(elasticnet) alpha(0.5) loss(mse) lags(1) `f' `opt'
cn_check, case(dy2012_elastic_a05_mse) tol(1e-8) prefix(R_)
connectedness `v' in 1/400, model(lasso) lags(1) window(200) `f' `opt'
cn_check, case(dy2012_lasso_roll200) tol(1e-8) prefix(R_)
connectedness `v' in 1/300, model(elasticnet) nalpha(4) lags(1) window(200) `f' `opt'
cn_check, case(dy2012_elastic4_roll200) tol(1e-8) prefix(R_)

use ../data/cn_gghm2022, clear
tsset t
quietly ds t date month, not
timer on 10
connectedness `r(varlist)', model(lasso) lags(1) `f' `opt' nofevd
timer off 10
cn_check, case(gghm2022_lasso_static) tol(1e-8) prefix(R_)

global cn_tests = $cn_tests + 1
capture noisily {
    use ../data/cn_dy2012, clear
    tsset t
    connectedness SP500 R_10Y DJUBSCOM USDX, model(lasso) seed(12345) `opt'
    * a scalar keeps full precision (a local holds a rounded copy in Stata 19)
    tempname t1
    scalar `t1' = e(TCI)
    connectedness SP500 R_10Y DJUBSCOM USDX, model(lasso) seed(12345) `opt'
    assert e(TCI) == scalar(`t1') & "`e(seed)'" == "12345" & "`e(folds)'" == "random"
    assert "`e(alphas)'" == "1" & e(nfolds) == 10 & "`e(loss)'" == "mae"
}
if _rc {
    global cn_fails = $cn_fails + 1
    di as err "  FAIL: elastic net with random folds (r(" _rc "))"
}
else di as txt "  random folds reproducible with seed(): pass"
cn_expect 198 connectedness SP500 R_10Y, model(lasso) alpha(0.5) `opt'
cn_expect 198 connectedness SP500 R_10Y, nfolds(5) `opt'
cn_expect 198 connectedness SP500 R_10Y, model(ridge) loss(rmse) `opt'
cn_expect 198 connectedness SP500 R_10Y, model(elasticnet) folds(blocked) `opt'
cn_expect 198 connectedness SP500 R_10Y, model(lasso) nfolds(2) `opt'

di _n as txt "{hline 78}" _n "13. DCC-GARCH volatility connectedness (0.7.0)" _n "{hline 78}"
* R_<case>_par.csv holds rmgarch's estimates and log likelihood. Each case is checked
* twice: with R's parameters fixed through dccparams() the volatility FEVD must agree
* exactly; with Stata's own maximum likelihood estimates the parameters must agree with
* R's to 1e-4 and the connectedness measures to 1e-3.
capture program drop cn_dcccase
program define cn_dcccase
    syntax varlist, case(string) [opts(string)]
    local f "fixtures/R_`case'_par.csv"
    capture confirm file "`f'"
    if _rc {
        di as txt "  R_`case': no fixture, skipped"
        exit
    }
    local k : word count `varlist'
    tempname P E
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`f'", clear varnames(1) asdouble
        mkmat _all, matrix(`P')
    }
    frame drop __fx
    local llh = `P'[1, 4*`k'+3]
    matrix __rpar = `P'[1, 1..4*`k'+2]
    connectedness `varlist', model(dccgarch) dccparams(__rpar) `opts' frame(cn_test) replace notable
    cn_check, case(`case') tol(1e-8) prefix(R_)
    global cn_tests = $cn_tests + 1
    local d = reldif(e(ll), `llh')
    local ok = (`d' < 1e-10)
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-26s "R_`case'" " log likelihood " as res %9.2e `d' as txt "  " cond(`ok', "pass", "FAIL")
    connectedness `varlist', model(dccgarch) `opts' frame(cn_test) replace notable
    matrix `E' = e(garch)'
    matrix `E' = vec(`E')', e(dcc)
    global cn_tests = $cn_tests + 1
    local d = mreldif(`E', __rpar)
    local ok = (`d' < 1e-4) & (e(converged) == 1)
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-26s "R_`case'" " estimates " as res %9.2e `d' as txt "  " cond(`ok', "pass", "FAIL")
    cn_check, case(`case') tol(1e-3) prefix(R_)
    matrix drop __rpar
end

use ../data/cn_g2020, clear
tsset t
timer on 11
cn_dcccase EUR GBP CHF JPY, case(g2020_dcc)
timer off 11
cn_dcccase EUR GBP CHF JPY, case(g2020_dcc_std_corr) opts(horizon(20) standardize corrected)
use ../data/cn_cgp2024, clear
tsset t
cn_dcccase GRE SOL WIN GEO, case(cgp2024_dcc)

use ../data/cn_g2020, clear
tsset t
local v EUR GBP CHF JPY
cn_expect 198 connectedness `v', model(dccgarch) approach(frequency) `opt'
cn_expect 198 connectedness `v', model(dccgarch) window(500) `opt'
cn_expect 198 connectedness `v', model(dccgarch) orthogonalized `opt'
cn_expect 198 connectedness `v', standardize `opt'
cn_expect 2001 connectedness `v' in 1/80, model(dccgarch) `opt'
matrix __bad = (1, 2, 3)
cn_expect 503 connectedness `v', model(dccgarch) dccparams(__bad) `opt'
matrix drop __bad
cn_expect 111 connectedness `v', model(dccgarch) dccparams(__nosuch) `opt'

di _n as txt "{hline 78}" _n "14. Portfolios and hedging: cnportfolio (0.8.0)" _n "{hline 78}"
* R_pf_<case>.csv holds, per asset or ordered pair, the summary of the weights, hedge
* weights or hedge ratios, the hedging effectiveness, the p-values of five variance
* tests and the annualized return, risk and Sharpe ratio; _ret.csv and _w.csv hold the
* portfolio returns and the weights per period (see test/r/make_fixtures_portfolio.R).
capture mata: mata drop cn_pfdiff()
mata:
real scalar cn_pfdiff(string scalar resframe, string scalar fxframe, string scalar vars)
{
    string scalar  cur
    real colvector per
    real matrix    A, B

    cur = st_framecurrent()
    st_framecurrent(fxframe)
    per = st_data(., "period")
    A   = st_data(., (2..st_nvar()))
    st_framecurrent(resframe)
    B   = st_data(per, tokens(vars))
    st_framecurrent(cur)
    if (cols(A) != cols(B)) return(.)
    return(mreldif(B, A))
}
end
capture program drop cn_pfcase
program define cn_pfcase
    gettoken case 0 : 0
    gettoken tol 0 : 0
    local f "fixtures/R_pf_`case'"
    capture confirm file "`f'.csv"
    if _rc {
        di as txt "  R_pf_`case': no fixture, skipped"
        exit
    }
    quietly cnportfolio `0' saving(cn_pftest, replace) notable
    tempname F T P
    matrix `T' = r(table)
    matrix `P' = r(pvalues)
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`f'.csv", clear varnames(1) asdouble
        mkmat _all, matrix(`F')
    }
    tempname A B
    matrix `A' = `T'[1..., 1..5]
    matrix `B' = `F'[1..., 1..5]
    local d1 = mreldif(`A', `B')
    matrix `A' = `T'[1..., 7..9]
    matrix `B' = `F'[1..., 11..13]
    local d1 = max(`d1', mreldif(`A', `B'))
    matrix `B' = `F'[1..., 6..10]
    local d2 = mreldif(`P', `B')
    local rets
    local ws
    frame cn_pftest {
        foreach v of varlist _all {
            if substr("`v'", 1, 3) == "ret" local rets `rets' `v'
            if inlist(substr("`v'", 1, 2), "w_") | inlist(substr("`v'", 1, 3), "hr_") | substr("`v'", 1, 5) == "beta_" local ws `ws' `v'
        }
    }
    frame __fx: quietly import delimited using "`f'_ret.csv", clear varnames(1) asdouble
    mata: st_local("d3", strofreal(cn_pfdiff("cn_pftest", "__fx", "`rets'"), "%12.0g"))
    frame __fx: quietly import delimited using "`f'_w.csv", clear varnames(1) asdouble
    mata: st_local("d4", strofreal(cn_pfdiff("cn_pftest", "__fx", "`ws'"), "%12.0g"))
    frame drop __fx
    local d3 = max(`d3', `d4')
    local ok = (`d1' < `tol') & (`d2' < `tol') & (`d3' < `tol')
    global cn_tests = $cn_tests + 1
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-22s "R_pf_`case'" " table " as res %8.2e `d1' as txt " p " as res %8.2e `d2' ///
        as txt " series " as res %8.2e `d3' as txt "  " cond(`ok', "pass", "FAIL")
end

* A. static VAR(1), g2020 (decimal returns, as in the R examples)
use ../data/cn_g2020, clear
tsset t
foreach v in EUR GBP CHF JPY {
    quietly generate double r_`v' = `v'/100
}
local r r_EUR r_GBP r_CHF r_JPY
connectedness EUR GBP CHF JPY, lags(1) horizon(10) frame(cn_test) replace notable
timer on 12
cn_pfcase var_mcp_pci        1e-7 `r', method(mcp)
timer off 12
cn_pfcase var_mcp_cov_short  1e-7 `r', method(mcp) matrix(cov) nolong
cn_pfcase var_mcp_cor        1e-7 `r', method(mcp) matrix(cor)
cn_pfcase var_rpp_cov        1e-6 `r', method(rpp)
cn_pfcase var_ewp            1e-7 `r', method(ewp)
cn_pfcase var_bivariate_cov  1e-7 `r', method(bivariate)
cn_pfcase var_hedge_cov      1e-7 `r', method(hedge) cumulative(prod)
cn_pfcase var_mhp_cov        1e-7 `r', method(mhp)
cn_expect 198 cnportfolio r_EUR r_GBP, method(mcp)
cn_expect 198 cnportfolio `r', method(markowitz)
cn_expect 198 cnportfolio `r', method(mcp) test(kruskal)
cn_expect 498 cnportfolio `r', method(mcp) saving(cn_test)
cn_expect 111 cnportfolio `r', method(mcp) frame(cn_nosuch)
connectedness EUR GBP CHF JPY, lags(1) horizon(10) nocovariance frame(cn_test) replace notable
cn_expect 498 cnportfolio `r', method(mhp)
connectedness EUR GBP CHF JPY, approach(r2) lags(1) frame(cn_test) replace notable
cn_expect 498 cnportfolio `r', method(rpp)

* B. TVP-VAR(1), acg2020 (monthly): time-varying covariances and PCI
use ../data/cn_acg2020, clear
tsset month
foreach v in EUR GBP JPY CHF {
    quietly generate double r_`v' = `v'/100
}
local r r_EUR r_GBP r_JPY r_CHF
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) ///
    frame(cn_test) replace notable
cn_pfcase tvp_mcp_pci        1e-6 `r', method(mcp)
cn_pfcase tvp_rpp_cov        1e-6 `r', method(rpp)
cn_pfcase tvp_bivariate_cov  1e-6 `r', method(bivariate) cumulative(prod)
cn_pfcase tvp_hedge_cov      1e-6 `r', method(hedge)
cn_pfcase tvp_mhp_cov        1e-6 `r', method(mhp)

* C. DCC-GARCH(1,1), g2020, with rmgarch's estimates: conditional covariances
use ../data/cn_g2020, clear
tsset t
foreach v in EUR GBP CHF JPY {
    quietly generate double r_`v' = `v'/100
}
local r r_EUR r_GBP r_CHF r_JPY
capture confirm file "fixtures/R_pf_dcc_par.csv"
if !_rc {
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "fixtures/R_pf_dcc_par.csv", clear varnames(1) asdouble
        mkmat _all, matrix(__rpar)
    }
    frame drop __fx
    connectedness EUR GBP CHF JPY, model(dccgarch) dccparams(__rpar) frame(cn_test) replace notable
    matrix drop __rpar
    cn_pfcase dcc_mhp_cov        1e-7 `r', method(mhp)
    cn_pfcase dcc_bivariate_cov  1e-7 `r', method(bivariate)
    cn_pfcase dcc_hedge_cov      1e-7 `r', method(hedge)
    cn_pfcase dcc_mcp_cor        1e-7 `r', method(mcp) matrix(cor)
}
else di as txt "  R_pf_dcc: no fixture, skipped"
capture frame drop cn_pftest

di _n as txt "{hline 78}" _n "15. DCC-GARCH R2 connectedness and estat correlation (0.8.1)" _n "{hline 78}"
* R2 decomposed connectedness of the DCC-GARCH conditional correlations, as R's
* R2Correlations(): exact with rmgarch's estimates (dccparams()), and to 1e-3 with
* Stata's own estimates. estat correlation: conditional and partial correlations of the
* stored covariances, as ConditionalCorrelation() and PartialCorrelations().
capture program drop cn_dccr2case
program define cn_dccr2case
    syntax varlist, case(string)
    local f "fixtures/R_dccr2_`case'_par.csv"
    capture confirm file "`f'"
    if _rc {
        di as txt "  R_dccr2_`case': no fixture, skipped"
        exit
    }
    local k : word count `varlist'
    tempname P
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`f'", clear varnames(1) asdouble
        mkmat _all, matrix(`P')
    }
    frame drop __fx
    matrix __rpar = `P'[1, 1..4*`k'+2]
    connectedness `varlist', model(dccgarch) approach(r2) dccparams(__rpar) frame(cn_test) replace notable
    cn_check, case(dccr2_`case') tol(1e-8) prefix(R_)
    connectedness `varlist', model(dccgarch) approach(r2) frame(cn_test) replace notable
    cn_check, case(dccr2_`case') tol(1e-3) prefix(R_)
    matrix drop __rpar
end
capture program drop cn_corrcase
program define cn_corrcase
    gettoken case 0 : 0
    local f "fixtures/R_corr_`case'"
    capture confirm file "`f'_ts.csv"
    if _rc {
        di as txt "  R_corr_`case': no fixture, skipped"
        exit
    }
    local k = e(k)
    tempname C P M A B
    quietly estat correlation, frame(cn_ctest) replace notable
    matrix `C' = r(mean)
    quietly estat correlation, partial frame(cn_ptest) replace notable
    matrix `P' = r(mean)
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`f'_mean.csv", clear varnames(nonames) asdouble
        mkmat _all, matrix(`M')
    }
    matrix `A' = `M'[1..`k', 1...]
    matrix `B' = `M'[`k'+1..2*`k', 1...]
    local d1 = max(mreldif(`C', `A'), mreldif(`P', `B'))
    frame __fx {
        quietly import delimited using "`f'_ts.csv", clear varnames(1) asdouble
        keep period cor_*
    }
    mata: st_local("d2", strofreal(cn_tsdiff("cn_ctest", "__fx"), "%12.0g"))
    frame __fx {
        quietly import delimited using "`f'_ts.csv", clear varnames(1) asdouble
        keep period pcor_*
    }
    mata: st_local("d3", strofreal(cn_tsdiff("cn_ptest", "__fx"), "%12.0g"))
    frame drop __fx
    local d2 = max(`d2', `d3')
    local ok = (`d1' < 1e-8) & (`d2' < 1e-8)
    global cn_tests = $cn_tests + 1
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-22s "R_corr_`case'" " means " as res %8.2e `d1' as txt "  series " ///
        as res %8.2e `d2' as txt "  " cond(`ok', "pass", "FAIL")
end

use ../data/cn_g2020, clear
tsset t
local v EUR GBP CHF JPY
timer on 13
cn_dccr2case `v', case(g2020)
timer off 13
capture confirm file "fixtures/R_dccr2_g2020_par.csv"
if !_rc {
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "fixtures/R_dccr2_g2020_par.csv", clear varnames(1) asdouble
        mkmat _all, matrix(__rpar)
    }
    frame drop __fx
    matrix __rpar = __rpar[1, 1..18]
    connectedness `v', model(dccgarch) dccparams(__rpar) frame(cn_test) replace notable
    matrix drop __rpar
    cn_corrcase dcc_g2020
}
connectedness `v', lags(1) horizon(10) frame(cn_test) replace notable
cn_corrcase var_g2020
cn_expect 198 connectedness `v', model(dccgarch) approach(r2) cormethod(kendall) `opt'
cn_expect 198 connectedness `v', model(dccgarch) approach(r2) standardize `opt'
cn_expect 198 connectedness `v', model(dccgarch) approach(modelfree) `opt'
connectedness `v', lags(1) nocovariance frame(cn_test) replace notable
cn_expect 498 estat correlation
connectedness `v', approach(r2) lags(1) frame(cn_test) replace notable
cn_expect 498 estat correlation, partial
use ../data/cn_cgp2024, clear
tsset t
cn_dccr2case GRE SOL WIN GEO, case(cgp2024)
use ../data/cn_acg2020, clear
tsset month
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) ///
    frame(cn_test) replace notable
cn_corrcase tvp_acg2020
capture frame drop cn_ctest
capture frame drop cn_ptest

di _n as txt "{hline 78}" _n "16. Normalizations, adaptive elastic net and SAM (0.8.2)" _n "{hline 78}"
* Scalar normalizations of the generalized FEVD (Caloia, Cipollini and Muzzioli 2019),
* adaptive LASSO, ridge and elastic-net VARs (Demirer et al. 2018) and the spillover
* asymmetry measure of user-supplied positive and negative series (Barunik, Kocenda
* and Vacha 2016).
capture program drop cn_assert
program define cn_assert
    args name d tol
    global cn_tests = $cn_tests + 1
    local ok = (`d' < `tol')
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-36s "`name'" " diff " as res %9.2e `d' as txt "  " cond(`ok', "pass", "FAIL")
end
capture program drop cn_pubcheck
program define cn_pubcheck
    * largest absolute difference between TO, FROM, NET and TCI of e(table) and the
    * values published by Caloia, Cipollini and Muzzioli (2019, Table 11), in shares
    args name to from net tci
    tempname T
    matrix `T' = e(table)
    local k = e(k)
    local d = abs(`T'[`k'+3, `k'+1]/100 - `tci')
    forvalues i = 1/`k' {
        local d = max(`d', abs(`T'[`k'+1, `i']/100 - real(word("`to'", `i'))))
        local d = max(`d', abs(`T'[`i', `k'+1]/100 - real(word("`from'", `i'))))
        local d = max(`d', abs(`T'[`k'+3, `i']/100 - real(word("`net'", `i'))))
    }
    cn_assert "`name'" `d' 0.00051
end

use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
foreach n in maxrow maxcol spectral {
    connectedness `v', lags(4) horizon(10) normalize(`n') `opt'
    cn_check, case(dy2012_norm_`n') tol(1e-8) prefix(R_)
}
connectedness `v', lags(4) horizon(10) normalize(maxrow) `opt'
cn_pubcheck "published, max row sum (Table 11)" ".161 .169 .045 .109" ".107 .186 .055 .137" ///
    ".054 -.016 -.010 -.028" .121
connectedness `v', lags(4) horizon(10) `opt'
tempname F0
matrix `F0' = e(fevd)
cn_pubcheck "published, row (Table 11)" ".164 .180 .046 .114" ".112 .186 .063 .143" ///
    ".051 -.005 -.017 -.029" .126
connectedness `v', lags(4) horizon(10) normalize(row) `opt'
tempname F1 FX
matrix `F1' = e(fevd)
cn_assert "normalize(row) equals the default" `=mreldif(`F1', `F0')' 1e-15
cn_expect 198 connectedness `v', lags(4) normalize(maxrow) orthogonalized `opt'
cn_expect 198 connectedness `v', lags(4) normalize(maxrow) approach(frequency) `opt'
cn_expect 198 connectedness `v', lags(4) normalize(columns) `opt'
connectedness `v', lags(4) horizon(10) normalize(maxcol) `opt'
cn_expect 198 estat group, type(conditional) groups(SP500 R_10Y)

use ../data/cn_acg2020, clear
tsset month
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) ///
    normalize(maxrow) `opt'
cn_check, case(acg2020_tvp_norm_maxrow) tol(1e-6) prefix(R_)
use ../data/cn_g2020, clear
tsset t
cn_expect 198 connectedness EUR GBP CHF JPY, model(dccgarch) normalize(maxrow) `opt'

* adaptive elastic net: glmnet with penalty.factor = 1/|b_OLS|, interleaved folds
use ../data/cn_dy2012, clear
tsset t
local f folds(interleaved)
timer on 14
connectedness `v', model(lasso) adaptive lags(2) `f' `opt'
timer off 14
cn_check, case(dy2012_alasso_static) tol(1e-8) prefix(R_)
connectedness `v', model(ridge) adaptive lags(2) `f' `opt'
cn_check, case(dy2012_aridge_static) tol(1e-8) prefix(R_)
connectedness `v', model(elasticnet) alpha(0.5) adaptive lags(2) `f' `opt'
cn_check, case(dy2012_aenet_a05) tol(1e-8) prefix(R_)
connectedness `v', model(elasticnet) nalpha(4) loss(mse) adaptive lags(1) `f' `opt'
cn_check, case(dy2012_aenet4_mse) tol(1e-8) prefix(R_)
connectedness `v' in 1/300, model(lasso) adaptive lags(1) window(200) `f' `opt'
cn_check, case(dy2012_alasso_roll200) tol(1e-8) prefix(R_)
cn_expect 198 connectedness `v', adaptive `opt'

* SAM of positive and negative series: squared positive and negative daily returns
* stand in for realized semivariances
use ../data/cn_g2020, clear
tsset t
foreach x in EUR GBP CHF JPY {
    quietly generate double rv_`x' = `x'^2
    quietly generate double rp_`x' = cond(`x' > 0, `x'^2, 0)
    quietly generate double rn_`x' = cond(`x' < 0, `x'^2, 0)
    quietly generate double p_`x'  = max(`x', 0)
    quietly generate double n_`x'  = min(`x', 0)
}
local rv rv_EUR rv_GBP rv_CHF rv_JPY
local rp rp_EUR rp_GBP rp_CHF rp_JPY
local rn rn_EUR rn_GBP rn_CHF rn_JPY
connectedness `rv', positive(`rp') negative(`rn') window(500) `opt'
local sam = e(SAM)
local tp  = e(TCI_pos)
local tn  = e(TCI_neg)
tempname FP FN
matrix `FP' = e(fevd_pos)
matrix `FN' = e(fevd_neg)
cn_assert "SAM = TCI positive - TCI negative" `=abs(`sam' - (`tp' - `tn'))' 1e-12
* per-period SAM series against the positive and negative frames
frame cn_test {
    quietly frlink 1:1 t, frame(cn_test_pos) generate(__lp)
    quietly frget __tp = tci, from(__lp)
    quietly frget __top = to_1, from(__lp)
    quietly frget __fp = from_1, from(__lp)
    quietly frlink 1:1 t, frame(cn_test_neg) generate(__ln)
    quietly frget __tn = tci, from(__ln)
    quietly frget __ton = to_1, from(__ln)
    quietly frget __fn = from_1, from(__ln)
    quietly generate double __d = max(abs(sam - (__tp - __tn)), abs(sam_to_1 - (__top - __ton)), ///
        abs(sam_from_1 - (__fp - __fn)), abs(sam_net_1 - ((__top - __fp) - (__ton - __fn))))
    quietly summarize __d, meanonly
    local d = r(max)
    drop __lp __ln __tp __tn __top __ton __fp __fn __d
}
cn_assert "SAM series against the part frames" `d' 1e-10
connectedness `rp', window(500) frame(cn_test2) replace notable
matrix `FX' = e(fevd)
cn_assert "positive part = its own estimation" `=max(mreldif(`FX', `FP'), abs(e(TCI) - `tp'))' 1e-12
connectedness `rn', window(500) frame(cn_test2) replace notable
matrix `FX' = e(fevd)
cn_assert "negative part = its own estimation" `=max(mreldif(`FX', `FN'), abs(e(TCI) - `tn'))' 1e-12
capture frame drop cn_test2
* positive() and negative() with the positive and negative changes = asymmetric
connectedness EUR GBP CHF JPY, asymmetric window(500) `opt'
matrix `FP' = e(fevd_pos)
matrix `FN' = e(fevd_neg)
local sam = e(SAM)
connectedness EUR GBP CHF JPY, positive(p_EUR p_GBP p_CHF p_JPY) negative(n_EUR n_GBP n_CHF n_JPY) ///
    window(500) `opt'
tempname GP GN
matrix `GP' = e(fevd_pos)
matrix `GN' = e(fevd_neg)
cn_assert "positive()/negative() = asymmetric" ///
    `=max(mreldif(`GP', `FP'), mreldif(`GN', `FN'), abs(e(SAM) - `sam'))' 1e-12
cn_expect 198 connectedness `rv', positive(`rp') `opt'
cn_expect 198 connectedness `rv', positive(`rp') negative(rn_EUR rn_GBP) `opt'
cn_expect 184 connectedness `rv', positive(`rp') negative(`rn') asymmetric `opt'
cn_expect 198 connectedness `rv', positive(`rp') negative(`rn') approach(r2) `opt'

di _n as txt "{hline 78}" _n "17. Clustered network connectedness (0.9.0)" _n "{hline 78}"
* Buchwalter, Diebold and Yilmaz (2026a). R_<case> fixtures come from the independent
* R implementation in test/r/make_fixtures_090.R; the nested cases are checked
* against the generalized and the orthogonalized decompositions.
capture program drop cn_clcheck
program define cn_clcheck
    syntax , case(string) tol(real)
    cn_check, case(`case') tol(`tol') prefix(R_)
    local f "fixtures/R_`case'"
    capture confirm file "`f'_clts.csv"
    if _rc {
        di as txt "  R_`case' clusters: no fixture, skipped"
        exit
    }
    local C = e(nclusters)
    tempname S T A B S1
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`f'_cl.csv", clear varnames(nonames) asdouble
        mkmat _all, matrix(`S')
        quietly import delimited using "`f'_ctg.csv", clear varnames(nonames) asdouble
        mkmat _all, matrix(`T')
    }
    matrix `A' = e(cl_shares)
    matrix `B' = `A'[1..`C', 1...]
    matrix `S1' = `S'[1..`C', 1...]
    local d1 = mreldif(`B', `S1')
    matrix `B' = `A'[`C'+1, 1..3]
    matrix `S1' = `S'[`C'+1, 1..3]
    local d1 = max(`d1', mreldif(`B', `S1'))
    matrix `B' = e(cl_contagion)
    local d1 = max(`d1', mreldif(`B', `T'))
    frame __fx: quietly import delimited using "`f'_clts.csv", clear varnames(1) asdouble
    mata: st_local("d2", strofreal(cn_tsdiff("`e(frame)'", "__fx"), "%12.0g"))
    frame drop __fx
    local ok = (`d1' < `tol') & (`d2' < `tol')
    global cn_tests = $cn_tests + 1
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-26s "R_`case' clusters" " tables " as res %9.2e `d1' ///
        as txt "  series " as res %9.2e `d2' as txt "  " cond(`ok', "pass", "FAIL")
    cn_cl093, case(`case') tol(`tol')
end

* 0.9.3: cluster table with TO and NET rows (Buchwalter, Diebold and Yilmaz 2026b),
* within/cross split of the TCI (2026a, Fig. 5), against make_fixtures_093.R
capture program drop cn_cl093
program define cn_cl093
    syntax , case(string) tol(real)
    local f "fixtures/R_`case'"
    capture confirm file "`f'_cltab.csv"
    if _rc {
        di as txt "  R_`case' cluster table: no fixture, skipped"
        exit
    }
    local C = e(nclusters)
    tempname R A B Y A1 R1 A2 R2
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`f'_cltab.csv", clear varnames(nonames) asdouble
        mkmat _all, matrix(`R')
        quietly import delimited using "`f'_clsys.csv", clear varnames(nonames) asdouble
        mkmat _all, matrix(`Y')
    }
    matrix `A' = e(cl_table)
    matrix `A1' = `A'[1..., 1..`C']
    matrix `R1' = `R'[1..., 1..`C']
    matrix `A2' = `A'[1..`C', `C'+1]
    matrix `R2' = `R'[1..`C', `C'+1]
    local d1 = max(mreldif(`A1', `R1'), mreldif(`A2', `R2'))
    matrix `B' = (e(cl_within), e(cl_cross))
    local d1 = max(`d1', mreldif(`B', `Y'))
    local d3 = abs(e(cl_within) + e(cl_cross) - e(TCI))
    frame __fx: quietly import delimited using "`f'_clts093.csv", clear varnames(1) asdouble
    mata: st_local("d2", strofreal(cn_tsdiff("`e(frame)'", "__fx"), "%12.0g"))
    frame drop __fx
    local ok = (`d1' < `tol') & (`d2' < `tol') & (`d3' < 1e-10)
    global cn_tests = $cn_tests + 1
    if !`ok' global cn_fails = $cn_fails + 1
    di as txt "  " %-26s "R_`case' cl. table" " tables " as res %9.2e `d1' ///
        as txt "  series " as res %9.2e `d2' as txt "  " cond(`ok', "pass", "FAIL")
end

use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
connectedness `v', lags(4) horizon(10) clusters("A: SP500 R_10Y" "B: DJUBSCOM USDX") `opt'
cn_clcheck, case(dy2012_clust2) tol(1e-8)
tempname S C0
matrix `S' = e(cl_shares)
matrix `C0' = e(cl_contagion)
local d 0
forvalues c = 1/2 {
    local d = max(`d', abs(`S'[`c', 1] + `S'[`c', 2] + `S'[`c', 3] - 100), ///
        abs(`C0'[`c', 1] + `C0'[`c', 2] - 100))
}
local d = max(`d', abs(`S'[3, 1] + `S'[3, 2] + `S'[3, 3] - 100))
cn_assert "own + co-movement + contagion = 100" `d' 1e-10

* nested cases: one cluster = generalized; one variable per cluster = Cholesky,
* in the given order or averaged over all orderings
connectedness `v', lags(4) horizon(10) `opt'
tempname G F
matrix `G' = e(fevd)
connectedness `v', lags(4) horizon(10) clusters("all: `v'") `opt'
matrix `F' = e(fevd)
cn_assert "one cluster = generalized FEVD" `=mreldif(`F', `G')' 1e-12
connectedness `v', lags(4) horizon(10) orthogonalized `opt'
matrix `G' = e(fevd)
connectedness `v', lags(4) horizon(10) clusters("a: SP500" "b: R_10Y" "c: DJUBSCOM" "d: USDX") ///
    clorder(given) `opt'
matrix `F' = e(fevd)
cn_assert "singletons, given order = Cholesky" `=mreldif(`F', `G')' 1e-12
mata: cn__P = cn_perms(4); cn__A = J(4, 4, 0); cn__T = tokens("`v'")
forvalues r = 1/24 {
    mata: st_local("pv", invtokens(cn__T[cn__P[`r', .]]))
    quietly connectedness `pv', lags(4) horizon(10) orthogonalized frame(cn_test) replace notable
    mata: cn__F = st_matrix("e(fevd)"); cn__A[cn__P[`r', .], cn__P[`r', .]] = cn__A[cn__P[`r', .], cn__P[`r', .]] + cn__F
}
mata: st_matrix("`G'", cn__A/24)
mata: mata drop cn__P cn__A cn__F cn__T
connectedness `v', lags(4) horizon(10) clusters("a: SP500") `opt'
matrix `F' = e(fevd)
cn_assert "singletons, averaged = mean Cholesky" `=mreldif(`F', `G')' 1e-12

use ../data/cn_dy2009, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
local reg `""Americas: US ARG BRA CHL MEX" "Europe: UK FRA GER TUR" "AsiaPacific: HKG JPN AUS IDN KOR MYS PHL SGP TAI THA""'
connectedness `v', lags(2) horizon(10) clusters(`reg') `opt'
cn_clcheck, case(dy2009_regions) tol(1e-8)
connectedness `v', lags(2) horizon(10) clusters(`reg') clorder(given) `opt'
cn_clcheck, case(dy2009_regions_given) tol(1e-8)
cn_expect 198 connectedness `v', lags(2) clusters("a: US") `opt'

use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
connectedness `v' in 1/400, lags(2) horizon(10) window(200) clusters("A: SP500" "B: R_10Y DJUBSCOM USDX") `opt'
cn_clcheck, case(dy2012_clust_roll200) tol(1e-8)
cn_expect 198 connectedness `v', clusters("a: SP500 R_10Y") orthogonalized `opt'
cn_expect 198 connectedness `v', clusters("a: SP500 R_10Y") approach(frequency) `opt'
cn_expect 198 connectedness `v', clusters("a: SP500 R_10Y") normalize(maxrow) `opt'
cn_expect 198 connectedness `v', clusters("a: SP500 R_10Y" "b: R_10Y DJUBSCOM") `opt'
cn_expect 198 connectedness `v', clusters("a: SP500 FOO") `opt'
cn_expect 198 connectedness `v', clorder(given) `opt'
cn_expect 198 connectedness `v', clusters("a: SP500 R_10Y") clorder(random) `opt'

use ../data/cn_acg2020, clear
tsset month
connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) ///
    clusters("Europe: EUR GBP CHF" "Japan: JPY") `opt'
cn_clcheck, case(acg2020_tvp_clust) tol(1e-6)

use ../data/cn_cgs2022, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
timer on 15
connectedness `v', lags(1) horizon(10) clusters("CHF: X2YCHF X5YCHF X10YCHF" ///
    "DKK: X2YDKK X5YDKK X10YDKK" "EUR: X2YEUR X5YEUR X10YEUR" "GBP: X2YGBP X5YGBP X10YGBP" ///
    "NOK: X2YNOK X5YNOK X10YNOK" "SEK: X2YSEK X5YSEK X10YSEK") `opt'
timer off 15
cn_clcheck, case(cgs2022_countries) tol(1e-8)

use ../data/cn_g2020, clear
tsset t
connectedness EUR GBP CHF JPY, window(500) asymmetric clusters("Europe: EUR GBP CHF" "Japan: JPY") `opt'
capture frame cn_test_neg: confirm variable cl_contag_2 cl_2_1
cn_assert "clusters with asymmetric" `=_rc' 0.5
cn_expect 198 connectedness EUR GBP CHF JPY, model(dccgarch) clusters("a: EUR GBP") `opt'

di _n as txt "{hline 78}" _n "18. Bootstrap inference and event tests (0.9.1)" _n "{hline 78}"
* Greenwood-Nimmo, Kocenda and Nguyen (2024). Exact checks against the independent R
* implementation in test/r/make_fixtures_091.R, which writes the resampled rows it uses
* (R_boot_idx_*.csv) so that Stata can use the same ones (option testindices());
* statistical checks against the probabilities published by the authors; and their
* simulation with a break.
capture program drop cn_loadidx
program define cn_loadidx
    args file name
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`file'", clear varnames(nonames) asdouble
        mkmat _all, matrix(`name')
    }
    frame drop __fx
end
capture program drop cn_loadfx
program define cn_loadfx
    args file name
    capture frame drop __fx
    frame create __fx
    frame __fx {
        quietly import delimited using "`file'", clear varnames(nonames) asdouble
        mkmat _all, matrix(`name')
    }
    frame drop __fx
end

capture confirm file "fixtures/R_boot_idx_static.csv"
if !_rc {
    use ../data/cn_dy2012, clear
    tsset t
    local v SP500 R_10Y DJUBSCOM USDX
    tempname IDX RT A
    cn_loadidx "fixtures/R_boot_idx_static.csv" `IDX'
    connectedness `v' in 1/800, lags(2) horizon(10) `opt'
    timer on 16
    quietly estat bootstrap, reps(40) level(90) testindices(`IDX') nodots
    timer off 16
    matrix `A' = r(table)
    cn_loadfx "fixtures/R_boot_static_gen.csv" `RT'
    cn_assert "R_boot static, generalized" `=mreldif(`A', `RT')' 1e-8
    connectedness `v' in 1/800, lags(2) horizon(10) orthogonalized `opt'
    quietly estat bootstrap, reps(30) testindices(`IDX') nodots
    matrix `A' = r(table)
    cn_loadfx "fixtures/R_boot_static_orth.csv" `RT'
    cn_assert "R_boot static, orthogonalized" `=mreldif(`A', `RT')' 1e-8

    cn_loadidx "fixtures/R_boot_idx_roll.csv" `IDX'
    connectedness `v' in 1/230, lags(2) horizon(10) window(202) `opt'
    quietly estat bootstrap, reps(25) testindices(`IDX') detect(0) frame(cn_btest) replace nodots
    capture frame drop __fx
    frame create __fx
    frame __fx: quietly import delimited using "fixtures/R_boot_roll_ts.csv", clear varnames(1) asdouble
    mata: st_local("d2", strofreal(cn_tsdiff("cn_btest", "__fx"), "%12.0g"))
    frame drop __fx
    cn_assert "R_boot rolling bands and detect" `d2' 1e-8
    cn_loadfx "fixtures/R_boot_roll_event.csv" `RT'
    tempname E1 E2 R1 R2
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) nodots
    matrix `E1' = r(prob)
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) reference(3) ///
        threshold(2) decrease nodots
    matrix `E2' = r(prob)
    matrix `R1' = `RT'[1, 1...]
    matrix `R2' = `RT'[2, 1...]
    cn_assert "R_boot event probabilities" `=max(mreldif(`E1', `R1'), mreldif(`E2', `R2'))' 1e-8
    * a fixed seed reproduces the draws
    quietly estat bootstrap, reps(20) seed(11) nodots frame(cn_btest) replace every(10)
    frame cn_btest: summarize tci_lo, meanonly
    local s1 = r(mean)
    quietly estat bootstrap, reps(20) seed(11) nodots frame(cn_btest) replace every(10)
    frame cn_btest: summarize tci_lo, meanonly
    cn_assert "seed() reproduces the bootstrap" `=abs(r(mean) - `s1')' 1e-15
    cn_expect 198 estat bootstrap, reps(1)
    cn_expect 198 estat bootstrap, leads(0 1)
    cn_expect 198 estat bootstrap, detect(5) every(2)
    cn_expect 198 estat bootstrap, event(01jan2050)
    connectedness `v' in 1/800, lags(2) `opt'
    cn_expect 198 estat bootstrap, event(10mar1999)
    connectedness `v' in 1/300, model(tvpvar) `opt'
    cn_expect 198 estat bootstrap, reps(10)
    connectedness `v' in 1/300, lags(1) approach(frequency) `opt'
    cn_expect 198 estat bootstrap, reps(10)
    capture frame drop cn_btest
}
else di as txt "  R_boot: no fixture, skipped"

* the authors' daily data: probabilities published for four events (generalized
* FEVD, VAR(2), horizon 10, 200-day windows, 1,000 draws per round); bootstrap
* probabilities differ by Monte Carlo error only
use ../data/cn_gkn2024, clear
tsset t
connectedness US-TUR, lags(2) horizon(10) window(202) datevar(date) `opt'
timer on 17
quietly estat bootstrap, event(17oct1997 23jan2001 11sep2001 09aug2007) leads(0 1 5 10 22) ///
    seed(1357) nodots
timer off 17
tempname G P
matrix `G' = r(prob)
matrix `P' = (97.1, 84.9, 100, 100, 100 \ 100, 0.1, 0.1, 0.7, 0.6 \ ///
    99.2, 100, 100, 100, 100 \ 48.3, 57.3, 47.7, 100, 100)
local dmax 0
local dsum 0
forvalues i = 1/4 {
    forvalues j = 1/5 {
        local d = abs(`G'[`i', `j'] - `P'[`i', `j'])
        local dmax = max(`dmax', `d')
        local dsum = `dsum' + `d'
    }
}
di as txt "  published event probabilities: largest difference " %4.1f `dmax' ///
    ", mean " %4.2f `dsum'/20 " percentage points"
cn_assert "published probabilities (max |diff| < 8)" `dmax' 8
cn_assert "published probabilities (mean |diff| < 3)" `=`dsum'/20' 3

* the authors' simulation (19 variables, VAR(1), break at observation 301): the
* rise is detected on impact with the permanent and the temporary shift
foreach s in perm temp {
    capture frame drop __sim
    import delimited using "fixtures/gkn_simul_`s'.csv", clear varnames(1) asdouble
    tsset t
    connectedness y1-y19, lags(1) horizon(10) window(201) `opt'
    quietly estat bootstrap, event(301) leads(0) reps(500) seed(2468) nodots
    tempname S
    matrix `S' = r(prob)
    di as txt "  simulation (`s' shift): probability on impact " %5.1f `S'[1, 1]
    cn_assert "simulated break detected (`s')" `=cond(`S'[1, 1] >= 90, 0, 1)' 0.5
}

di _n as txt "{hline 78}" _n "19. Bootstrap of pairwise and NPDC measures; stat() (0.9.2)" _n "{hline 78}"
* Exact checks against the independent R implementation in test/r/make_fixtures_092.R
* (the same resampled rows as section 18): bands for the pairwise table and the NPDC
* matrix of ConnectednessApproach, static and rolling, and event and detection
* probabilities for TO, NET, pairwise and NPDC measures; identities that hold draw by
* draw (FROM, TO and NET from the pairwise table, NPDC = CT - CT').
capture confirm file "fixtures/R_boot_static_pw.csv"
if !_rc {
    use ../data/cn_dy2012, clear
    tsset t
    local v SP500 R_10Y DJUBSCOM USDX
    tempname IDX RT A Tb pw_b pw_bias pw_bc pw_se pw_lo pw_hi npdc_b npdc_bias npdc_bc ///
        npdc_se npdc_lo npdc_hi
    cn_loadidx "fixtures/R_boot_idx_static.csv" `IDX'
    connectedness `v' in 1/800, lags(2) horizon(10) `opt'
    quietly estat bootstrap, reps(40) level(90) testindices(`IDX') pairwise npdc nodots
    matrix `Tb' = r(table)
    foreach blk in pw npdc {
        foreach s in b bias bc se lo hi {
            matrix ``blk'_`s'' = r(`blk'_`s')
        }
    }
    cn_loadfx "fixtures/R_boot_static_gen.csv" `RT'
    cn_assert "R_boot static table, with pairwise" `=mreldif(`Tb', `RT')' 1e-8
    foreach blk in pw npdc {
        mata: st_matrix("`A'", (vec(st_matrix("``blk'_b'")'), vec(st_matrix("``blk'_bias'")'), ///
            vec(st_matrix("``blk'_bc'")'), vec(st_matrix("``blk'_se'")'), ///
            vec(st_matrix("``blk'_lo'")'), vec(st_matrix("``blk'_hi'")')))
        cn_loadfx "fixtures/R_boot_static_`blk'.csv" `RT'
        cn_assert "R_boot static `blk'" `=mreldif(`A', `RT')' 1e-8
    }
    * identities of the bias-corrected estimates (linear in the draws)
    mata: cn__T = st_matrix("`Tb'")[., 3]; cn__P = st_matrix("`pw_bc'"); cn__N = st_matrix("`npdc_bc'")
    mata: st_local("d1", strofreal(max(abs(rowsum(cn__P) - diagonal(cn__P) - cn__T[|6 \ 9|])), "%12.0g"))
    mata: st_local("d2", strofreal(max(abs(colsum(cn__P)' - diagonal(cn__P) - cn__T[|2 \ 5|])), "%12.0g"))
    mata: st_local("d3", strofreal(max(abs(colsum(cn__N)' - cn__T[|10 \ 13|])), "%12.0g"))
    mata: st_local("d4", strofreal(max(abs(cn__N - (cn__P - cn__P'))), "%12.0g"))
    mata: st_local("d5", strofreal(max(abs(st_matrix("`npdc_lo'") + st_matrix("`npdc_hi'")')), "%12.0g"))
    mata: mata drop cn__T cn__P cn__N
    cn_assert "pairwise rows less own = FROM" `d1' 1e-10
    cn_assert "pairwise columns less own = TO" `d2' 1e-10
    cn_assert "NPDC columns = NET" `d3' 1e-10
    cn_assert "NPDC = pairwise - pairwise'" `d4' 1e-10
    cn_assert "NPDC intervals antisymmetric" `d5' 1e-10

    cn_loadidx "fixtures/R_boot_idx_roll.csv" `IDX'
    connectedness `v' in 1/230, lags(2) horizon(10) window(202) `opt'
    quietly estat bootstrap, reps(25) testindices(`IDX') pairwise npdc frame(cn_btest) replace nodots
    cn_assert "default pairwise frame" `=("`r(pwframe)'" != "cn_btest_pw")' 0.5
    frame cn_btest_pw {
        mkmat pw-npdc_hi, matrix(`A')
        local okn = (from[1] == "SP500" & to[1] == "R_10Y" & from[_N] == "USDX" & ///
            to[_N] == "DJUBSCOM" & _N == 31*12)
    }
    cn_loadfx "fixtures/R_boot_roll_pw.csv" `RT'
    matrix `RT' = `RT'[1..., 4...]
    cn_assert "R_boot rolling pairwise and NPDC" `=mreldif(`A', `RT')' 1e-8
    cn_assert "pairwise frame layout" `=1 - `okn'' 0.5

    tempname E E2
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) ///
        stat(npdc SP500 R_10Y) nodots
    matrix `E' = r(prob)
    local st `r(stat)'
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) ///
        stat(pairwise USDX SP500) reference(3) threshold(2) nodots
    matrix `E' = `E' \ r(prob)
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) ///
        stat(net DJUBSCOM) decrease nodots
    matrix `E' = `E' \ r(prob)
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) ///
        stat(to R_10Y) threshold(1) nodots
    matrix `E' = `E' \ r(prob)
    cn_loadfx "fixtures/R_boot_roll_statevent.csv" `RT'
    cn_assert "R_boot events, stat()" `=mreldif(`E', `RT')' 1e-8
    cn_assert "r(stat)" `=("`st'" != "npdc SP500 R_10Y")' 0.5
    quietly estat bootstrap, reps(25) testindices(`IDX') detect(0) stat(npdc SP500 R_10Y) ///
        frame(cn_btest) replace nodots
    frame cn_btest: mkmat prob, matrix(`E')
    quietly estat bootstrap, reps(25) testindices(`IDX') detect(2) stat(from USDX) ///
        frame(cn_btest) replace nodots
    frame cn_btest: mkmat prob, matrix(`E2')
    matrix `E' = `E', `E2'
    cn_loadfx "fixtures/R_boot_roll_statdetect.csv" `RT'
    cn_assert "R_boot detect, stat()" `=mreldif(`E', `RT')' 1e-8

    cn_expect 198 estat bootstrap, reps(5) stat(net SP500)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(npdc SP500 R_10Y) threshold(2)
    cn_expect 198 estat bootstrap, reps(5) detect(2) stat(net SP500)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(level SP500)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(npdc SP500 SP500)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(npdc SP500)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(to XYZ)
    cn_expect 198 estat bootstrap, reps(5) event(209) pairwise
    capture frame drop cn_btest cn_btest_pw

    * clustered identification: identities hold, NPDC table displayed
    connectedness `v' in 1/800, lags(2) clusters("A: SP500 R_10Y" "B: DJUBSCOM USDX") `opt'
    capture noisily estat bootstrap, reps(10) seed(5) pairwise npdc nodots
    local rc = _rc
    if !`rc' {
        mata: cn__T = st_matrix("r(table)")[., 3]; cn__N = st_matrix("r(npdc_bc)")
        mata: st_local("d3", strofreal(max(abs(colsum(cn__N)' - cn__T[|10 \ 13|])), "%12.0g"))
        mata: mata drop cn__T cn__N
    }
    else local d3 = 1
    cn_assert "clusters: NPDC columns = NET" `d3' 1e-10
    cn_expect 198 estat bootstrap, reps(5) pwframe(x) pairwise
}
else di as txt "  R_boot pairwise: no fixture, skipped"

di _n as txt "{hline 78}" _n "20. Cluster table and within/cross split of the TCI (0.9.3)" _n "{hline 78}"
* Identities of the measures added to align with the published version of
* Buchwalter, Diebold and Yilmaz (2026a) and their follow-up (2026b); the R
* comparisons run with each case of section 17 (R_<case> cl. table).
use ../data/cn_dy2009, clear
tsset t
quietly ds t date, not
local v `r(varlist)'
local reg `""Americas: US ARG BRA CHL MEX" "Europe: UK FRA GER TUR" "AsiaPacific: HKG JPN AUS IDN KOR MYS PHL SGP TAI THA""'
connectedness `v', lags(2) horizon(10) clusters(`reg') `opt'
tempname S T I
matrix `S' = e(cl_shares)
matrix `T' = e(cl_table)
matrix `I' = e(cl_id)
mata: cn__S = st_matrix("`S'"); cn__T = st_matrix("`T'"); cn__I = st_matrix("`I'")'
mata: cn__n = J(3, 1, .); for (cn__c=1; cn__c<=3; cn__c++) cn__n[cn__c] = sum(cn__I :== cn__c)
mata: st_local("d1", strofreal(abs(sum(cn__n :* cn__S[(1..3), 2])/19 - st_numscalar("e(cl_within)")) + abs(sum(cn__n :* cn__S[(1..3), 3])/19 - st_numscalar("e(cl_cross)")), "%12.0g"))
mata: st_local("d2", strofreal(max(abs(rowsum(cn__T[(1..3), (1..3)]) :- 100)), "%12.0g"))
mata: st_local("d3", strofreal(max(abs(cn__T[(1..3), 4] - cn__S[(1..3), 3])) + abs(sum(cn__T[5, (1..3)])), "%12.0g"))
mata: st_local("d4", strofreal(max(abs(cn__T[5, (1..3)] - (cn__T[4, (1..3)] - cn__T[(1..3), 4]'))), "%12.0g"))
mata: mata drop cn__S cn__T cn__I cn__n cn__c
cn_assert "within/cross = size-weighted shares" `d1' 1e-10
cn_assert "cluster table rows sum to 100" `d2' 1e-10
cn_assert "FROM = contagion; NET sums to zero" `d3' 1e-10
cn_assert "NET = TO - FROM" `d4' 1e-10
* rolling: within + cross = TCI in every period
use ../data/cn_dy2012, clear
tsset t
connectedness SP500 R_10Y DJUBSCOM USDX in 1/400, lags(2) window(202) ///
    clusters("A: SP500" "B: R_10Y DJUBSCOM USDX") `opt'
frame `e(frame)' {
    tempvar dd
    quietly generate double `dd' = abs(cl_within + cl_cross - tci)
    quietly summarize `dd'
    local d5 = r(max)
}
cn_assert "within + cross = TCI, every period" `d5' 1e-10
* one cluster: everything is within; one variable per cluster: everything is cross
connectedness SP500 R_10Y DJUBSCOM USDX, lags(2) clusters("all: SP500 R_10Y DJUBSCOM USDX") `opt'
cn_assert "one cluster: cross = 0, within = TCI" `=abs(e(cl_cross)) + abs(e(cl_within) - e(TCI))' 1e-10
connectedness SP500 R_10Y DJUBSCOM USDX, lags(2) clusters("a: SP500" "b: R_10Y" "c: DJUBSCOM" "d: USDX") `opt'
cn_assert "singletons: within = 0, cross = TCI" `=abs(e(cl_within)) + abs(e(cl_cross) - e(TCI))' 1e-10

di _n as txt "{hline 78}" _n "21. Bootstrap of the spillover asymmetry measure (0.9.4)" _n "{hline 78}"
* Exact checks against the independent R implementation in test/r/make_fixtures_094.R,
* which bootstraps the total, positive and negative VARs jointly (the same resampled
* rows in every system, a draw discarded if any VAR is unstable): static SAM table,
* rolling bands with positive() and negative(), event and detection probabilities of
* SAM measures. Identities: point estimates against connectedness, bias correction of
* the difference, the total system unchanged by the joint draws, and an exact mirror
* case (negative = - positive), where every draw of the SAM is zero.
capture confirm file "fixtures/R_sam_static_sam.csv"
if !_rc {
    use ../data/cn_g2020, clear
    tsset t
    local v EUR GBP CHF JPY
    tempname IDX RT A Tb Sm Sp E E2
    cn_loadidx "fixtures/R_boot_idx_sam_static.csv" `IDX'
    connectedness `v' in 1/600, lags(2) horizon(10) asymmetric `opt'
    local esam = e(SAM)
    local etp  = e(TCI_pos)
    local etn  = e(TCI_neg)
    local etci = e(TCI)
    timer on 18
    quietly estat bootstrap, sam reps(40) level(90) testindices(`IDX') nodots
    timer off 18
    matrix `Tb' = r(table)
    matrix `Sm' = r(sam)
    matrix `Sp' = r(tci_parts)
    local sp = r(sam_prob)
    local nd = r(ndiscard)
    cn_loadfx "fixtures/R_sam_static_tot.csv" `RT'
    cn_assert "R_sam static, total system" `=mreldif(`Tb', `RT')' 1e-8
    cn_loadfx "fixtures/R_sam_static_sam.csv" `RT'
    cn_assert "R_sam static, SAM table" `=mreldif(`Sm', `RT')' 1e-8
    cn_loadfx "fixtures/R_sam_static_parts.csv" `RT'
    cn_assert "R_sam static, TCI of each part" `=mreldif(`Sp', `RT')' 1e-8
    cn_assert "estimates = connectedness (SAM, TCI+, TCI-, TCI)" ///
        `=abs(`Sm'[1, 1] - `esam') + abs(`Sp'[1, 1] - `etp') + abs(`Sp'[2, 1] - `etn') + abs(`Tb'[1, 1] - `etci')' 1e-10
    cn_assert "bias and corrected SAM = differences of the parts" ///
        `=abs(`Sm'[1, 2] - (`Sp'[1, 2] - `Sp'[2, 2])) + abs(`Sm'[1, 3] - (`Sp'[1, 3] - `Sp'[2, 3]))' 1e-10
    mata: cn__S = st_matrix("`Sm'")
    mata: st_local("d1", strofreal(abs(sum(cn__S[(2..5), 3]) - 4*cn__S[1, 3]) + abs(sum(cn__S[(6..9), 3]) - 4*cn__S[1, 3]) + abs(sum(cn__S[(10..13), 3])), "%12.0g"))
    mata: st_local("d2", strofreal(max(abs(cn__S[(10..13), 3] - (cn__S[(2..5), 3] - cn__S[(6..9), 3]))), "%12.0g"))
    mata: mata drop cn__S
    cn_assert "SAM TO and FROM sum to k SAM; NET to zero" `d1' 1e-10
    cn_assert "SAM NET = SAM TO - SAM FROM" `d2' 1e-10
    cn_assert "r(sam_prob)" `=abs(`sp' - `Sm'[1, 7])' 1e-12
    * the total system from the joint draws = its own bootstrap on the same rows (no
    * draw discarded)
    quietly estat bootstrap, reps(40) level(90) testindices(`IDX') nodots
    matrix `A' = r(table)
    cn_assert "total system: joint = alone (no discards)" `=cond(`nd' == 0 & r(ndiscard) == 0, mreldif(`A', `Tb'), 1)' 1e-12
    * the same with positive() and negative() given as the parts
    foreach x of local v {
        quietly generate double p_`x' = max(`x', 0)
        quietly generate double n_`x' = min(`x', 0)
        quietly generate double m_`x' = -p_`x'
    }
    connectedness `v' in 1/600, lags(2) horizon(10) positive(p_*) negative(n_*) `opt'
    quietly estat bootstrap, sam reps(40) level(90) testindices(`IDX') nodots
    matrix `A' = r(sam)
    cn_assert "positive()/negative() = asymmetric" `=mreldif(`A', `Sm')' 1e-12

    * mirror case: the negative series is minus the positive series, so the two systems
    * coincide up to sign and every draw of SAM and its TO, FROM and NET is zero
    connectedness `v' in 1/400, lags(1) horizon(10) positive(p_*) negative(m_*) `opt'
    quietly estat bootstrap, sam reps(30) seed(4) nodots
    matrix `A' = r(sam)
    matrix `A' = `A'[1..., 1..6]
    mata: st_local("d1", strofreal(max(abs(st_matrix("`A'"))), "%12.0g"))
    cn_assert "mirror: SAM bands exactly zero" `d1' 1e-10
    matrix `A' = r(tci_parts)
    matrix `E' = `A'[1, 1..6]
    matrix `E2' = `A'[2, 1..6]
    cn_assert "mirror: TCI+ = TCI-" `=mreldif(`E', `E2')' 1e-12
    connectedness `v' in 1/500, lags(1) horizon(10) window(301) positive(p_*) negative(m_*) `opt'
    quietly estat bootstrap, sam all reps(20) every(50) seed(4) frame(cn_btest) replace nodots
    local ns = r(nsam)
    frame cn_btest {
        tempvar dd
        quietly generate double `dd' = max(abs(sam), abs(sam_bc), abs(sam_lo), abs(sam_hi), ///
            abs(sam_net_2_lo), abs(sam_to_4_hi), abs(sam_prob))
        quietly summarize `dd'
        local d1 = r(max)
    }
    cn_assert "mirror, rolling: SAM bands exactly zero" `d1' 1e-10
    cn_assert "mirror, rolling: r(nsam)" `ns' 0.5

    * rolling, positive() and negative(): squared returns and squared positive and
    * negative returns
    foreach x of local v {
        quietly generate double rv_`x' = `x'^2
        quietly generate double rp_`x' = cond(`x' > 0, `x'^2, 0)
        quietly generate double rn_`x' = cond(`x' < 0, `x'^2, 0)
    }
    cn_loadidx "fixtures/R_boot_idx_sam_roll.csv" `IDX'
    connectedness rv_* in 1/240, lags(1) horizon(10) window(201) positive(rp_*) negative(rn_*) `opt'
    quietly estat bootstrap, sam all reps(25) testindices(`IDX') frame(cn_btest) replace nodots
    frame cn_btest: mkmat tci-tci_hi sam-tci_neg_hi ndiscard stable, matrix(`A')
    cn_loadfx "fixtures/R_sam_roll.csv" `RT'
    cn_assert "R_sam rolling bands" `=mreldif(`A', `RT')' 1e-8
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) stat(sam) nodots
    matrix `E' = r(prob)
    local st `r(stat)'
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) ///
        stat(sam net rv_CHF) reference(3) decrease nodots
    matrix `E' = `E' \ r(prob)
    local st2 `r(stat)'
    quietly estat bootstrap, reps(25) testindices(`IDX') event(209) leads(0 1 5) ///
        stat(sam to rv_EUR) nodots
    matrix `E' = `E' \ r(prob)
    cn_loadfx "fixtures/R_sam_roll_event.csv" `RT'
    cn_assert "R_sam events, stat(sam)" `=mreldif(`E', `RT')' 1e-8
    cn_assert "r(stat) with sam" `=("`st'" != "sam") + ("`st2'" != "sam net rv_CHF")' 0.5
    quietly estat bootstrap, reps(25) testindices(`IDX') detect(0) stat(sam) ///
        frame(cn_btest) replace nodots
    frame cn_btest: mkmat prob, matrix(`E')
    quietly estat bootstrap, reps(25) testindices(`IDX') detect(0) stat(sam from rv_GBP) ///
        frame(cn_btest) replace nodots
    frame cn_btest: mkmat prob, matrix(`E2')
    matrix `E' = `E', `E2'
    cn_loadfx "fixtures/R_sam_roll_detect.csv" `RT'
    cn_assert "R_sam detect, stat(sam)" `=mreldif(`E', `RT')' 1e-8

    cn_expect 198 estat bootstrap, reps(5) event(209) stat(sam) threshold(2)
    cn_expect 198 estat bootstrap, reps(5) detect(2) stat(sam)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(sam level rv_EUR)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(sam net XYZ)
    cn_expect 198 estat bootstrap, reps(5) event(209) stat(sam net)
    capture frame drop cn_btest

    * clustered and orthogonalized decompositions: estimates as connectedness, bias of
    * the difference
    foreach o in `"clusters("E: EUR GBP CHF" "J: JPY")"' orthogonalized {
        connectedness `v' in 1/600, lags(2) asymmetric `o' `opt'
        local esam = e(SAM)
        capture noisily estat bootstrap, sam reps(10) seed(5) nodots
        if !_rc {
            matrix `Sm' = r(sam)
            matrix `Sp' = r(tci_parts)
            local d1 = abs(`Sm'[1, 1] - `esam') + abs(`Sm'[1, 3] - (`Sp'[1, 3] - `Sp'[2, 3]))
        }
        else local d1 = 1
        gettoken o1 : o, parse("(")
        cn_assert "`o1': SAM estimate and bias" `d1' 1e-10
    }
    use ../data/cn_dy2012, clear
    tsset t
    connectedness SP500 R_10Y DJUBSCOM USDX in 1/300, lags(2) `opt'
    cn_expect 198 estat bootstrap, reps(5) sam
}
else di as txt "  R_sam: no fixture, skipped"

di _n as txt "{hline 78}" _n "22. Control cluster, ordering diagnostic, simultaneous orthogonalization (0.9.5)" _n "{hline 78}"
* Exact checks against the independent R implementation in test/r/make_fixtures_095.R
* (control cluster by Schur complement, Buchwalter, Diebold and Yilmaz 2026b, eq. 7;
* net connectedness under every ordering and the first-versus-last diagnostic, eq. 4;
* simultaneous clustered orthogonalization, 2026c, eq. 24), and identities.
capture program drop cn_cl095
program define cn_cl095
    syntax , case(string) tol(real)
    local f "fixtures/R_`case'"
    local C = e(nclusters)
    local nfix = (e(cl_ncontrols) == 1)
    local Cp = `C' - `nfix'
    tempname R A B R1 A1
    capture confirm file "`f'_clord.csv"
    if !_rc {
        * net connectedness by ordering, matched on the positions of the clusters
        capture frame drop __fx
        frame copy `e(frame_clord)' __fx
        frame __fx {
            sort pos_*
            mkmat pos_* net_*, matrix(`A')
            quietly import delimited using "`f'_clord.csv", clear varnames(nonames) asdouble
            local sv
            forvalues c = 1/`C' {
                local sv `sv' v`c'
            }
            sort `sv'
            mkmat _all, matrix(`R')
        }
        frame drop __fx
        local d1 = mreldif(`A', `R')
        quietly estat clorder
        matrix `A' = r(net)
        matrix `B' = r(bypos)
        matrix `A1' = r(dfl)
        capture frame drop __fx
        frame create __fx
        frame __fx: quietly import delimited using "`f'_cldiag.csv", clear varnames(nonames) asdouble
        frame __fx: mkmat _all, matrix(`R')
        frame drop __fx
        local r1 = `C' + 1
        local r2 = `C' + `Cp'
        local r3 = `C' + `Cp' + 1
        local r4 = `C' + 2*`Cp'
        matrix `R1' = `R'[1..`C', 1..5]
        local d2 = mreldif(`A', `R1')
        matrix `R1' = `R'[`r1'..`r2', 1..`Cp']
        local d2 = max(`d2', mreldif(`B', `R1'))
        * pairs: sd missing with a single pair
        matrix `R1' = `R'[`r3'..`r4', 1..5]
        mata: cn__a = st_matrix("`A1'"); cn__r = st_matrix("`R1'")
        mata: st_local("d3", strofreal(any(missing(cn__a) :!= missing(cn__r)) + max(reldif(editmissing(cn__a, 0), editmissing(cn__r, 0))), "%12.0g"))
        mata: mata drop cn__a cn__r
        local ok = (`d1' < `tol') & (`d2' < `tol') & (`d3' < `tol')
        global cn_tests = $cn_tests + 1
        if !`ok' global cn_fails = $cn_fails + 1
        di as txt "  " %-26s "R_`case' orderings" " nets " as res %9.2e `d1' ///
            as txt "  clorder " as res %9.2e `=max(`d2', `d3')' as txt "  " cond(`ok', "pass", "FAIL")
    }
    capture confirm file "`f'_clts095.csv"
    if !_rc {
        capture frame drop __fx
        frame create __fx
        frame __fx: quietly import delimited using "`f'_clts095.csv", clear varnames(1) asdouble
        frame __fx {
            mata: cn__per = st_data(., "period"); cn__r = st_data(., "cl_dfl_* cl_nsd_*")
        }
        frame `e(frame)': mata: cn__a = st_data(cn__per, "cl_dfl_* cl_nsd_*")
        frame drop __fx
        mata: st_local("d1", strofreal(any(missing(cn__a) :!= missing(cn__r)) + max(reldif(editmissing(cn__a, 0), editmissing(cn__r, 0))), "%12.0g"))
        mata: mata drop cn__per cn__r cn__a
        local ok = (`d1' < `tol')
        global cn_tests = $cn_tests + 1
        if !`ok' global cn_fails = $cn_fails + 1
        di as txt "  " %-26s "R_`case' ordering series" " diff " as res %9.2e `d1' ///
            as txt "  " cond(`ok', "pass", "FAIL")
    }
end

capture confirm file "fixtures/R_dy2009_ctl_clord.csv"
if !_rc {
    use ../data/cn_dy2009, clear
    tsset t
    quietly ds t date, not
    local v `r(varlist)'
    local reg3 `""Americas: ARG BRA CHL MEX" "Europe: UK FRA GER TUR" "AsiaPacific: HKG JPN AUS IDN KOR MYS PHL SGP TAI THA""'
    connectedness `v', lags(2) horizon(10) clcontrols(US) clusters(`reg3') `opt'
    cn_assert "control cluster: 6 orderings, e(clcontrols)" ///
        `=(e(norderings) != 6) + ("`e(clcontrols)'" != "US") + (e(cl_ncontrols) != 1) + ("`: word 1 of `e(clusters)''" != "Controls")' 0.5
    cn_clcheck, case(dy2009_ctl) tol(1e-8)
    cn_cl095, case(dy2009_ctl) tol(1e-8)
    * identities of the orderings frame: nets sum to zero, mean over orderings =
    * averaged decomposition, (Cp - 1)! pairs
    frame `e(frame_clord)' {
        tempvar s
        quietly egen double `s' = rowtotal(net_*)
        quietly replace `s' = abs(`s')
        quietly summarize `s'
        local d1 = r(max)
        local no = _N
    }
    cn_assert "nets of an ordering sum to zero" `d1' 1e-10
    quietly estat clorder, frame(cn_ctest) replace
    tempname N D
    matrix `N' = r(net)
    matrix `D' = r(dfl)
    mata: st_local("d2", strofreal(max(abs(st_matrix("`N'")[., 1] - st_matrix("`N'")[., 2])), "%12.0g"))
    cn_assert "mean net over orderings = net of the average" `d2' 1e-10
    tempname D5 J2
    matrix `D5' = `D'[1..., 5]
    matrix `J2' = J(3, 1, 2)
    cn_assert "pairs: (Cp - 1)! per cluster" `=mreldif(`D5', `J2')' 1e-12
    frame cn_ctest {
        quietly count if dfl < .
        local nd = r(N)
        quietly summarize dfl if cluster == "Europe", meanonly
        local md = r(mean)
        local nl = _N
    }
    cn_assert "long frame: rows and paired differences" `=(`nl' != 6*4) + (`nd' != 3*2)' 0.5
    cn_assert "long frame: mean paired difference" `=abs(`md' - `D'[2, 1])' 1e-10
    capture frame drop cn_ctest
    * control first with clorder(given) = the control cluster listed first, given order
    connectedness `v', lags(2) horizon(10) clcontrols(US) clusters(`reg3') clorder(given) `opt'
    cn_clcheck, case(dy2009_ctl_given) tol(1e-8)
    tempname F1 F2
    matrix `F1' = e(fevd)
    connectedness `v', lags(2) horizon(10) clusters("Controls: US" `reg3') clorder(given) `opt'
    matrix `F2' = e(fevd)
    cn_assert "clcontrols(), given = control listed first" `=mreldif(`F1', `F2')' 1e-12
    cn_expect 198 estat clorder
    * three regions without controls, ordering diagnostic
    local reg `""Americas: US ARG BRA CHL MEX" "Europe: UK FRA GER TUR" "AsiaPacific: HKG JPN AUS IDN KOR MYS PHL SGP TAI THA""'
    connectedness `v', lags(2) horizon(10) clusters(`reg') `opt'
    cn_clcheck, case(dy2009_reg_diag) tol(1e-8)
    cn_cl095, case(dy2009_reg_diag) tol(1e-8)
    * simultaneous orthogonalization of the three regions
    connectedness `v', lags(2) horizon(10) clusters(`reg') clorder(simultaneous) `opt'
    cn_assert "simultaneous: one ordering, no orderings frame" `=(e(norderings) != 1) + ("`e(frame_clord)'" != "")' 0.5
    cn_clcheck, case(dy2009_sco) tol(1e-8)
    cn_expect 198 estat clorder
    * seven clusters, 5,040 orderings
    local reg7 `""USA: US" "Latam: ARG BRA CHL MEX" "UK: UK" "Euro: FRA GER TUR" "Japan: JPN" "EastAsia: HKG KOR TAI SGP" "RestAP: AUS IDN MYS PHL THA""'
    timer on 19
    connectedness `v', lags(2) horizon(10) clusters(`reg7') `opt'
    timer off 19
    cn_assert "seven clusters: 5,040 orderings" `=(e(norderings) != 5040)' 0.5
    cn_clcheck, case(dy2009_seven) tol(1e-8)
    cn_cl095, case(dy2009_seven) tol(1e-8)
    * cap: nine clusters to order are refused, eight and a control cluster are not
    cn_expect 198 connectedness `v', lags(2) clusters("a: US" "b: ARG" "c: BRA" "d: CHL" "e: MEX" ///
        "f: UK" "g: FRA" "h: GER" "i: TUR JPN HKG AUS IDN KOR MYS PHL SGP TAI THA") `opt'
    cn_expect 198 connectedness `v', lags(2) clcontrols(US) clusters(`reg3') clorder(simultaneous) `opt'
    cn_expect 198 connectedness `v', lags(2) clcontrols(US) clusters("Americas: US ARG BRA") `opt'
    cn_expect 198 connectedness `v', lags(2) clcontrols(`v') `opt'
    cn_expect 198 connectedness `v', lags(2) clcontrols(US FOO) clusters(`reg3') `opt'
    cn_expect 198 connectedness `v', lags(2) clusters(`reg') clorder(sequential) `opt'

    * cgs2022: six clusters (720 orderings) and simultaneous orthogonalization
    use ../data/cn_cgs2022, clear
    tsset t
    quietly ds t date, not
    local v `r(varlist)'
    local cc `""CHF: X2YCHF X5YCHF X10YCHF" "DKK: X2YDKK X5YDKK X10YDKK" "EUR: X2YEUR X5YEUR X10YEUR" "GBP: X2YGBP X5YGBP X10YGBP" "NOK: X2YNOK X5YNOK X10YNOK" "SEK: X2YSEK X5YSEK X10YSEK""'
    connectedness `v', lags(1) horizon(10) clusters(`cc') `opt'
    cn_clcheck, case(cgs2022_diag) tol(1e-8)
    cn_cl095, case(cgs2022_diag) tol(1e-8)
    connectedness `v', lags(1) horizon(10) clusters(`cc') clorder(simultaneous) `opt'
    cn_clcheck, case(cgs2022_sco) tol(1e-8)

    * rolling dy2012: SP500 as the control cluster, four singletons, simultaneous
    use ../data/cn_dy2012, clear
    tsset t
    local v SP500 R_10Y DJUBSCOM USDX
    connectedness `v' in 1/400, lags(2) horizon(10) window(200) clcontrols(SP500) `opt'
    cn_clcheck, case(dy2012_ctl_roll) tol(1e-8)
    cn_cl095, case(dy2012_ctl_roll) tol(1e-8)
    quietly estat clorder
    tempname S
    matrix `S' = r(time)
    frame `e(frame)': quietly summarize cl_dfl_2, meanonly
    cn_assert "estat clorder, rolling: mean over time" `=abs(`S'[2, 1] - r(mean))' 1e-10
    connectedness `v' in 1/400, lags(2) horizon(10) window(200) clusters("a: SP500" "b: R_10Y" ///
        "c: DJUBSCOM" "d: USDX") `opt'
    cn_clcheck, case(dy2012_four_roll) tol(1e-8)
    cn_cl095, case(dy2012_four_roll) tol(1e-8)
    connectedness `v' in 1/400, lags(2) horizon(10) window(200) clusters("A: SP500 R_10Y" ///
        "B: DJUBSCOM USDX") clorder(simultaneous) `opt'
    cn_clcheck, case(dy2012_sco_roll) tol(1e-8)

    * TVP-VAR acg2020: simultaneous; GBP as the control cluster
    use ../data/cn_acg2020, clear
    tsset month
    connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) ///
        clusters("Europe: EUR CHF GBP" "Japan: JPY") clorder(simultaneous) `opt'
    cn_clcheck, case(acg2020_sco) tol(1e-6)
    connectedness EUR GBP JPY CHF, model(tvpvar) horizon(12) window(200) kappa2(0.96) ///
        clcontrols(GBP) clusters("Europe: EUR CHF" "Japan: JPY") `opt'
    cn_clcheck, case(acg2020_ctl) tol(1e-6)
    cn_cl095, case(acg2020_ctl) tol(1e-6)
}
else di as txt "  R 0.9.5: no fixture, skipped"

* identities of simultaneous orthogonalization: one cluster = generalized; invariance
* to the order of the variables and of the clusters
use ../data/cn_dy2012, clear
tsset t
local v SP500 R_10Y DJUBSCOM USDX
tempname G F1 F2 F3
connectedness `v', lags(4) horizon(10) `opt'
matrix `G' = e(fevd)
connectedness `v', lags(4) horizon(10) clusters("all: `v'") clorder(simultaneous) `opt'
matrix `F1' = e(fevd)
cn_assert "simultaneous, one cluster = generalized" `=mreldif(`F1', `G')' 1e-10
connectedness `v', lags(4) horizon(10) clusters("A: SP500 DJUBSCOM" "B: R_10Y USDX") clorder(simultaneous) `opt'
matrix `F1' = e(fevd)
connectedness USDX DJUBSCOM R_10Y SP500, lags(4) horizon(10) clusters("B: USDX R_10Y" "A: DJUBSCOM SP500") ///
    clorder(simultaneous) `opt'
matrix `F2' = e(fevd)
mata: cn__p = (4, 3, 2, 1); st_local("d1", strofreal(mreldif(st_matrix("`F2'")[cn__p, cn__p], st_matrix("`F1'")), "%12.0g"))
mata: mata drop cn__p
cn_assert "simultaneous: invariant to variable and cluster order" `d1' 1e-10
* control cluster averaged over two orderings = mean of the two given orderings
connectedness `v', lags(4) horizon(10) clcontrols(SP500) clusters("A: R_10Y" "B: DJUBSCOM USDX") `opt'
matrix `F1' = e(fevd)
connectedness `v', lags(4) horizon(10) clcontrols(SP500) clusters("A: R_10Y" "B: DJUBSCOM USDX") clorder(given) `opt'
matrix `F2' = e(fevd)
connectedness `v', lags(4) horizon(10) clcontrols(SP500) clusters("B: DJUBSCOM USDX" "A: R_10Y") clorder(given) `opt'
matrix `F3' = e(fevd)
matrix `F2' = (`F2' + `F3')/2
cn_assert "control cluster: average = mean of the given orderings" `=mreldif(`F1', `F2')' 1e-12
* bootstrap point estimates with the control cluster and with simultaneous
* orthogonalization equal the results of connectedness
connectedness `v' in 1/800, lags(2) horizon(10) clcontrols(SP500) clusters("A: R_10Y" "B: DJUBSCOM USDX") `opt'
local t1 = e(TCI)
tempname BT
capture noisily estat bootstrap, reps(5) seed(3) nodots
local d1 = 1
if !_rc {
    matrix `BT' = r(table)
    local d1 = abs(`BT'[1, 1] - `t1')
}
connectedness `v' in 1/800, lags(2) horizon(10) clusters("A: SP500 R_10Y" "B: DJUBSCOM USDX") clorder(simultaneous) `opt'
local t1 = e(TCI)
capture noisily estat bootstrap, reps(5) seed(3) nodots
local d2 = 1
if !_rc {
    matrix `BT' = r(table)
    local d2 = abs(`BT'[1, 1] - `t1')
}
cn_assert "bootstrap estimate, control cluster and simultaneous" `=`d1' + `d2'' 1e-10
* asymmetric with a control cluster: orderings frame of the overall system only
use ../data/cn_g2020, clear
tsset t
connectedness EUR GBP CHF JPY in 1/600, asymmetric clcontrols(EUR) clusters("A: GBP CHF" "B: JPY") `opt'
local fo `e(frame_clord)'
mata: st_local("ok", strofreal(st_frameexists("`fo'") & !st_frameexists("`fo'_pos")))
quietly estat clorder
cn_assert "asymmetric with clcontrols(): orderings frame" `=1 - `ok' + (r(norderings) != 2)' 0.5
connectedness EUR GBP CHF JPY in 1/600, clusters("A: EUR GBP" "B: CHF JPY") `opt'
cn_expect 198 estat clorder, digits(9)
connectedness EUR GBP CHF JPY in 1/600, `opt'
cn_expect 198 estat clorder

di _n as txt "{hline 78}"
quietly timer list
di as txt "Timing: rolling VAR, k = 19, 632 windows: " as res %6.1f r(t1) as txt " s"
di as txt "Timing: TVP-VAR, k = 11, T = 3,915:        " as res %6.1f r(t2) as txt " s"
di as txt "Timing: TVP-VAR, k = 18, T = 5,558:        " as res %6.1f r(t3) as txt " s"
di as txt "Timing: frequency, rolling VAR, H = 50:    " as res %6.1f r(t4) as txt " s"
di as txt "Timing: frequency, TVP-VAR, H = 100:       " as res %6.1f r(t5) as txt " s"
di as txt "Timing: extended joint, TVP-VAR, k = 12:   " as res %6.1f r(t6) as txt " s"
di as txt "Timing: QVAR, rolling, 501 windows:        " as res %6.1f r(t7) as txt " s"
di as txt "Timing: R2 (Kendall), rolling, 2,571 windows:" as res %6.1f r(t8) as txt " s"
di as txt "Timing: elastic net, 11 alphas, k = 4:     " as res %6.1f r(t9) as txt " s"
di as txt "Timing: LASSO, k = 51, static:             " as res %6.1f r(t10) as txt " s"
di as txt "Timing: DCC-GARCH case, k = 4, T = 4,200:  " as res %6.1f r(t11) as txt " s"
di as txt "Timing: portfolio (MCP), T = 4,200:         " as res %6.1f r(t12) as txt " s"
di as txt "Timing: DCC-GARCH R2, k = 4, T = 4,200:     " as res %6.1f r(t13) as txt " s"
di as txt "Timing: adaptive LASSO, k = 4, 2 lags:      " as res %6.1f r(t14) as txt " s"
di as txt "Timing: clusters, k = 18, 720 orderings:   " as res %6.1f r(t15) as txt " s"
di as txt "Timing: bootstrap, static, 80 draws:        " as res %6.1f r(t16) as txt " s"
di as txt "Timing: bootstrap, 4 events, k = 19:        " as res %6.1f r(t17) as txt " s"
di as txt "Timing: bootstrap of SAM, static, 80 draws: " as res %6.1f r(t18) as txt " s"
di as txt "Timing: clusters, k = 19, 5,040 orderings: " as res %6.1f r(t19) as txt " s"
di as txt "Tests run: " as res $cn_tests as txt "   failures: " as res $cn_fails
assert $cn_fails == 0
di as txt "ALL TESTS PASSED"
