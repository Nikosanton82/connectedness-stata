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
    syntax , case(string) tol(real) [prefix(string) MATrix(string) FRame(string)]
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
cn_expect 198 connectedness SP500 R_10Y, model(qvar) `opt'
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

di _n as txt "{hline 78}"
quietly timer list
di as txt "Timing: rolling VAR, k = 19, 632 windows: " as res %6.1f r(t1) as txt " s"
di as txt "Timing: TVP-VAR, k = 11, T = 3,915:        " as res %6.1f r(t2) as txt " s"
di as txt "Timing: TVP-VAR, k = 18, T = 5,558:        " as res %6.1f r(t3) as txt " s"
di as txt "Tests run: " as res $cn_tests as txt "   failures: " as res $cn_fails
assert $cn_fails == 0
di as txt "ALL TESTS PASSED"
