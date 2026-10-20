*! version 0.5.0  20oct2026
*! connectedness: VAR, QVAR, LAD-VAR and TVP-VAR connectedness measures (time, frequency, joint)
*! Port of the ConnectednessApproach R package (D. Gabauer, version 1.0.4, GPL-3)
*! Author: Nikolaos Antonakakis.  License: GPL-3
program define connectedness, eclass
    version 16.1
    if replay() {
        if "`e(cmd)'" != "connectedness" error 301
        Display `0'
        exit
    }
    Estimate `0'
end

program define Estimate, eclass
    syntax varlist(numeric ts min=2) [if] [in] [,               ///
        Model(string) Lags(integer 1) Horizon(integer 10)         ///
        Window(integer 0) ORThogonalized GENeralized CORRected    ///
        KAPPA1(real 0.99) KAPPA2(real 0.99) PRior(string)         ///
        GAMma(real 0.01) FRame(string) REPLACE SAVEFevd noFEVD    ///
        ASYMmetric DATEvar(varname numeric) noTABle DIGits(integer 2) ///
        APProach(string) BANDs(numlist >1 sort) PARTition(numlist) SCENario(string) ///
        QUANTile(numlist >0 <1) CORMethod(string) ]

    // ---------------------------------------------------------- time series
    quietly tsset
    local timevar  `r(timevar)'
    local panelvar `r(panelvar)'
    local delta    = r(tdelta)
    if "`panelvar'" != "" {
        di as err "panel data are not supported; {cmd:tsset} the data with a time variable only"
        exit 459
    }

    marksample touse
    tsrevar `varlist'
    local tsvars `r(varlist)'
    markout `touse' `tsvars'

    quietly count if `touse'
    local T = r(N)
    if `T' == 0 error 2000
    quietly summarize `timevar' if `touse', meanonly
    if round((r(max) - r(min))/`delta') + 1 != `T' {
        di as err "the estimation sample must be contiguous in time;"
        di as err "gaps or missing values were found within the sample"
        exit 498
    }
    local k : word count `varlist'

    // ---------------------------------------------------------- options
    local model = lower(strtrim(`"`model'"'))
    if inlist("`model'", "", "var") {
        local model var
    }
    else if inlist("`model'", "tvpvar", "tvp-var", "tvp") {
        local model tvpvar
    }
    else if inlist("`model'", "qvar", "quantile") {
        local model qvar
    }
    else if inlist("`model'", "lad", "ladvar") {
        local model lad
    }
    else if inlist("`model'", "lasso", "ridge", "elastic", "elasticnet", "dccgarch", "dcc-garch") {
        di as err "model(`model') is planned for a later release of {cmd:connectedness}"
        exit 198
    }
    else {
        di as err "model(`model') not allowed; specify {cmd:var}, {cmd:tvpvar}, {cmd:qvar} or {cmd:lad}"
        exit 198
    }

    local approach = lower(strtrim(`"`approach'"'))
    if inlist("`approach'", "", "time") {
        local approach time
    }
    else if inlist("`approach'", "frequency", "freq") {
        local approach frequency
    }
    else if "`approach'" == "joint" {
        local approach joint
    }
    else if inlist("`approach'", "extjoint", "extendedjoint", "extended joint", "ejoint") {
        local approach extjoint
    }
    else if inlist("`approach'", "r2", "r2decomposed") {
        local approach r2
    }
    else if inlist("`approach'", "modelfree", "model-free", "unconditional") {
        local approach modelfree
    }
    else {
        di as err "approach(`approach') not allowed; specify {cmd:time}, {cmd:frequency}, " ///
            "{cmd:joint}, {cmd:extjoint}, {cmd:r2} or {cmd:modelfree}"
        exit 198
    }
    if "`approach'" != "frequency" & `"`bands'`partition'`scenario'"' != "" {
        di as err "options bands(), partition() and scenario() require approach(frequency)"
        exit 198
    }
    // R2 and model-free connectedness work on the data directly
    local cormethod = lower(strtrim("`cormethod'"))
    if !inlist("`approach'", "r2", "modelfree") & "`cormethod'" != "" {
        di as err "option cormethod() requires approach(r2) or approach(modelfree)"
        exit 198
    }
    if inlist("`approach'", "r2", "modelfree") {
        if "`cormethod'" == "" local cormethod pearson
        if !inlist("`cormethod'", "pearson", "spearman", "kendall") {
            di as err "cormethod() must be pearson, spearman or kendall"
            exit 198
        }
        if "`model'" != "var" {
            di as err "approach(`approach') does not use an estimated model; model() not allowed"
            exit 198
        }
        if "`orthogonalized'" != "" | "`asymmetric'" != "" {
            di as err "options orthogonalized and asymmetric are not allowed with approach(`approach')"
            exit 198
        }
        if "`approach'" == "modelfree" local lags 0
    }
    if inlist("`approach'", "joint", "extjoint") {
        if "`orthogonalized'" != "" {
            di as err "approach(`approach') is based on the generalized FEVD; " ///
                "option orthogonalized is not allowed"
            exit 198
        }
        if "`corrected'" != "" {
            di as err "option corrected is not available with approach(`approach')"
            exit 198
        }
    }

    // frequency bands: partition in radians, in decreasing order (R's partition)
    local partstr
    local blabels
    local nbands 0
    local scenario = lower(strtrim("`scenario'"))
    if "`approach'" == "frequency" {
        if "`bands'" != "" & "`partition'" != "" {
            di as err "options bands() and partition() may not be combined"
            exit 184
        }
        if `horizon' < 2 {
            di as err "approach(frequency) requires horizon() of at least 2"
            exit 198
        }
        if "`scenario'" == "" local scenario abs
        if !inlist("`scenario'", "abs", "wth") {
            di as err "scenario() must be abs or wth"
            exit 198
        }
        if "`bands'" != "" {
            local partstr : display %21.17g c(pi) + 0.00001
            foreach b of local bands {
                local pv : display %21.17g c(pi)/`b'
                local partstr `partstr' `pv'
            }
            local partstr `partstr' 0
        }
        else if "`partition'" != "" {
            local prev .
            foreach v of local partition {
                if `v' < 0 | `v' >= `prev' {
                    di as err "partition() must list the band limits in radians in decreasing " ///
                        "order, ending at zero or above"
                    exit 198
                }
                local prev `v'
                local pv : display %21.17g `v'
                local partstr `partstr' `pv'
            }
            if `: word count `partstr'' < 2 {
                di as err "partition() needs at least two limits"
                exit 198
            }
        }
        else {
            local p1 : display %21.17g c(pi)
            local p2 : display %21.17g c(pi)/2
            local partstr `p1' `p2' 0
        }
        local nbands = `: word count `partstr'' - 1
        forvalues b = 1/`nbands' {
            local hi : word `b' of `partstr'
            local lo : word `=`b'+1' of `partstr'
            local phi = round(c(pi)/`hi')
            local plo = cond(`lo' == 0, "Inf", string(round(c(pi)/`lo')))
            local blabels `blabels' `phi'-`plo'
        }
    }

    // quantiles for model(qvar): one value, or one per variable
    local taus
    if "`quantile'" != "" & "`model'" != "qvar" {
        di as err "option quantile() is allowed only with model(qvar)"
        exit 198
    }
    if "`model'" == "qvar" {
        if "`quantile'" == "" local quantile 0.5
        local nq : word count `quantile'
        if `nq' != 1 & `nq' != `k' {
            di as err "quantile() must contain one value or one value per variable (`k')"
            exit 198
        }
        forvalues i = 1/`k' {
            local q : word `=cond(`nq' == 1, 1, `i')' of `quantile'
            local taus `taus' `q'
        }
    }
    if "`model'" != "tvpvar" & (`kappa1' != 0.99 | `kappa2' != 0.99 | `gamma' != 0.01) {
        di as err "options kappa1(), kappa2() and gamma() are allowed only with model(tvpvar)"
        exit 198
    }

    // components stored with suffixes: frequency bands, or R2 parts
    local csufs
    if "`approach'" == "frequency" {
        forvalues b = 1/`nbands' {
            local csufs `csufs' f`b'
        }
    }
    if "`approach'" == "r2" & `lags' > 0 {
        local csufs c l
        local blabels Contemporaneous Lagged
    }
    local ncomp : word count `csufs'

    local prior = lower(strtrim(`"`prior'"'))
    if "`model'" != "tvpvar" & "`prior'" != "" {
        di as err "option prior() is allowed only with model(tvpvar)"
        exit 198
    }
    if inlist("`prior'", "", "bayes", "bayesprior") {
        local prior bayes
    }
    else if inlist("`prior'", "minnesota", "minnesotaprior") {
        local prior minnesota
    }
    else if inlist("`prior'", "uninformative", "uninformativeprior", "flat") {
        local prior uninformative
    }
    else {
        di as err "prior(`prior') not allowed; specify {cmd:bayes}, {cmd:minnesota} or {cmd:uninformative}"
        exit 198
    }

    if "`orthogonalized'" != "" & "`generalized'" != "" {
        di as err "options orthogonalized and generalized may not be combined"
        exit 184
    }
    local gen  = ("`orthogonalized'" == "")
    local corr = ("`corrected'" != "")

    if `lags' < 1 & !inlist("`approach'", "r2", "modelfree") {
        di as err "lags() must be a positive integer"
        exit 198
    }
    if `lags' < 0 {
        di as err "lags() must be zero or a positive integer"
        exit 198
    }
    if `horizon' < 1 {
        di as err "horizon() must be a positive integer"
        exit 198
    }
    if `window' < 0 {
        di as err "window() must be zero (static or full-sample prior) or a positive integer"
        exit 198
    }
    if `window' > `T' {
        di as err "window(`window') exceeds the number of observations (`T')"
        exit 498
    }
    if `kappa1' <= 0 | `kappa1' >= 1 | `kappa2' <= 0 | `kappa2' >= 1 {
        di as err "kappa1() and kappa2() must lie strictly between 0 and 1"
        exit 198
    }
    if `gamma' <= 0 {
        di as err "gamma() must be positive"
        exit 198
    }
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }

    // observations per regression: rolling or prior sample minus lags
    local ws = cond(`window' == 0, `T', `window' - `lags')
    local nreg = `ws' - `lags'
    if "`model'" == "tvpvar" & "`prior'" != "bayes" {
        local nreg = `T' - `lags'
    }
    if inlist("`approach'", "r2", "modelfree") {
        // R2 decomposition: each window holds window - 2*lags rows of (x_t, ..., x_t-p)
        local nreg = cond(`window' == 0, `T' - `lags', `window' - 2*`lags')
        if `nreg' <= `k'*(`lags' + 1) {
            di as err "too few observations per window: `nreg' rows for " `k'*(`lags'+1) " series"
            exit 2001
        }
    }
    else if `nreg' <= `k'*`lags' + 1 {
        if "`model'" != "tvpvar" {
            di as err "too few observations per VAR estimation: `nreg' observations for " ///
                `k'*`lags' + 1 " coefficients per equation"
        }
        else {
            di as err "too few observations for the Bayes prior: `nreg' observations for " ///
                `k'*`lags' + 1 " coefficients per equation; increase window()"
        }
        exit 2001
    }

    // number of variables in the results frame
    local np = `k'*(`k'-1)/2
    local nvars = cond("`approach'" == "joint", 1 + 3*`k', 2 + 4*`k' + 3*`np')
    if "`fevd'" == "" local nvars = `nvars' + `k'*`k'
    local nvars = (`nvars')*(`ncomp' + 1) + 2
    if `nvars' > c(maxvar) - 10 {
        di as err "the results frame would need `nvars' variables, more than maxvar (" c(maxvar) ") allows;"
        di as err "specify {cmd:nofevd}, reduce the number of variables or frequency bands, or {help set maxvar}"
        exit 900
    }

    // ---------------------------------------------------------- results frame
    if `"`frame'"' == "" local frame cn_results
    confirm name `frame'
    local parts all
    local frames `frame'
    if "`asymmetric'" != "" {
        if length("`frame'") > 28 {
            di as err "frame() name must not exceed 28 characters with asymmetric"
            exit 198
        }
        local parts all positive negative
        local frames `frame' `frame'_pos `frame'_neg
    }
    foreach f of local frames {
        if "`f'" == "`c(frame)'" {
            di as err "frame `f' is the current frame; choose another name in frame()"
            exit 498
        }
        mata: st_local("fexist", strofreal(st_frameexists("`f'")))
        if `fexist' & "`replace'" == "" {
            di as err "frame {bf:`f'} already exists; specify {cmd:replace} or another {cmd:frame()}"
            exit 110
        }
    }
    if "`datevar'" != "" & "`datevar'" == "`timevar'" local datevar

    // ---------------------------------------------------------- estimation
    capture mata: mata describe cn_main()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }

    // positive and negative parts of each series (Adekoya et al. 2022)
    if "`asymmetric'" != "" {
        local posvars
        local negvars
        foreach v of local tsvars {
            tempvar p n
            quietly generate double `p' = cond(missing(`v'), ., max(`v', 0)) if `touse'
            quietly generate double `n' = cond(missing(`v'), ., min(`v', 0)) if `touse'
            local posvars `posvars' `p'
            local negvars `negvars' `n'
        }
    }

    local names `varlist'
    local tfmt : format `timevar'
    local dfmt
    if "`datevar'" != "" local dfmt : format `datevar'
    local savef = ("`fevd'" == "")
    local i 0
    foreach part of local parts {
        local ++i
        local f : word `i' of `frames'
        if "`part'" == "all"      local vars `tsvars'
        if "`part'" == "positive" local vars `posvars'
        if "`part'" == "negative" local vars `negvars'
        tempname tab`i' fev`i' sc`i'
        mata: cn_main("`vars'", "`touse'", "`timevar'", "`model'", `lags', `horizon', ///
            `window', `gen', `corr', "`prior'", `kappa1', `kappa2', `gamma',         ///
            "`tab`i''", "`fev`i''", "`sc`i''", "`names'", "`datevar'",               ///
            "`approach'", strtoreal(tokens("`partstr'")), "`scenario'", "`blabels'", ///
            strtoreal(tokens("`taus'")), "`csufs'", "`cormethod'")
        capture frame drop `f'
        frame create `f'
        frame `f' {
            mata: cn_postframe("`timevar'", "`tfmt'", `savef', "`datevar'", "`dfmt'")
            quietly tsset `timevar', delta(`delta')
            label data "connectedness results (`part')"
            char _dta[cn_kind]      "main"
            char _dta[cn_part]      "`part'"
            char _dta[cn_names]     "`names'"
            char _dta[cn_k]         "`k'"
            char _dta[cn_timevar]   "`timevar'"
            char _dta[cn_datevar]   "`datevar'"
            char _dta[cn_corrected] "`corrected'"
            char _dta[cn_fevd]      "`savef'"
            char _dta[cn_approach]  "`approach'"
            char _dta[cn_nbands]    "`nbands'"
            char _dta[cn_bands]     "`blabels'"
            char _dta[cn_comps]     "`csufs'"
        }
    }
    local sc `sc1'

    // sample period, in dates when datevar() is given
    local first  : display `tfmt' `sc'[1,5]
    local last   : display `tfmt' `sc'[1,6]
    local firstm : display `tfmt' `sc'[1,7]
    if "`datevar'" != "" {
        foreach w in first last firstm {
            local col = cond("`w'" == "first", 5, cond("`w'" == "last", 6, 7))
            quietly summarize `datevar' if `touse' & `timevar' == `sc'[1,`col'], meanonly
            local `w' : display `dfmt' r(min)
        }
    }

    // ---------------------------------------------------------- stored results
    ereturn post, esample(`touse')
    ereturn scalar N       = `sc'[1,3]
    ereturn scalar T0      = `sc'[1,4]
    ereturn scalar k       = `k'
    ereturn scalar lags    = `lags'
    ereturn scalar horizon = `horizon'
    ereturn scalar window  = `window'
    ereturn scalar TCI     = `sc'[1,1]
    ereturn scalar cTCI    = `sc'[1,2]
    if "`approach'" == "joint" ereturn scalar cTCI = .
    if "`model'" == "qvar" {
        local tl : list uniq taus
        if `: word count `tl'' == 1 ereturn scalar quantile = `tl'
        ereturn local quantiles "`taus'"
    }
    if "`model'" == "tvpvar" {
        ereturn scalar kappa1 = `kappa1'
        ereturn scalar kappa2 = `kappa2'
        if "`prior'" == "minnesota" ereturn scalar gamma = `gamma'
        ereturn local prior "`prior'"
    }
    if "`approach'" == "frequency" {
        forvalues b = 1/`nbands' {
            ereturn scalar TCI_f`b'  = `sc'[1, 6 + 2*`b']
            ereturn scalar cTCI_f`b' = `sc'[1, 7 + 2*`b']
            ereturn matrix table_f`b' = `tab1'f`b'
            ereturn matrix fevd_f`b'  = `fev1'f`b'
        }
        ereturn scalar nbands = `nbands'
        ereturn local bands     "`blabels'"
        ereturn local partition "`partstr'"
        ereturn local scenario  "`scenario'"
    }
    if "`approach'" == "r2" & `lags' > 0 {
        local j 0
        foreach c in c l {
            local ++j
            ereturn scalar TCI_`c'  = `sc'[1, 6 + 2*`j']
            ereturn scalar cTCI_`c' = `sc'[1, 7 + 2*`j']
            ereturn matrix table_`c' = `tab1'`c'
            ereturn matrix fevd_`c'  = `fev1'`c'
        }
    }
    if inlist("`approach'", "r2", "modelfree") ereturn local cormethod "`cormethod'"
    ereturn matrix table = `tab1'
    ereturn matrix fevd  = `fev1'
    if "`asymmetric'" != "" {
        ereturn scalar TCI_pos  = `sc2'[1,1]
        ereturn scalar cTCI_pos = `sc2'[1,2]
        ereturn scalar TCI_neg  = `sc3'[1,1]
        ereturn scalar cTCI_neg = `sc3'[1,2]
        ereturn matrix table_pos = `tab2'
        ereturn matrix fevd_pos  = `fev2'
        ereturn matrix table_neg = `tab3'
        ereturn matrix fevd_neg  = `fev3'
        ereturn local frame_pos "`frame'_pos"
        ereturn local frame_neg "`frame'_neg"
        ereturn local asymmetric "asymmetric"
    }
    ereturn local first    "`=strtrim("`first'")'"
    ereturn local last     "`=strtrim("`last'")'"
    ereturn local firstm   "`=strtrim("`firstm'")'"
    ereturn local datevar  "`datevar'"
    ereturn local timevar  "`timevar'"
    ereturn local frame    "`frame'"
    ereturn local names    "`names'"
    ereturn local depvar   "`varlist'"
    ereturn local corrected "`corrected'"
    ereturn local fevdtype = cond(`gen', "generalized", "orthogonalized")
    if "`approach'" == "r2" ereturn local fevdtype "R2 decomposition (`cormethod')"
    if "`approach'" == "modelfree" ereturn local fevdtype "squared correlations (`cormethod')"
    if inlist("`approach'", "r2", "modelfree") ereturn scalar horizon = .
    ereturn local approach "`approach'"
    ereturn local model    "`model'"
    ereturn local title    "Connectedness approach"
    ereturn local estat_cmd "connectedness_estat"
    ereturn local cmd      "connectedness"

    Display, digits(`digits') `table'
end

program define Display
    syntax [, DIGits(integer 2) PARt(string) BANd(integer 0) COMPonent(string) noTABle]
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }
    local part = lower("`part'")
    if inlist("`part'", "", "all") local suf
    else if inlist("`part'", "pos", "positive") local suf _pos
    else if inlist("`part'", "neg", "negative") local suf _neg
    else {
        di as err "part() must be all, positive or negative"
        exit 198
    }
    if "`suf'" != "" & "`e(asymmetric)'" == "" {
        di as err "part() requires results estimated with the asymmetric option"
        exit 198
    }
    local approach "`e(approach)'"
    if "`approach'" == "" local approach time
    if `band' != 0 {
        if "`approach'" != "frequency" {
            di as err "band() requires results estimated with approach(frequency)"
            exit 198
        }
        if `band' < 0 | `band' > e(nbands) {
            di as err "band() must be between 1 and " e(nbands)
            exit 198
        }
        if "`suf'" != "" {
            di as err "band() may not be combined with part()"
            exit 184
        }
        local suf _f`band'
    }
    local component = lower("`component'")
    if "`component'" != "" {
        if "`approach'" != "r2" | "`e(TCI_c)'" == "" {
            di as err "component() requires results estimated with approach(r2) and lags() > 0"
            exit 198
        }
        if "`suf'" != "" {
            di as err "component() may not be combined with part() or band()"
            exit 184
        }
        if inlist("`component'", "c", "cont", "contemporaneous") local suf _c
        else if inlist("`component'", "l", "lag", "lagged") local suf _l
        else {
            di as err "component() must be contemporaneous or lagged"
            exit 198
        }
    }
    local w = `digits' + 6
    local mname "VAR"
    if "`e(model)'" == "qvar" local mname "QVAR"
    if "`e(model)'" == "lad"  local mname "LAD-VAR"
    local wdesc = cond(e(window) > 0, "rolling window", "static")
    if "`approach'" == "r2" local mdesc "R2, `e(lags)' lag(s), `wdesc'"
    else if "`approach'" == "modelfree" local mdesc "model-free, `wdesc'"
    else if "`e(model)'" == "tvpvar" local mdesc "TVP-VAR(`e(lags)'), `e(prior)' prior"
    else if e(window) > 0      local mdesc "`mname'(`e(lags)'), rolling window"
    else                       local mdesc "`mname'(`e(lags)'), static"

    di _n as txt "`e(title)'" _col(49) "Number of variables = " as res %8.0g e(k)
    if inlist("`approach'", "r2", "modelfree") {
        di as txt "Model:  " as res "`mdesc'" _col(49) as txt "Lags                = " as res %8.0g e(lags)
    }
    else {
        di as txt "Model:  " as res "`mdesc'" _col(49) as txt "Forecast horizon    = " as res %8.0g e(horizon)
    }
    di as txt "FEVD:   " as res "`e(fevdtype)'" _col(49) as txt "Observations used   = " as res %8.0g e(N)
    di as txt "Sample: " as res "`e(first)' to `e(last)'" _col(49) as txt "Periods measured    = " as res %8.0g e(T0)
    if "`e(model)'" == "tvpvar" {
        local wtxt = cond(e(window) > 0, string(e(window)), "full sample")
        if "`e(prior)'" == "bayes" di as txt "Prior sample: " as res "`wtxt'" _c
        if "`e(prior)'" == "minnesota" di as txt "Prior tightness: " as res %4.2f e(gamma) _c
        di _col(49) as txt "kappa1, kappa2      = " as res %4.2f e(kappa1) as txt ", " as res %4.2f e(kappa2)
    }
    else if e(window) > 0 {
        di as txt "Window: " as res e(window) as txt " observations" _c
    }
    if "`e(model)'" == "qvar" {
        local qtxt = cond("`e(quantile)'" != "", strofreal(e(quantile)), "`e(quantiles)'")
        di _col(49) as txt "Quantile (tau)      = " as res %8s "`qtxt'"
    }
    else if "`e(model)'" != "tvpvar" & e(window) > 0 di
    if "`approach'" == "frequency" {
        local bl "`e(bands)'"
        local bl : subinstr local bl " " ", ", all
        di as txt "Approach: " as res "frequency" as txt " (bands in periods: " as res "`bl'" ///
            as txt "; " as res strupper("`e(scenario)'") as txt ")"
    }
    else if "`approach'" == "joint" {
        di as txt "Approach: " as res "joint" as txt " (Lastrapes and Wiesen 2021)"
    }
    else if "`approach'" == "extjoint" {
        di as txt "Approach: " as res "extended joint" as txt " (Balcilar, Gabauer and Umar 2021)"
    }
    else if "`approach'" == "r2" {
        di as txt "Approach: " as res "R2 decomposed connectedness" as txt " (" as res "`e(cormethod)'" ///
            as txt " correlations)"
    }
    else if "`approach'" == "modelfree" {
        di as txt "Approach: " as res "model-free connectedness" as txt " (" as res "`e(cormethod)'" ///
            as txt " correlations)"
    }
    local fr "`e(frame)'"
    if inlist("`suf'", "_pos", "_neg") local fr "`e(frame`suf')'"
    di as txt "Results frame: " as res "`fr'" as txt " (first period " as res "`e(firstm)'" as txt ")"
    if "`e(asymmetric)'" != "" {
        local ptxt = cond("`suf'" == "_pos", "positive", cond("`suf'" == "_neg", "negative", "all"))
        di as txt "Asymmetric connectedness: " as res "`ptxt'" as txt " changes shown; " ///
            "see {cmd:connectedness, part()}"
    }
    if `band' {
        di as txt "Frequency band: " as res "`: word `band' of `e(bands)''" as txt " periods"
    }
    else if inlist("`suf'", "_c", "_l") {
        di as txt "Component: " as res cond("`suf'" == "_c", "contemporaneous", "lagged")
    }
    else if "`approach'" == "r2" & "`e(TCI_c)'" != "" {
        di as txt "Table: " as res "overall" as txt "; see {cmd:connectedness, component()}"
    }
    else if "`approach'" == "frequency" {
        di as txt "Table: " as res "total over the bands" as txt "; see {cmd:connectedness, band()}"
    }
    di
    if "`table'" == "" {
        matlist e(table`suf'), format(%`w'.`digits'f) border(rows) lines(oneline) nodotz ///
            rowtitle(" ") twidth(12)
        if "`approach'" == "joint" {
            di as txt "Rows: generalized FEVD (percent); TO, FROM and NET are joint spillovers."
        }
        else if "`approach'" == "r2" {
            di as txt "Rows: share of the row variable's variance explained by the column variable (R2, percent)."
        }
        else {
            di as txt "Rows: variance of the row variable explained by the column variable (percent)."
        }
    }
    if "`approach'" == "joint" {
        di as txt "Joint TCI = " as res %6.`digits'f e(TCI)
        exit
    }
    di as txt "TCI = " as res %6.`digits'f e(TCI`suf') as txt "    corrected TCI = " ///
        as res %6.`digits'f e(cTCI`suf') _c
    if "`e(corrected)'" != "" di as txt "    (NET row reports the corrected TCI)"
    else di
    if "`e(asymmetric)'" != "" & "`suf'" == "" {
        di as txt "Positive changes: TCI = " as res %6.`digits'f e(TCI_pos) ///
            as txt "    Negative changes: TCI = " as res %6.`digits'f e(TCI_neg)
    }
    if "`approach'" == "r2" & "`e(TCI_c)'" != "" & "`suf'" == "" {
        di _n as txt %-16s "Component" %12s "TCI" %16s "corrected TCI"
        di as txt "{hline 44}"
        di as txt %-16s "Contemporaneous" as res %12.`digits'f e(TCI_c) %16.`digits'f e(cTCI_c)
        di as txt %-16s "Lagged" as res %12.`digits'f e(TCI_l) %16.`digits'f e(cTCI_l)
        di as txt "{hline 44}"
        di as txt %-16s "Overall" as res %12.`digits'f e(TCI) %16.`digits'f e(cTCI)
    }
    if "`approach'" == "frequency" & `band' == 0 {
        di _n as txt %-16s "Band (periods)" %12s "TCI" %16s "corrected TCI"
        di as txt "{hline 44}"
        forvalues b = 1/`=e(nbands)' {
            di as txt %-16s "`: word `b' of `e(bands)''" as res %12.`digits'f e(TCI_f`b') ///
                %16.`digits'f e(cTCI_f`b')
        }
        di as txt "{hline 44}"
        di as txt %-16s "Total" as res %12.`digits'f e(TCI) %16.`digits'f e(cTCI)
    }
end
