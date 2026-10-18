*! version 0.2.1  18oct2026
*! connectedness: VAR and TVP-VAR connectedness measures
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
        ASYMmetric DATEvar(varname numeric) noTABle DIGits(integer 2) ]

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
    else if inlist("`model'", "qvar", "lad", "lasso", "ridge", "elastic", "elasticnet", "dccgarch", "dcc-garch") {
        di as err "model(`model') is planned for a later release of {cmd:connectedness}"
        exit 198
    }
    else {
        di as err "model(`model') not allowed; specify {cmd:var} or {cmd:tvpvar}"
        exit 198
    }

    local prior = lower(strtrim(`"`prior'"'))
    if "`model'" == "var" & "`prior'" != "" {
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

    if `lags' < 1 {
        di as err "lags() must be a positive integer"
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
    if `nreg' <= `k'*`lags' + 1 {
        if "`model'" == "var" {
            di as err "too few observations per VAR estimation: `nreg' observations for " ///
                `k'*`lags' + 1 " coefficients per equation"
        }
        else {
            di as err "too few observations for the Bayes prior: `nreg' observations for " ///
                `k'*`lags' + 1 " coefficients per equation; increase window()"
        }
        exit 2001
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
            "`tab`i''", "`fev`i''", "`sc`i''", "`names'", "`datevar'")
        capture frame drop `f'
        frame create `f'
        frame `f' {
            mata: cn_postframe("`timevar'", "`tfmt'", "`names'", `savef', "`datevar'", "`dfmt'")
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
    if "`model'" == "tvpvar" {
        ereturn scalar kappa1 = `kappa1'
        ereturn scalar kappa2 = `kappa2'
        if "`prior'" == "minnesota" ereturn scalar gamma = `gamma'
        ereturn local prior "`prior'"
    }
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
    ereturn local approach "time"
    ereturn local model    "`model'"
    ereturn local title    "Connectedness approach"
    ereturn local estat_cmd "connectedness_estat"
    ereturn local cmd      "connectedness"

    Display, digits(`digits') `table'
end

program define Display
    syntax [, DIGits(integer 2) PARt(string) noTABle]
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
    local w = `digits' + 6
    if "`e(model)'" == "tvpvar" local mdesc "TVP-VAR(`e(lags)'), `e(prior)' prior"
    else if e(window) > 0      local mdesc "VAR(`e(lags)'), rolling window"
    else                       local mdesc "VAR(`e(lags)'), static"

    di _n as txt "`e(title)'" _col(49) "Number of variables = " as res %8.0g e(k)
    di as txt "Model:  " as res "`mdesc'" _col(49) as txt "Forecast horizon    = " as res %8.0g e(horizon)
    di as txt "FEVD:   " as res "`e(fevdtype)'" _col(49) as txt "Observations used   = " as res %8.0g e(N)
    di as txt "Sample: " as res "`e(first)' to `e(last)'" _col(49) as txt "Periods measured    = " as res %8.0g e(T0)
    if "`e(model)'" == "tvpvar" {
        local wtxt = cond(e(window) > 0, string(e(window)), "full sample")
        if "`e(prior)'" == "bayes" di as txt "Prior sample: " as res "`wtxt'" _c
        if "`e(prior)'" == "minnesota" di as txt "Prior tightness: " as res %4.2f e(gamma) _c
        di _col(49) as txt "kappa1, kappa2      = " as res %4.2f e(kappa1) as txt ", " as res %4.2f e(kappa2)
    }
    else if e(window) > 0 {
        di as txt "Window: " as res e(window) as txt " observations"
    }
    local fr "`e(frame)'"
    if "`suf'" != "" local fr "`e(frame`suf')'"
    di as txt "Results frame: " as res "`fr'" as txt " (first period " as res "`e(firstm)'" as txt ")"
    if "`e(asymmetric)'" != "" {
        local ptxt = cond("`suf'" == "_pos", "positive", cond("`suf'" == "_neg", "negative", "all"))
        di as txt "Asymmetric connectedness: " as res "`ptxt'" as txt " changes shown; " ///
            "see {cmd:connectedness, part()}"
    }
    di
    if "`table'" == "" {
        matlist e(table`suf'), format(%`w'.`digits'f) border(rows) lines(oneline) nodotz ///
            rowtitle(" ") twidth(12)
        di as txt "Rows: variance of the row variable explained by the column variable (percent)."
    }
    di as txt "TCI = " as res %6.`digits'f e(TCI`suf') as txt "    corrected TCI = " ///
        as res %6.`digits'f e(cTCI`suf') _c
    if "`e(corrected)'" != "" di as txt "    (NET row reports the corrected TCI)"
    else di
    if "`e(asymmetric)'" != "" & "`suf'" == "" {
        di as txt "Positive changes: TCI = " as res %6.`digits'f e(TCI_pos) ///
            as txt "    Negative changes: TCI = " as res %6.`digits'f e(TCI_neg)
    }
end
