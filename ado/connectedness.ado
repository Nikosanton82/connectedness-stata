*! version 0.9.3  25oct2026
*! connectedness: VAR-type, TVP-VAR and DCC-GARCH connectedness measures (time, frequency, joint, R2)
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
        QUANTile(numlist >0 <1) CORMethod(string)                 ///
        ALPha(numlist max=1 >=0 <=1) NALPha(integer 10) NFOLDs(integer 10) ///
        LOSS(string) FOLDs(string) SEED(string)                   ///
        STANDardize DCCParams(name) noCOVariance                  ///
        NORMalize(string) ADAPTive                                ///
        POSitive(varlist numeric ts) NEGative(varlist numeric ts) ///
        CLUSters(string asis) CLORDer(string) ]

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
    local k : word count `varlist'

    // asymmetric connectedness of user-supplied positive and negative series, such as
    // realized semivariances (Barunik, Kocenda and Vacha 2016)
    local semi = ("`positive'`negative'" != "")
    if `semi' {
        if "`positive'" == "" | "`negative'" == "" {
            di as err "options positive() and negative() must be specified together"
            exit 198
        }
        if "`asymmetric'" != "" {
            di as err "option asymmetric may not be combined with positive() and negative()"
            exit 184
        }
        if `: word count `positive'' != `k' | `: word count `negative'' != `k' {
            di as err "positive() and negative() must each contain `k' variables, one per variable in varlist"
            exit 198
        }
        tsrevar `positive'
        local posvars `r(varlist)'
        tsrevar `negative'
        local negvars `r(varlist)'
        markout `touse' `posvars' `negvars'
        local asymmetric asymmetric
    }

    quietly count if `touse'
    local T = r(N)
    if `T' == 0 error 2000
    quietly summarize `timevar' if `touse', meanonly
    if round((r(max) - r(min))/`delta') + 1 != `T' {
        di as err "the estimation sample must be contiguous in time;"
        di as err "gaps or missing values were found within the sample"
        exit 498
    }

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
    else if inlist("`model'", "lasso", "ridge") {
        local model `model'
    }
    else if inlist("`model'", "elastic", "elasticnet", "enet") {
        local model elasticnet
    }
    else if inlist("`model'", "dccgarch", "dcc-garch", "dcc") {
        local model dccgarch
    }
    else {
        di as err "model(`model') not allowed; specify {cmd:var}, {cmd:tvpvar}, {cmd:qvar}, " ///
            "{cmd:lad}, {cmd:lasso}, {cmd:ridge}, {cmd:elasticnet} or {cmd:dccgarch}"
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
    // R2 connectedness of DCC-GARCH conditional correlations (Cocca et al. 2024)
    local dccr2 = ("`model'" == "dccgarch" & "`approach'" == "r2")
    if `dccr2' {
        if !inlist("`cormethod'", "", "pearson") {
            di as err "option cormethod() is not allowed with model(dccgarch); " ///
                "the conditional correlations are used"
            exit 198
        }
        local cormethod dcc
        local lags 0
    }
    if inlist("`approach'", "r2", "modelfree") & !`dccr2' {
        if "`cormethod'" == "" local cormethod pearson
        if !inlist("`cormethod'", "pearson", "spearman", "kendall") {
            di as err "cormethod() must be pearson, spearman or kendall"
            exit 198
        }
        if "`model'" != "var" {
            di as err "approach(`approach') does not use an estimated model; " ///
                "model() not allowed (except model(dccgarch) with approach(r2))"
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

    // DCC-GARCH: full-sample estimation, time-domain volatility FEVD (as in R)
    local dcc = ("`model'" == "dccgarch")
    if !`dcc' & ("`standardize'" != "" | "`dccparams'" != "") {
        di as err "options standardize and dccparams() are allowed only with model(dccgarch)"
        exit 198
    }
    if `dcc' {
        if !inlist("`approach'", "time", "r2") {
            di as err "model(dccgarch) is available with approach(time) or approach(r2) only"
            exit 198
        }
        if `dccr2' & "`standardize'" != "" {
            di as err "option standardize is not allowed with approach(r2)"
            exit 198
        }
        if `window' != 0 {
            di as err "model(dccgarch) is estimated on the full sample; option window() not allowed"
            exit 198
        }
        if "`orthogonalized'" != "" | "`asymmetric'" != "" {
            di as err "options orthogonalized and asymmetric are not allowed with model(dccgarch)"
            exit 198
        }
        if `T' < 100 {
            di as err "model(dccgarch) requires at least 100 observations, as rmgarch does"
            exit 2001
        }
        if "`dccparams'" != "" confirm matrix `dccparams'
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
    // elastic net, LASSO and ridge VARs (R ElasticNetVAR)
    local enet = inlist("`model'", "lasso", "ridge", "elasticnet")
    local loss = lower(strtrim("`loss'"))
    local folds = lower(strtrim("`folds'"))
    if !`enet' & ("`alpha'`loss'`folds'`seed'" != "" | `nalpha' != 10 | `nfolds' != 10) {
        di as err "options alpha(), nalpha(), nfolds(), loss(), folds() and seed() are allowed " ///
            "only with model(lasso), model(ridge) or model(elasticnet)"
        exit 198
    }
    local alphas
    if `enet' {
        if "`loss'" == "" local loss mae
        if !inlist("`loss'", "mae", "mse") {
            di as err "loss() must be mae or mse"
            exit 198
        }
        if "`folds'" == "" local folds random
        if !inlist("`folds'", "random", "interleaved") {
            di as err "folds() must be random or interleaved"
            exit 198
        }
        if `nfolds' < 3 {
            di as err "nfolds() must be at least 3"
            exit 198
        }
        if "`alpha'" != "" & "`model'" != "elasticnet" {
            di as err "option alpha() is allowed only with model(elasticnet)"
            exit 198
        }
        if `nalpha' < 1 {
            di as err "nalpha() must be a positive integer"
            exit 198
        }
        if "`model'" == "lasso" local alphas 1
        else if "`model'" == "ridge" local alphas 0
        else if "`alpha'" != "" local alphas `alpha'
        else {
            // R: seq(0, 1, 1/n_alpha), that is j*(1/n_alpha), at full precision
            forvalues j = 0/`nalpha' {
                local a : display %21.17g `j'*(1/`nalpha')
                local alphas `alphas' `=strtrim("`a'")'
            }
        }
        if "`seed'" != "" set seed `seed'
    }
    if "`adaptive'" != "" & !`enet' {
        di as err "option adaptive is allowed only with model(lasso), model(ridge) or model(elasticnet)"
        exit 198
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

    // normalization of the generalized FEVD: row (Diebold and Yilmaz 2012) or scalar
    // (Caloia, Cipollini and Muzzioli 2019)
    local normalize = lower(strtrim("`normalize'"))
    local normalize : subinstr local normalize " " "", all
    if inlist("`normalize'", "", "row") local normalize row
    else if inlist("`normalize'", "maxrow", "rowmax", "maxrowsum") local normalize maxrow
    else if inlist("`normalize'", "maxcol", "colmax", "maxcolumn", "maxcolsum", "maxcolumnsum") ///
        local normalize maxcol
    else if inlist("`normalize'", "spectral", "spectralradius") local normalize spectral
    else {
        di as err "normalize() must be row, maxrow, maxcol or spectral"
        exit 198
    }
    if "`normalize'" != "row" {
        if "`approach'" != "time" | `dcc' {
            di as err "normalize(`normalize') requires approach(time) and a VAR-type or TVP-VAR model"
            exit 198
        }
        if !`gen' {
            di as err "normalize(`normalize') applies to the generalized FEVD; " ///
                "option orthogonalized is not allowed"
            exit 198
        }
        local gen = 1 + ("`normalize'" == "maxrow") + 2*("`normalize'" == "maxcol") + ///
            3*("`normalize'" == "spectral")
    }

    // clustered identification (Buchwalter, Diebold and Yilmaz 2026a): clusters listed
    // as in estat group, "name: var var ..."; variables not listed form clusters of
    // their own
    local clorder = lower(strtrim("`clorder'"))
    local ncl 0
    if `"`clusters'"' == "" & "`clorder'" != "" {
        di as err "option clorder() requires clusters()"
        exit 198
    }
    if `"`clusters'"' != "" {
        if "`approach'" != "time" | `dcc' {
            di as err "clusters() requires approach(time) and a VAR-type or TVP-VAR model"
            exit 198
        }
        if !`gen' {
            di as err "clusters() defines its own identification; option orthogonalized is not allowed"
            exit 198
        }
        if "`normalize'" != "row" {
            di as err "clusters() may not be combined with normalize(`normalize')"
            exit 198
        }
        if "`clorder'" == "" local clorder average
        if !inlist("`clorder'", "average", "given") {
            di as err "clorder() must be average or given"
            exit 198
        }
        tempname CLID
        matrix `CLID' = J(`k', 1, 0)
        local clnames
        local m : word count `clusters'
        forvalues g = 1/`m' {
            local spec : word `g' of `clusters'
            local gname "C`g'"
            local pos = strpos(`"`spec'"', ":")
            if `pos' {
                local gname = strtrim(substr(`"`spec'"', 1, `pos' - 1))
                local spec  = substr(`"`spec'"', `pos' + 1, .)
                if "`gname'" == "" local gname "C`g'"
            }
            local gname = strtoname("`gname'")
            if `: list gname in clnames' {
                di as err "cluster name `gname' is used twice"
                exit 198
            }
            local clnames `clnames' `gname'
            local nv 0
            foreach v of local spec {
                local p : list posof "`v'" in varlist
                if `p' == 0 {
                    di as err "`v' in clusters() is not one of the variables: `varlist'"
                    exit 198
                }
                if `CLID'[`p', 1] != 0 {
                    di as err "`v' appears in more than one cluster"
                    exit 198
                }
                matrix `CLID'[`p', 1] = `g'
                local ++nv
            }
            if `nv' == 0 {
                di as err "cluster `g' is empty"
                exit 198
            }
        }
        // variables not listed: one cluster each, in varlist order
        local ncl `m'
        forvalues i = 1/`k' {
            if `CLID'[`i', 1] == 0 {
                local ++ncl
                matrix `CLID'[`i', 1] = `ncl'
                local vn : word `i' of `varlist'
                local vn = strtoname("`vn'")
                if `: list vn in clnames' local vn "C`ncl'"
                local clnames `clnames' `vn'
            }
        }
        if "`clorder'" == "average" & `ncl' > 6 {
            di as err "averaging over all orderings is limited to 6 clusters (720 orderings); " ///
                "specify clorder(given) or fewer clusters"
            exit 198
        }
        local nord = cond("`clorder'" == "average", exp(lnfactorial(`ncl')), 1)
        local nord = round(`nord')
        local gen 5
    }

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
    else if `dcc' {
        // no VAR: the lags() option is not used
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
    local savecov = ("`covariance'" == "" & (!inlist("`approach'", "r2", "modelfree") | `dcc'))
    local nvars = (`nvars')*(`ncomp' + 1) + 2
    if "`asymmetric'" != "" local nvars = `nvars' + 1 + 3*`k'
    if `ncl' local nvars = `nvars' + (5*`ncl' + 3 + `ncl'*(`ncl'-1))*(1 + 2*("`asymmetric'" != ""))
    if `savecov' {
        // the covariances are left out, with a note, when they would not fit
        if `nvars' + `k'*(`k'+1)/2 > c(maxvar) - 10 {
            local savecov 0
            local covnote 1
        }
        else local nvars = `nvars' + `k'*(`k'+1)/2
    }
    if `nvars' > c(maxvar) - 10 {
        di as err "the results frame would need `nvars' variables, more than maxvar (" c(maxvar) ") allows;"
        di as err "specify {cmd:nofevd} or {cmd:nocovariance}, reduce the number of variables or frequency bands, or {help set maxvar}"
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
    // without datevar(), a numeric variable named date with a date format is used
    // when the time variable has no date format (for example an observation index
    // t); datevar() set to the time variable itself turns this off
    local dvauto = ("`datevar'" == "")
    if "`datevar'" != "" & "`datevar'" == "`timevar'" local datevar
    local tfmt0 : format `timevar'
    if `dvauto' & substr("`tfmt0'", 1, 2) != "%t" {
        capture confirm numeric variable date, exact
        if !_rc {
            local dfmt0 : format date
            if substr("`dfmt0'", 1, 2) == "%t" local datevar date
        }
    }

    // ---------------------------------------------------------- estimation
    capture mata: mata describe cn_main()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }

    // positive and negative parts of each series (Adekoya et al. 2022)
    if "`asymmetric'" != "" & !`semi' {
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
    if `ncl' {
        mata: cn__clid = st_matrix("`CLID'")
        if "`clorder'" == "average" mata: cn__clperm = cn_perms(`ncl')
        else mata: cn__clperm = (1..`ncl')
    }
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
            strtoreal(tokens("`taus'")), "`csufs'", "`cormethod'",                ///
            strtoreal(tokens("`alphas'")), `nfolds', "`loss'", "`folds'",         ///
            ("`standardize'" != ""), "`dccparams'", ("`adaptive'" != ""))
        capture frame drop `f'
        frame create `f'
        frame `f' {
            local covlab = cond("`model'" == "dccgarch", "Conditional", "Residual")
            mata: cn_postframe("`timevar'", "`tfmt'", `savef', "`datevar'", "`dfmt'", ///
                `savecov', "`names'", "`covlab'")
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
            char _dta[cn_model]     "`model'"
            char _dta[cn_cov]       "`savecov'"
            if `ncl' {
                tempname clsh`i' clct`i'
                mata: cn_postclust("`clnames'", "`fev`i''", "`clsh`i''", "`clct`i''")
                char _dta[cn_clusters] "`clnames'"
            }
        }
    }
    if `ncl' capture mata: mata drop cn__clid cn__clperm
    local sc `sc1'

    // spillover asymmetry measure: positive minus negative (Barunik, Kocenda and Vacha
    // 2016), in the main results frame
    if "`asymmetric'" != "" {
        local samv tci
        local samn sam
        forvalues j = 1/`k' {
            local samv `samv' to_`j' from_`j' net_`j'
            local samn `samn' sam_to_`j' sam_from_`j' sam_net_`j'
        }
        frame `frame'_pos: mata: cn__sp = st_data(., tokens("`samv'"))
        frame `frame'_neg: mata: cn__sn = st_data(., tokens("`samv'"))
        frame `frame' {
            local j 0
            foreach v of local samn {
                local ++j
                quietly generate double `v' = .
            }
            mata: st_store(., tokens("`samn'"), cn__sp - cn__sn)
            label variable sam "SAM: TCI of positive minus TCI of negative"
            forvalues j = 1/`k' {
                local nj : word `j' of `names'
                label variable sam_to_`j'   "SAM: TO, `nj'"
                label variable sam_from_`j' "SAM: FROM, `nj'"
                label variable sam_net_`j'  "SAM: NET, `nj'"
            }
        }
        mata: mata drop cn__sp cn__sn
    }

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
    if `enet' {
        ereturn scalar nfolds = `nfolds'
        ereturn local alphas "`alphas'"
        ereturn local loss   "`loss'"
        ereturn local folds  "`folds'"
        ereturn local seed   "`seed'"
        ereturn local adaptive "`adaptive'"
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
        ereturn scalar SAM = `sc2'[1,1] - `sc3'[1,1]
        if `semi' {
            ereturn local positive "`positive'"
            ereturn local negative "`negative'"
        }
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
    if `gen' & "`approach'" == "time" & !`dcc' ereturn local normalize "`normalize'"
    if `gen' > 1 ereturn local fevdtype = "generalized, " + cond(`gen' == 2, "max row sum", ///
        cond(`gen' == 3, "max column sum", "spectral radius"))
    if "`approach'" == "r2" ereturn local fevdtype = cond(`dccr2', ///
        "R2 decomposition (DCC correlations)", "R2 decomposition (`cormethod')")
    if "`approach'" == "modelfree" ereturn local fevdtype "squared correlations (`cormethod')"
    if inlist("`approach'", "r2", "modelfree") ereturn scalar horizon = .
    if `dcc' {
        tempname gp dp
        matrix `gp' = `tab1'_g
        matrix `dp' = `tab1'_d
        matrix rownames `gp' = `names'
        matrix colnames `gp' = mu omega alpha beta
        ereturn matrix garch = `gp'
        matrix `dp' = `dp'[1, 1..2]
        matrix colnames `dp' = a b
        ereturn scalar ll        = `tab1'_d[1, 3]
        ereturn scalar converged = `tab1'_d[1, 4]
        ereturn matrix dcc = `dp'
        ereturn scalar lags = .
        ereturn local standardize "`standardize'"
        ereturn local dccparams   "`dccparams'"
        ereturn local fevdtype = cond("`standardize'" != "", "volatility FEVD (standardized)", ///
            "volatility FEVD")
        if `dccr2' ereturn local fevdtype "R2 decomposition (DCC correlations)"
        matrix drop `tab1'_g `tab1'_d
    }
    if `ncl' {
        local rn `clnames' average
        matrix rownames `clsh1' = `rn'
        matrix colnames `clsh1' = own comove contag to net
        matrix rownames `clct1' = `clnames'
        matrix colnames `clct1' = `clnames'
        matrix rownames `clct1'_avg = `clnames'
        matrix colnames `clct1'_avg = `clnames'
        matrix rownames `clct1'_tab = `clnames' TO NET
        matrix colnames `clct1'_tab = `clnames' FROM
        ereturn scalar cl_within = `clct1'_sys[1, 1]
        ereturn scalar cl_cross  = `clct1'_sys[1, 2]
        matrix drop `clct1'_sys
        ereturn matrix cl_shares = `clsh1'
        ereturn matrix cl_contagion = `clct1'
        ereturn matrix cl_contagion_avg = `clct1'_avg
        ereturn matrix cl_table = `clct1'_tab
        ereturn scalar nclusters  = `ncl'
        ereturn scalar norderings = `nord'
        ereturn local clusters "`clnames'"
        ereturn local clorder  "`clorder'"
        tempname cm
        matrix `cm' = `CLID''
        matrix colnames `cm' = `names'
        ereturn matrix cl_id = `cm'
        ereturn local fevdtype "clustered (`ncl' clusters)"
        if `semi' | "`asymmetric'" != "" {
            capture matrix drop `clsh2' `clct2' `clct2'_avg `clsh3' `clct3' `clct3'_avg
            capture matrix drop `clct2'_tab `clct2'_sys `clct3'_tab `clct3'_sys
        }
    }
    ereturn local approach "`approach'"
    ereturn local model    "`model'"
    ereturn scalar covariance = `savecov'
    ereturn local title    "Connectedness approach"
    ereturn local estat_cmd "connectedness_estat"
    ereturn local cmd      "connectedness"

    Display, digits(`digits') `table'
    if `dcc' & e(converged) == 0 {
        di as txt "{p 0 6 2}note: the maximum likelihood optimization did not report convergence; " ///
            "check the estimates in {cmd:e(garch)} and {cmd:e(dcc)}{p_end}"
    }
    if "`r2fallback'" != "" {
        if `r2fallback' > 0 {
            di as txt "{p 0 6 2}note: in the first period, the R2 shares of `r2fallback' " ///
                "variable(s) sum to more than 1; they are kept (R's R2Correlations() stops here)," ///
                " while later periods take the previous period's shares, as in R{p_end}"
        }
    }
    if "`covnote'" != "" {
        di as txt "{p 0 6 2}note: the covariances were not stored in the results frame, " ///
            "because they would exceed {help maxvar}; {cmd:cnportfolio} can still use the PCI{p_end}"
    }
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
    if "`e(model)'" == "lasso" local mname "LASSO-VAR"
    if "`e(model)'" == "ridge" local mname "ridge VAR"
    if "`e(model)'" == "elasticnet" local mname "elastic-net VAR"
    if "`e(adaptive)'" != "" local mname "adaptive `mname'"
    local wdesc = cond(e(window) > 0, "rolling window", "static")
    if "`approach'" == "r2" & "`e(model)'" == "dccgarch" local mdesc "DCC-GARCH(1,1), full sample"
    else if "`approach'" == "r2" local mdesc "R2, `e(lags)' lag(s), `wdesc'"
    else if "`approach'" == "modelfree" local mdesc "model-free, `wdesc'"
    else if "`e(model)'" == "tvpvar" local mdesc "TVP-VAR(`e(lags)'), `e(prior)' prior"
    else if "`e(model)'" == "dccgarch" local mdesc "DCC-GARCH(1,1), full sample"
    else if e(window) > 0      local mdesc "`mname'(`e(lags)'), rolling window"
    else                       local mdesc "`mname'(`e(lags)'), static"
    if length("`mdesc'") > 39 local mdesc : subinstr local mdesc ", rolling window" ", rolling"

    di _n as txt "`e(title)'" _col(49) "Number of variables = " as res %8.0g e(k)
    if "`approach'" == "r2" & "`e(model)'" == "dccgarch" {
        di as txt "Model:  " as res "`mdesc'" _col(49) as txt "Lags                = " as res %8.0g 0
    }
    else if inlist("`approach'", "r2", "modelfree") {
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
    else if "`e(model)'" == "dccgarch" {
        di as txt "DCC(1,1): a = " as res %6.4f el(e(dcc), 1, 1) as txt ", b = " ///
            as res %6.4f el(e(dcc), 1, 2) _col(49) as txt "Log likelihood = " as res %11.3f e(ll)
    }
    else if inlist("`e(model)'", "lasso", "ridge", "elasticnet") {
        local na : word count `e(alphas)'
        local atxt = cond(`na' == 1, "`e(alphas)'", "`na' values")
        di _col(49) as txt "CV: " as res "`e(nfolds)'" as txt " `e(folds)' folds, alpha " as res "`atxt'"
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
    else if "`approach'" == "r2" & "`e(cormethod)'" == "dcc" {
        di as txt "Approach: " as res "R2 decomposed connectedness" as txt " (" ///
            as res "DCC-GARCH conditional correlations" as txt ")"
    }
    else if "`approach'" == "r2" {
        di as txt "Approach: " as res "R2 decomposed connectedness" as txt " (" as res "`e(cormethod)'" ///
            as txt " correlations)"
    }
    else if "`e(model)'" == "dccgarch" {
        di as txt "Approach: " as res "volatility connectedness" as txt " (Gabauer 2020)"
    }
    else if "`approach'" == "modelfree" {
        di as txt "Approach: " as res "model-free connectedness" as txt " (" as res "`e(cormethod)'" ///
            as txt " correlations)"
    }
    local fr "`e(frame)'"
    if inlist("`suf'", "_pos", "_neg") local fr "`e(frame`suf')'"
    di as txt "Results frame: " as res "`fr'" as txt " (first period " as res "`e(firstm)'" as txt ")"
    if "`e(clusters)'" != "" {
        local otxt = cond("`e(clorder)'" == "average", ///
            "averaged over " + string(e(norderings)) + " ordering" + cond(e(norderings) > 1, "s", ""), ///
            "ordered as listed")
        di as txt "Clusters: " as res "`e(clusters)'" as txt " (`otxt')"
    }
    if "`e(asymmetric)'" != "" {
        local ptxt = cond("`suf'" == "_pos", "positive", cond("`suf'" == "_neg", "negative", "all"))
        local stxt = cond("`e(positive)'" != "", "series", "changes")
        di as txt "Asymmetric connectedness: " as res "`ptxt'" as txt " `stxt' shown; " ///
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
            di as txt "Rows: R2 share of the row variable's variance explained by the column" _n ///
                "variable (percent)."
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
        local stxt = cond("`e(positive)'" != "", "series", "changes")
        di as txt "Positive `stxt': TCI = " as res %6.`digits'f e(TCI_pos) ///
            as txt "    Negative `stxt': TCI = " as res %6.`digits'f e(TCI_neg)
        di as txt "Spillover asymmetry measure: SAM = " as res %6.`digits'f e(SAM)
    }
    if "`approach'" == "r2" & "`e(TCI_c)'" != "" & "`suf'" == "" {
        di _n as txt %-16s "Component" %12s "TCI" %16s "corrected TCI"
        di as txt "{hline 44}"
        di as txt %-16s "Contemporaneous" as res %12.`digits'f e(TCI_c) %16.`digits'f e(cTCI_c)
        di as txt %-16s "Lagged" as res %12.`digits'f e(TCI_l) %16.`digits'f e(cTCI_l)
        di as txt "{hline 44}"
        di as txt %-16s "Overall" as res %12.`digits'f e(TCI) %16.`digits'f e(cTCI)
    }
    if "`e(clusters)'" != "" & "`suf'" == "" & "`table'" == "" {
        di _n as txt "Clustered connectedness (percent): own, co-movement (same cluster) and"
        di as txt "contagion (other clusters) shares, averaged over the variables of each cluster"
        matlist e(cl_shares), format(%`w'.`digits'f) border(rows) lines(oneline) nodotz ///
            rowtitle("Cluster") twidth(12)
        di as txt "TO and NET: contagion sent to other clusters, per variable of the cluster"
        di _n as txt "{p 0 0 2}Cluster table: share of the row cluster's variance due to the column " ///
            "cluster; FROM: row total, TO: column total, outside the own cluster{p_end}"
        matlist e(cl_table), format(%`w'.`digits'f) border(rows) lines(oneline) nodotz ///
            rowtitle("Cluster") twidth(12)
        di as txt "System-wide: TCI " as res %6.`digits'f (e(cl_within) + e(cl_cross)) as txt " = within " ///
            as res %6.`digits'f e(cl_within) as txt " + cross " as res %6.`digits'f e(cl_cross)
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
