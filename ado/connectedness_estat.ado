*! version 0.9.1  24oct2026
*! postestimation statistics after -connectedness-
program define connectedness_estat
    version 16.1
    if "`e(cmd)'" != "connectedness" error 301
    gettoken sub rest : 0, parse(" ,")
    local lsub = length("`sub'")
    if `lsub' >= 3 & "`sub'" == substr("table", 1, `lsub') {
        connectedness `rest'
    }
    else if `lsub' >= 3 & "`sub'" == substr("group", 1, `lsub') {
        Group `rest'
    }
    else if `lsub' >= 4 & "`sub'" == substr("correlation", 1, `lsub') {
        Correlation `rest'
    }
    else if `lsub' >= 4 & "`sub'" == substr("bootstrap", 1, `lsub') {
        Bootstrap `rest'
    }
    else {
        estat_default `0'
    }
end

// Group measures: aggregated, conditional, inclusive, exclusive, internal, external
program define Group, rclass
    syntax , Type(string) GRoups(string asis) [FRame(name) REPLACE DIGits(integer 2) noTABle]

    local type = lower("`type'")
    local types aggregated conditional inclusive exclusive internal external
    local ltype = length("`type'")
    local found
    foreach t of local types {
        if `ltype' >= 3 & "`type'" == substr("`t'", 1, `ltype') local found `t'
    }
    if "`found'" == "" {
        di as err "type() must be one of: `types'"
        exit 198
    }
    local type `found'
    local single = inlist("`type'", "conditional", "inclusive", "exclusive")
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }

    if !inlist("`e(approach)'", "", "time", "modelfree") {
        di as err "estat group requires results estimated with approach(time) or approach(modelfree)"
        exit 198
    }
    if !inlist("`e(normalize)'", "", "row") {
        di as err "estat group requires the row-normalized decomposition; " ///
            "the results use normalize(`e(normalize)')"
        exit 198
    }

    // results frame and its FEVD
    local fr `e(frame)'
    mata: st_local("fexist", strofreal(st_frameexists("`fr'")))
    if !`fexist' {
        di as err "results frame {bf:`fr'} not found; rerun {cmd:connectedness}"
        exit 111
    }
    frame `fr': local hasfevd : char _dta[cn_fevd]
    if "`hasfevd'" != "1" {
        di as err "estat group needs the decomposition for every period;"
        di as err "rerun {cmd:connectedness} without the {cmd:nofevd} option"
        exit 198
    }
    local names `e(names)'
    local k = e(k)

    // parse groups: "name: var var ..." "name: var ..." ; one group for single types
    if `single' {
        local allv
        foreach w of local groups {
            local allv `allv' `w'
        }
        local groups `""`allv'""'
    }
    local m : word count `groups'
    tempname G
    matrix `G' = J(`k', 1, 0)
    local gnames
    forvalues g = 1/`m' {
        local spec : word `g' of `groups'
        local gname "G`g'"
        local pos = strpos(`"`spec'"', ":")
        if `pos' {
            local gname = strtrim(substr(`"`spec'"', 1, `pos' - 1))
            local spec  = substr(`"`spec'"', `pos' + 1, .)
            if "`gname'" == "" local gname "G`g'"
        }
        local gname = strtoname("`gname'")
        local gnames `gnames' `gname'
        local nv 0
        foreach v of local spec {
            local p : list posof "`v'" in names
            if `p' == 0 {
                di as err "`v' is not one of the variables: `names'"
                exit 198
            }
            if `G'[`p', 1] != 0 {
                di as err "`v' appears in more than one group"
                exit 198
            }
            matrix `G'[`p', 1] = `g'
            local ++nv
        }
        if `nv' == 0 {
            di as err "group `g' is empty"
            exit 198
        }
    }
    if `single' {
        local m 1
        if "`type'" == "conditional" & `nv' < 2 {
            di as err "type(conditional) needs at least two variables"
            exit 198
        }
    }
    else if `m' < 2 {
        di as err "type(`type') needs at least two groups, each in double quotes"
        exit 198
    }

    // names of the units in the output
    if "`type'" == "aggregated" {
        local unames `gnames'
    }
    else if "`type'" == "conditional" {
        local unames
        forvalues i = 1/`k' {
            if `G'[`i', 1] == 1 local unames `unames' `: word `i' of `names''
        }
    }
    else local unames `names'

    // output frame
    if "`frame'" == "" local frame cn_group
    if "`frame'" == "`c(frame)'" | "`frame'" == "`fr'" {
        di as err "frame(`frame') may not be the current or the results frame"
        exit 498
    }
    mata: st_local("gexist", strofreal(st_frameexists("`frame'")))
    if `gexist' & "`replace'" == "" {
        di as err "frame {bf:`frame'} already exists; specify {cmd:replace} or another {cmd:frame()}"
        exit 110
    }

    capture mata: mata describe cn_group()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }
    local corr = ("`e(corrected)'" != "")
    local tvar `e(timevar)'
    local dvar `e(datevar)'
    tempname tab sc gt
    frame `fr' {
        local tfmt : format `tvar'
        local dfmt
        if "`dvar'" != "" local dfmt : format `dvar'
        mata: cn__gtv = st_data(., "`tvar'")
        if "`dvar'" != "" mata: cn__gdv = st_data(., "`dvar'")
        else mata: cn__gdv = J(0, 1, .)
        mata: cn_group("`type'", st_matrix("`G'"), `m', `corr', "`tab'", "`sc'", "`gt'")
    }
    local n = `sc'[1,3]
    matrix rownames `tab' = `unames' TO Inc_Own NET NPT
    matrix colnames `tab' = `unames' FROM

    capture frame drop `frame'
    frame create `frame'
    frame `frame' {
        mata: cn_postgroup("`type'", "`unames'", cn__gtv, "`tvar'", "`tfmt'", ///
            cn__gdv, "`dvar'", "`dfmt'", "`gnames'")
        quietly tsset `tvar'
        label data "connectedness group measures (`type')"
        char _dta[cn_kind]      "group"
        char _dta[cn_gtype]     "`type'"
        char _dta[cn_names]     "`unames'"
        char _dta[cn_gnames]    "`gnames'"
        char _dta[cn_k]         "`n'"
        char _dta[cn_timevar]   "`tvar'"
        char _dta[cn_datevar]   "`dvar'"
        char _dta[cn_corrected] "`e(corrected)'"
    }
    mata: mata drop cn__gtv cn__gdv

    // display (notable suppresses only the table)
    {
        local w = `digits' + 6
        local ttl = upper(substr("`type'", 1, 1)) + substr("`type'", 2, .)
        di _n as txt "`ttl' connectedness" _col(49) "Units in table      = " as res %8.0g `n'
        forvalues g = 1/`m' {
            local gl
            forvalues i = 1/`k' {
                if `G'[`i', 1] == `g' local gl `gl' `: word `i' of `names''
            }
            local gn : word `g' of `gnames'
            if `single' di as txt "Group: " as res "`gl'"
            else di as txt "Group " as res "`gn'" as txt ": " as res "`gl'"
        }
        di as txt "Results frame: " as res "`frame'"
        di
        if "`table'" == "" {
            matlist `tab', format(%`w'.`digits'f) border(rows) lines(oneline) nodotz ///
                rowtitle(" ") twidth(12)
        }
        di as txt "TCI = " as res %6.`digits'f `sc'[1,1] as txt "    corrected TCI = " ///
            as res %6.`digits'f `sc'[1,2]
        if "`type'" == "aggregated" {
            di as txt "TCI_ext (ConnectednessApproach) = " as res %6.`digits'f `gt'[1,1] ///
                as txt "    corrected = " as res %6.`digits'f `gt'[1,2]
        }
        if inlist("`type'", "internal", "external") {
            forvalues g = 1/`m' {
                local gn : word `g' of `gnames'
                di as txt "Group TCI, `gn' = " as res %6.`digits'f `gt'[1,`g'] ///
                    as txt "    corrected = " as res %6.`digits'f `gt'[1,`m'+`g']
            }
        }
    }

    return scalar TCI  = `sc'[1,1]
    return scalar cTCI = `sc'[1,2]
    return scalar n    = `n'
    return local type  "`type'"
    return local names "`unames'"
    return local groups "`gnames'"
    return local frame "`frame'"
    if inlist("`type'", "aggregated", "internal", "external") return matrix gtci = `gt'
    return matrix table = `tab'
end

// Conditional or partial correlations of the covariance matrices stored in the
// results frame, as R's ConditionalCorrelation() and PartialCorrelations()
program define Correlation, rclass
    syntax [, PARtial FRame(name) REPLACE DIGits(integer 3) noTABle]
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }
    local fr `e(frame)'
    mata: st_local("fexist", strofreal(st_frameexists("`fr'")))
    if !`fexist' {
        di as err "results frame {bf:`fr'} not found; rerun {cmd:connectedness}"
        exit 111
    }
    frame `fr': local hascov : char _dta[cn_cov]
    if "`hascov'" != "1" {
        di as err "estat correlation needs the covariance matrices of the results frame;"
        if inlist("`e(approach)'", "r2", "modelfree") & "`e(model)'" != "dccgarch" ///
            di as err "approach(`e(approach)') does not estimate a covariance matrix"
        else di as err "rerun {cmd:connectedness} without the {cmd:nocovariance} option"
        exit 498
    }
    local names `e(names)'
    local k = e(k)
    local ispart = ("`partial'" != "")
    if "`frame'" == "" local frame cn_corr
    if "`frame'" == "`c(frame)'" | "`frame'" == "`fr'" {
        di as err "frame(`frame') may not be the current or the results frame"
        exit 498
    }
    mata: st_local("cexist", strofreal(st_frameexists("`frame'")))
    if `cexist' & "`replace'" == "" {
        di as err "frame {bf:`frame'} already exists; specify {cmd:replace} or another {cmd:frame()}"
        exit 110
    }
    capture mata: mata describe cn_corr()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }
    local tvar `e(timevar)'
    local dvar `e(datevar)'
    tempname tab
    frame `fr' {
        local tfmt : format `tvar'
        local dfmt
        if "`dvar'" != "" local dfmt : format `dvar'
        local T = _N
        mata: cn__ctv = st_data(., "`tvar'")
        if "`dvar'" != "" mata: cn__cdv = st_data(., "`dvar'")
        else mata: cn__cdv = J(0, 1, .)
        mata: cn_corr(`k', `ispart', "`tab'")
    }
    capture frame drop `frame'
    frame create `frame'
    frame `frame' {
        mata: cn_postcorr(cn__ctv, "`tvar'", "`tfmt'", cn__cdv, "`dvar'", "`dfmt'", "`names'", `ispart')
        quietly tsset `tvar'
        local dlab = cond(`ispart', "partial correlations", "conditional correlations")
        label data "connectedness: `dlab'"
        char _dta[cn_kind]    "correlation"
        char _dta[cn_names]   "`names'"
        char _dta[cn_timevar] "`tvar'"
        char _dta[cn_datevar] "`dvar'"
    }
    mata: mata drop cn__ctv cn__cdv cn__CR
    matrix rownames `tab' = `names'
    matrix colnames `tab' = `names'

    local what = cond(`ispart', "Partial correlations", "Conditional correlations")
    local src = cond("`e(model)'" == "dccgarch", "DCC-GARCH conditional covariances", ///
        cond(e(window) > 0 | "`e(model)'" == "tvpvar", "residual covariances of each period", ///
        "residual covariance"))
    di _n as txt "`what'" _col(49) "Number of variables = " as res %8.0g `k'
    di as txt "Source: " as res "`src'" as txt _col(49) "Periods averaged    = " as res %8.0g `T'
    di as txt "Saved in frame: " as res "`frame'"
    if "`table'" == "" {
        local w = `digits' + 6
        matlist `tab', format(%`w'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(12)
    }
    return scalar T = `T'
    return scalar k = `k'
    return local type = cond(`ispart', "partial", "conditional")
    return local frame "`frame'"
    return matrix mean = `tab'
end

// Bootstrap inference (Greenwood-Nimmo, Kocenda and Nguyen 2024): residual bootstrap
// of the VAR of each window with bias correction; bands, event probabilities and
// detection of significant rises
program define Bootstrap, rclass
    syntax [, Reps(integer 1000) noBIAS Level(cilevel) SEED(string) FRame(name) REPLACE ///
        EVERY(integer 1) ALL EVENt(string) LEADs(numlist integer >=0) REFerence(integer 1) ///
        THReshold(real 0) DECrease DETect(real -1) CUToff(real 90) noDOTs ///
        DIGits(integer 2) TESTINDices(name) ]

    if "`e(model)'" != "var" | !inlist("`e(approach)'", "", "time") {
        di as err "estat bootstrap requires results of {cmd:connectedness} with model(var) and approach(time)"
        exit 198
    }
    if `reps' < 2 {
        di as err "reps() must be at least 2"
        exit 198
    }
    if `every' < 1 {
        di as err "every() must be a positive integer"
        exit 198
    }
    if `reference' < 1 {
        di as err "reference() must be a positive integer"
        exit 198
    }
    if `threshold' < 0 {
        di as err "threshold() must not be negative"
        exit 198
    }
    if `cutoff' <= 0 | `cutoff' >= 100 {
        di as err "cutoff() must be between 0 and 100"
        exit 198
    }
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }
    local static = (e(window) == 0)
    local isev = (`"`event'"' != "")
    local isdet = (`detect' >= 0)
    if `isev' & `static' {
        di as err "event() requires rolling-window results"
        exit 198
    }
    if `isdet' & (`static' | `every' != 1 | `isev') {
        di as err "detect() requires rolling-window results, every(1) and no event()"
        exit 198
    }
    if !`isev' & ("`leads'" != "" | `reference' != 1 | `threshold' != 0 | "`decrease'" != "") {
        di as err "options leads(), reference(), threshold() and decrease require event()"
        exit 198
    }
    if "`leads'" == "" local leads 0 1 5 22

    // data of the estimation, as in connectedness
    quietly count if e(sample)
    if r(N) != e(N) {
        di as err "the estimation sample is not in memory; restore the data used by {cmd:connectedness}"
        exit 498
    }
    tsrevar `e(depvar)'
    local vars `r(varlist)'
    tempvar touse
    quietly generate byte `touse' = e(sample)
    local p   = e(lags)
    local H   = e(horizon)
    local w   = e(window)
    local T0  = e(T0)
    local k   = e(k)
    local names `e(names)'

    // decomposition: generalized (1), its scalar normalizations (2-4), clustered
    // (5), orthogonalized (0)
    local gen = ("`e(fevdtype)'" != "orthogonalized")
    if "`e(normalize)'" == "maxrow" local gen 2
    if "`e(normalize)'" == "maxcol" local gen 3
    if "`e(normalize)'" == "spectral" local gen 4
    if "`e(clusters)'" != "" {
        local gen 5
        tempname cid
        matrix `cid' = e(cl_id)
        mata: cn__clid = st_matrix("`cid'")'
        if "`e(clorder)'" == "average" mata: cn__clperm = cn_perms(`=e(nclusters)')
        else mata: cn__clperm = (1..`=e(nclusters)')
    }
    capture mata: mata describe cn_bootstrap()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }

    // windows of the results frame and their dates
    local fr `e(frame)'
    local tv `e(timevar)'
    local dv `e(datevar)'
    mata: st_local("fexist", strofreal(st_frameexists("`fr'")))
    if !`fexist' {
        di as err "results frame {bf:`fr'} not found; rerun {cmd:connectedness}"
        exit 111
    }
    frame `fr' {
        local tfmt : format `tv'
        mata: cn__btv = st_data(., "`tv'")
        if "`dv'" != "" {
            local dfmt : format `dv'
            mata: cn__bdv = st_data(., "`dv'")
        }
        else mata: cn__bdv = cn__btv
    }
    if "`dv'" == "" local dfmt `tfmt'

    // events: dates (read with the date variable, or the time variable if it holds
    // daily dates) or values of the time variable
    local evlab
    if `isev' {
        local nev 0
        mata: cn__bev = J(0, 1, .)
        foreach tok of local event {
            local ++nev
            local val = real("`tok'")
            local mv cn__btv
            if missing(`val') {
                local val = date("`tok'", "DMY")
                if missing(`val') {
                    di as err "event(): `tok' is neither a number nor a date such as 17oct1997"
                    exit 198
                }
                if substr("`dfmt'", 1, 3) != "%td" {
                    di as err "event(): dates require a daily date variable; give values of `tv' instead"
                    exit 198
                }
                local mv cn__bdv
            }
            mata: st_local("re", strofreal(min(selectindex(`mv' :>= `val'))))
            if "`re'" == "." | "`re'" == "" {
                di as err "event(): `tok' is after the last window"
                exit 198
            }
            if `re' <= `reference' {
                di as err "event(): `tok' leaves fewer than `reference' earlier window(s)"
                exit 198
            }
            mata: cn__bev = cn__bev \ `re'
            local evlab `evlab' `tok'
        }
        // windows needed: the reference windows and the leads of each event
        mata: cn__bw = J(0, 1, .)
        mata: cn__bl = strtoreal(tokens("`leads'"))
        mata: for (cn__i=1; cn__i<=rows(cn__bev); cn__i++) cn__bw = cn__bw \ ((cn__bev[cn__i] - `reference')::(cn__bev[cn__i] - 1)) \ (cn__bev[cn__i] :+ cn__bl')
        mata: cn__bw = uniqrows(select(cn__bw, cn__bw :<= `T0'))
    }
    else if `static' {
        mata: cn__bw = 1
    }
    else {
        mata: cn__bw = uniqrows(((1::`T0')[selectindex(mod((1::`T0') :- 1, `every') :== 0)]) \ `T0')
    }
    mata: st_local("nw", strofreal(rows(cn__bw)))

    if "`seed'" != "" set seed `seed'
    if "`testindices'" != "" {
        confirm matrix `testindices'
        mata: cn__bidx = st_matrix("`testindices'")
    }
    local dobias = ("`bias'" == "")
    local doall  = ("`all'" != "" | `static')
    di as txt _n "Bootstrap of `nw' window(s), `reps' draws" ///
        cond(`dobias', " per round, two rounds (bias correction)", ", no bias correction") ///
        cond("`dots'" == "", ":", "")
    capture noisily mata: cn_bootstrap("`vars'", "`touse'", `p', `H', `w', `gen', `reps', ///
        `dobias', `doall', cn__bw, `level', ("`dots'" == ""))
    local rc = _rc
    capture mata: mata drop cn__bidx
    if `gen' == 5 capture mata: mata drop cn__clid cn__clperm
    if `rc' exit `rc'
    mata: st_local("ndisc", strofreal(sum(cn__BN)))
    mata: st_local("ncarry", strofreal(sum(cn__BST :== 0)))
    mata: st_local("nfail", strofreal(missing(cn__BST)))

    // ------------------------------------------------------------ event test
    if `isev' {
        local nl : word count `leads'
        tempname P
        mata: cn__bpos = J(`T0', 1, 0); cn__bpos[cn__bw] = (1::rows(cn__bw))
        mata: st_matrix("`P'", cn_bootprob(cn__bev, cn__bl, `reference', `threshold', ("`decrease'" != ""), cn__bpos))
        local cn
        foreach l of local leads {
            local cn `cn' lead`l'
        }
        matrix colnames `P' = `cn'
        local rn
        forvalues e = 1/`nev' {
            local rn `rn' event`e'
        }
        matrix rownames `P' = `rn'
        local cmp = cond("`decrease'" != "", "below", "above")
        local reft = cond(`reference' == 1, "the window before the event", ///
            "the mean of the `reference' windows before the event")
        local thrt = cond(`threshold' > 0, " by at least `threshold' percent", "")
        di _n as txt "Probability (percent) that the bias-corrected TCI is `cmp' its level in"
        di as txt "`reft'`thrt', by lead (windows after the event window)"
        di as txt "{hline 26}{c TT}{hline `=10*`nl''}"
        di as txt %-12s "Event" %14s "Window" "{c |}" _c
        foreach l of local leads {
            di as txt %10s "`l'" _c
        }
        di _n as txt "{hline 26}{c +}{hline `=10*`nl''}"
        forvalues e = 1/`nev' {
            local tok : word `e' of `evlab'
            mata: st_local("re", strofreal(cn__bev[`e']))
            mata: st_local("wd", strofreal(cn__bdv[`re']))
            local wds : display `dfmt' `wd'
            di as txt %-12s "`tok'" %14s "`=strtrim("`wds'")'" "{c |}" _c
            forvalues j = 1/`nl' {
                local v = `P'[`e', `j']
                local star = cond(`v' >= `cutoff' & `v' < ., "*", " ")
                if `v' < . di as res %9.1f `v' "`star'" _c
                else di as txt %10s "." _c
            }
            di
        }
        di as txt "{hline 26}{c BT}{hline `=10*`nl''}"
        di as txt "* at least `cutoff' percent. Windows are dated by their last observation."
        return matrix prob = `P'
        return local events "`evlab'"
        return local leads "`leads'"
        return scalar reference = `reference'
        return scalar threshold = `threshold'
    }
    // ------------------------------------------------------------ static table
    else if `static' {
        tempname Tb
        mata: st_matrix("`Tb'", (cn__BP', cn__BB', cn__BC', cn__BS', cn__BL', cn__BU'))
        local rn TCI
        foreach s in TO FROM NET {
            forvalues i = 1/`k' {
                local nm : word `i' of `names'
                local nm = strtoname("`nm'")
                local rn `rn' `s':`nm'
            }
        }
        matrix rownames `Tb' = `rn'
        matrix colnames `Tb' = estimate bias corrected se lower upper
        local w2 = `digits' + 7
        di _n as txt "Bootstrap of the connectedness measures (percent): estimate, bias,"
        di as txt "bias-corrected estimate, standard error and `level'% interval"
        matlist `Tb', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(12)
        return matrix table = `Tb'
    }
    // ------------------------------------------------------------ bands frame
    else {
        if "`frame'" == "" local frame cn_boot
        if "`frame'" == "`c(frame)'" | "`frame'" == "`fr'" {
            di as err "frame(`frame') may not be the current or the results frame"
            exit 498
        }
        mata: st_local("bexist", strofreal(st_frameexists("`frame'")))
        if `bexist' & "`replace'" == "" {
            di as err "frame {bf:`frame'} already exists; specify {cmd:replace} or another {cmd:frame()}"
            exit 110
        }
        local mn tci
        if `doall' {
            foreach s in to from net {
                forvalues i = 1/`k' {
                    local mn `mn' `s'_`i'
                }
            }
        }
        if `isdet' {
            tempname PR
            mata: cn__bpos = (1::`T0')
            mata: cn__P = J(`T0', 1, .); for (cn__i=2; cn__i<=`T0'; cn__i++) cn__P[cn__i] = cn_bootprob(cn__i, 0, 1, `detect', 0, cn__bpos)
            mata: st_matrix("`PR'", cn__P)
        }
        capture frame drop `frame'
        frame create `frame'
        frame `frame' {
            quietly set obs `nw'
            mata: st_store(., st_addvar("double", "`tv'"), cn__btv[cn__bw])
            format `tv' `tfmt'
            if "`dv'" != "" {
                mata: st_store(., st_addvar("double", "`dv'"), cn__bdv[cn__bw])
                format `dv' `dfmt'
            }
            local j 0
            foreach m of local mn {
                local ++j
                mata: st_store(., st_addvar("double", ("`m'", "`m'_bias", "`m'_bc", "`m'_se", "`m'_lo", "`m'_hi")), (cn__BP[., `j'], cn__BB[., `j'], cn__BC[., `j'], cn__BS[., `j'], cn__BL[., `j'], cn__BU[., `j']))
            }
            mata: st_store(., st_addvar("double", "ndiscard"), cn__BN)
            mata: st_store(., st_addvar("byte", "stable"), cn__BST)
            label variable stable "VAR of the window stable (0: estimates of the last stable window)"
            label variable tci "TCI, estimate"
            label variable tci_bias "TCI, bootstrap bias"
            label variable tci_bc "TCI, bias-corrected"
            label variable tci_se "TCI, bootstrap standard error"
            label variable tci_lo "TCI, lower `level'% limit"
            label variable tci_hi "TCI, upper `level'% limit"
            label variable ndiscard "unstable draws discarded"
            if `isdet' {
                quietly generate double prob = .
                quietly generate byte flag = .
                mata: st_store(., "prob", cn__P)
                quietly replace flag = (prob >= `cutoff') if prob < .
                label variable prob "Probability (%) of a rise of at least `detect'% from the previous window"
                label variable flag "prob at least `cutoff'"
            }
            quietly tsset `tv'
            label data "connectedness: bootstrap (`reps' draws)"
            char _dta[cn_kind] "bootstrap"
            if `isdet' {
                quietly count if flag == 1
                local nflag = r(N)
            }
        }
        quietly frame `frame': summarize tci_bias, meanonly
        local mb = r(mean)
        di _n as txt "Bias-corrected TCI and `level'% bands saved in frame " as res "`frame'" ///
            as txt " (`nw' windows)"
        di as txt "Mean bootstrap bias of the TCI: " as res %6.`digits'f `mb'
        if `isdet' {
            di as txt "Windows with a probability of at least `cutoff'% of a rise of at least " ///
                "`detect'%: " as res `nflag'
            if `nflag' > 0 & `nflag' <= 30 {
                frame `frame': list `=cond("`dv'" != "", "`dv'", "`tv'")' prob if flag == 1, noobs clean
            }
            return scalar nflag = `nflag'
        }
        return local frame "`frame'"
    }
    di as txt "Unstable draws discarded: " as res `ndisc'
    if `ncarry' > 0 {
        di as txt "{p 0 6 2}note: `ncarry' window(s) with an unstable VAR use the estimates of the " ///
            "most recent window with a stable VAR, as Greenwood-Nimmo, Kocenda and Nguyen do{p_end}"
    }
    if `nfail' > 0 {
        di as txt "{p 0 6 2}note: `nfail' window(s) have no earlier window with a stable VAR; " ///
            "their results are missing{p_end}"
    }
    return scalar ncarried = `ncarry'
    return scalar reps = `reps'
    return scalar level = `level'
    return scalar nwindows = `nw'
    return scalar ndiscard = `ndisc'
    return local bias = cond(`dobias', "corrected", "none")
    capture mata: mata drop cn__BP cn__BB cn__BC cn__BS cn__BL cn__BU cn__BD cn__BN cn__BST cn__btv cn__bdv cn__bw
    capture mata: mata drop cn__bev cn__bl cn__bpos cn__P cn__i
end
