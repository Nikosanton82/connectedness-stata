*! version 0.9.5  25oct2026
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
    else if `lsub' >= 4 & "`sub'" == substr("clorder", 1, `lsub') {
        Clorder `rest'
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
        DIGits(integer 2) TESTINDices(name) PAIRwise NPDC STat(string) PWFRame(name) ///
        noTABle SAM ]

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
    local dopw = ("`pairwise'" != "")
    local donp = ("`npdc'" != "")
    if (`dopw' | `donp') & `isev' {
        di as err "options pairwise and npdc may not be combined with event()"
        exit 198
    }
    if "`pwframe'" != "" & ((!`dopw' & !`donp') | `static') {
        di as err "pwframe() requires pairwise or npdc with rolling-window results"
        exit 198
    }

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

    // measure tested by event() and detect(): stat(tci), stat(to|from|net name),
    // stat(pairwise|npdc name1 name2), from name1 to name2, or stat(sam) and
    // stat(sam to|from|net name). Indices into cn_bootext(): TCI, TO, FROM, NET,
    // pairwise table and NPDC matrix by rows (cn_bootfull(), lsam elements), then SAM,
    // its TO, FROM and NET, and the TCI of the positive and of the negative system
    local k2 = `k'*`k'
    local lsam = 1 + 3*`k' + 2*`k2'
    local sidx 1
    local slab TCI
    local stype tci
    local snames
    local ssam 0
    if `"`stat'"' != "" {
        if !`isev' & !`isdet' {
            di as err "stat() requires event() or detect()"
            exit 198
        }
        gettoken stype snames : stat
        local snames : list retokenize snames
        local stype = lower("`stype'")
        if "`stype'" == "pw" local stype pairwise
        // stat(sam [to|from|net name]): the measure of the spillover asymmetry
        if "`stype'" == "sam" {
            local ssam 1
            gettoken stype snames : snames
            local snames : list retokenize snames
            local stype = lower("`stype'")
            if "`stype'" == "" local stype tci
            else if !inlist("`stype'", "to", "from", "net") {
                di as err "stat(sam `stype'): only stat(sam) and stat(sam to|from|net name) are allowed"
                exit 198
            }
        }
        local nn : word count `snames'
        local need = cond("`stype'" == "tci", 0, cond(inlist("`stype'", "to", "from", "net"), 1, ///
            cond(inlist("`stype'", "pairwise", "npdc"), 2, -1)))
        if `need' < 0 {
            di as err "stat(): `stype' is not tci, to, from, net, pairwise or npdc"
            exit 198
        }
        if `nn' != `need' {
            di as err "stat(`stype') requires `need' variable name(s)"
            exit 198
        }
        local pos
        foreach nm of local snames {
            local q : list posof "`nm'" in names
            if `q' == 0 {
                di as err "stat(): `nm' is not one of the variables (`names')"
                exit 198
            }
            local pos `pos' `q'
        }
        local a : word 1 of `pos'
        local b : word 2 of `pos'
        local na : word 1 of `snames'
        local nb : word 2 of `snames'
        if "`stype'" == "to" {
            local sidx = 1 + `a'
            local slab "TO of `na'"
        }
        if "`stype'" == "from" {
            local sidx = 1 + `k' + `a'
            local slab "FROM of `na'"
        }
        if "`stype'" == "net" {
            local sidx = 1 + 2*`k' + `a'
            local slab "NET of `na'"
        }
        if inlist("`stype'", "pairwise", "npdc") {
            if `a' == `b' {
                di as err "stat(`stype') requires two different variables"
                exit 198
            }
            local sidx = 1 + 3*`k' + (`b' - 1)*`k' + `a'
            local slab "spillover from `na' to `nb'"
            if "`stype'" == "npdc" {
                local sidx = `sidx' + `k'*`k'
                local slab "net spillover from `na' to `nb'"
            }
        }
        if `ssam' {
            local sidx = `lsam' + `sidx'
            local slab = cond("`stype'" == "tci", "SAM", "SAM `slab'")
            local stype sam `stype'
            if "`stype'" == "sam tci" local stype sam
        }
        if (inlist("`stype'", "net", "npdc") | `ssam') & (`threshold' > 0 | `detect' > 0) {
            di as err "{p 0 0 2}threshold() and detect(#) with # > 0 are rises in percent, " ///
                "not defined for stat(`stype'), which may be zero or negative; " ///
                "use detect(0) for the probability of any rise{p_end}"
            exit 198
        }
    }
    // spillover asymmetry measure: the total, positive and negative systems are
    // bootstrapped jointly, with the same resampled rows
    local dosam = ("`sam'" != "" | `ssam')
    if `dosam' & "`e(asymmetric)'" == "" {
        di as err "{p 0 0 2}option sam and stat(sam) require results of {cmd:connectedness} " ///
            "with asymmetric or positive() and negative(){p_end}"
        exit 198
    }
    local nsys 1
    if `dosam' {
        if "`e(positive)'" != "" {
            tsrevar `e(positive)'
            local posvars `r(varlist)'
            tsrevar `e(negative)'
            local negvars `r(varlist)'
        }
        else {
            // positive and negative parts of each series, as in connectedness
            local posvars
            local negvars
            foreach v of local vars {
                tempvar pv nv
                quietly generate double `pv' = cond(missing(`v'), ., max(`v', 0)) if `touse'
                quietly generate double `nv' = cond(missing(`v'), ., min(`v', 0)) if `touse'
                local posvars `posvars' `pv'
                local negvars `negvars' `nv'
            }
        }
        local vars `vars' `posvars' `negvars'
        local nsys 3
    }

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
        local ncl = e(nclusters)
        local nctl = (e(cl_ncontrols) == 1)
        if "`e(clorder)'" == "average" {
            if `nctl' mata: cn__clperm = J(round(exp(lnfactorial(`ncl' - 1))), 1, 1), (cn_perms(`ncl' - 1) :+ 1)
            else mata: cn__clperm = cn_perms(`ncl')
        }
        else mata: cn__clperm = (1..`ncl')
        mata: cn__clsco = J(1, 1, ("`e(clorder)'" == "simultaneous"))
        mata: cn__clpairs = J(0, 0, .)
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
    // measures bootstrapped: TCI, optionally TO, FROM and NET (offset 2), the pairwise
    // table (offset opw) and the NPDC matrix (offset onp), and the tested measure
    local k2 = `k'*`k'
    mata: cn__bsel = 1
    if `doall' mata: cn__bsel = cn__bsel, (2..`=1+3*`k'')
    local opw 0
    local onp 0
    if `dopw' {
        mata: st_local("opw", strofreal(cols(cn__bsel) + 1))
        mata: cn__bsel = cn__bsel, (`=2+3*`k''..`=1+3*`k'+`k2'')
    }
    if `donp' {
        mata: st_local("onp", strofreal(cols(cn__bsel) + 1))
        mata: cn__bsel = cn__bsel, (`=2+3*`k'+`k2''..`=1+3*`k'+2*`k2'')
    }
    // SAM (offset osam), optionally its TO, FROM and NET, then the TCI of the positive
    // and of the negative system
    local osam 0
    if `dosam' {
        mata: st_local("osam", strofreal(cols(cn__bsel) + 1))
        mata: cn__bsel = cn__bsel, `=`lsam'+1'
        if `doall' mata: cn__bsel = cn__bsel, (`=`lsam'+2'..`=`lsam'+1+3*`k'')
        mata: cn__bsel = cn__bsel, `=`lsam'+2+3*`k'', `=`lsam'+3+3*`k''
    }
    mata: st_local("has", strofreal(anyof(cn__bsel, `sidx')))
    if !`has' mata: cn__bsel = cn__bsel, `sidx'
    mata: st_local("scol", strofreal(selectindex(cn__bsel :== `sidx')[1]))
    di as txt _n "Bootstrap of `nw' window(s), `reps' draws" ///
        cond(`dobias', " per round, two rounds (bias correction)", ", no bias correction") ///
        cond(`dosam', ", total, positive and negative systems jointly", "") ///
        cond("`dots'" == "", ":", "")
    capture noisily mata: cn_bootstrap("`vars'", "`touse'", `p', `H', `w', `gen', `reps', ///
        `dobias', cn__bsel, `scol', cn__bw, `level', ("`dots'" == ""), `nsys')
    local rc = _rc
    capture mata: mata drop cn__bidx
    if `gen' == 5 capture mata: mata drop cn__clid cn__clperm cn__clsco cn__clpairs cn__clnet
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
        di _n as txt "{p 0 0 2}Probability (percent) that the bias-corrected `slab' is `cmp' " ///
            "its level in `reft'`thrt', by lead (windows after the event window){p_end}"
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
        local m3 = 1 + 3*`k'
        mata: st_matrix("`Tb'", (cn__BP[|1, 1 \ 1, `m3'|]', cn__BB[|1, 1 \ 1, `m3'|]', ///
            cn__BC[|1, 1 \ 1, `m3'|]', cn__BS[|1, 1 \ 1, `m3'|]', cn__BL[|1, 1 \ 1, `m3'|]', ///
            cn__BU[|1, 1 \ 1, `m3'|]'))
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
        // pairwise table and NPDC matrix: k x k, row i and column j for the
        // (net) spillover from j to i
        local rnm
        foreach nm of local names {
            local rnm `rnm' `=strtoname("`nm'")'
        }
        foreach blk in pw npdc {
            local off = cond("`blk'" == "pw", `opw', `onp')
            if `off' == 0 continue
            foreach s in b bias bc se lo hi {
                local M = cond("`s'" == "b", "P", cond("`s'" == "bias", "B", cond("`s'" == "bc", ///
                    "C", cond("`s'" == "se", "S", cond("`s'" == "lo", "L", "U")))))
                tempname `blk'_`s'
                mata: st_matrix("``blk'_`s''", rowshape(cn__B`M'[1, (`off'..`off'+`k2'-1)], `k'))
                matrix rownames ``blk'_`s'' = `rnm'
                matrix colnames ``blk'_`s'' = `rnm'
            }
        }
        if `opw' & "`table'" == "" {
            di _n as txt "{p 0 0 2}Pairwise connectedness (percent), bias-corrected: row i, column j " ///
                "is the spillover from j to i{p_end}"
            matlist `pw_bc', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(12)
        }
        if `onp' & "`table'" == "" {
            di _n as txt "{p 0 0 2}Net pairwise directional connectedness (percent), bias-corrected: " ///
                "row i, column j is the net spillover from j to i{p_end}"
            CnStarMat `npdc_bc' `npdc_lo' `npdc_hi' `digits'
            di as txt "* `level'% interval excludes zero"
        }
        foreach blk in pw npdc {
            if `=cond("`blk'" == "pw", `opw', `onp')' == 0 continue
            foreach s in hi lo se bc bias b {
                return matrix `blk'_`s' = ``blk'_`s''
            }
        }
        // spillover asymmetry measure: SAM, its TO, FROM and NET, with the
        // probability that the measure is above zero; TCI of each part
        if `osam' {
            tempname Sm Sp
            local o1 = `osam'
            local o2 = `osam' + 3*`k'
            mata: st_matrix("`Sm'", (cn__BP[1, (`o1'..`o2')]', cn__BB[1, (`o1'..`o2')]', ///
                cn__BC[1, (`o1'..`o2')]', cn__BS[1, (`o1'..`o2')]', cn__BL[1, (`o1'..`o2')]', ///
                cn__BU[1, (`o1'..`o2')]', cn__BZ[1, (`o1'..`o2')]'))
            local o1 = `o2' + 1
            local o2 = `o2' + 2
            mata: st_matrix("`Sp'", (cn__BP[1, (`o1'..`o2')]', cn__BB[1, (`o1'..`o2')]', ///
                cn__BC[1, (`o1'..`o2')]', cn__BS[1, (`o1'..`o2')]', cn__BL[1, (`o1'..`o2')]', ///
                cn__BU[1, (`o1'..`o2')]'))
            local rn SAM
            foreach s in TO FROM NET {
                forvalues i = 1/`k' {
                    local nm : word `i' of `names'
                    local nm = strtoname("`nm'")
                    local rn `rn' `s':`nm'
                }
            }
            matrix rownames `Sm' = `rn'
            matrix colnames `Sm' = estimate bias corrected se lower upper prob
            matrix rownames `Sp' = TCI_pos TCI_neg
            matrix colnames `Sp' = estimate bias corrected se lower upper
            di _n as txt "{p 0 0 2}Spillover asymmetry measure (percent): SAM = TCI of the positive " ///
                "minus TCI of the negative system, and its TO, FROM and NET; prob: " ///
                "probability (percent) that the measure is above zero{p_end}"
            matlist `Sm', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(12)
            forvalues r = 1/2 {
                local pl = cond(`r' == 1, "positive", "negative")
                local e1 = strtrim(string(`Sp'[`r', 3], "%`w2'.`digits'f"))
                local e2 = strtrim(string(`Sp'[`r', 5], "%`w2'.`digits'f"))
                local e3 = strtrim(string(`Sp'[`r', 6], "%`w2'.`digits'f"))
                di as txt "TCI of the `pl' system, bias-corrected: " as res "`e1'" ///
                    as txt ", `level'% interval [" as res "`e2'" as txt ", " as res "`e3'" as txt "]"
            }
            return scalar sam_prob = `Sm'[1, 7]
            return matrix tci_parts = `Sp'
            return matrix sam = `Sm'
        }
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
            // spillover asymmetry measure, with the probability that it is above zero,
            // and the TCI of each part
            if `osam' {
                local samn sam
                if `doall' {
                    foreach s in to from net {
                        forvalues i = 1/`k' {
                            local samn `samn' sam_`s'_`i'
                        }
                    }
                }
                local samn `samn' tci_pos tci_neg
                local j = `osam' - 1
                foreach m of local samn {
                    local ++j
                    mata: st_store(., st_addvar("double", ("`m'", "`m'_bias", "`m'_bc", "`m'_se", "`m'_lo", "`m'_hi")), (cn__BP[., `j'], cn__BB[., `j'], cn__BC[., `j'], cn__BS[., `j'], cn__BL[., `j'], cn__BU[., `j']))
                    if substr("`m'", 1, 3) == "sam" mata: st_store(., st_addvar("double", "`m'_prob"), cn__BZ[., `j'])
                }
                local ml sam
                if `doall' {
                    forvalues i = 1/`k' {
                        local nj : word `i' of `names'
                        local ml `ml' sam_to_`i' sam_from_`i' sam_net_`i'
                        local lb_sam_to_`i' "SAM: TO, `nj'"
                        local lb_sam_from_`i' "SAM: FROM, `nj'"
                        local lb_sam_net_`i' "SAM: NET, `nj'"
                    }
                }
                local lb_sam "SAM"
                foreach m of local ml {
                    label variable `m' "`lb_`m'', estimate"
                    label variable `m'_bias "`lb_`m'', bootstrap bias"
                    label variable `m'_bc "`lb_`m'', bias-corrected"
                    label variable `m'_se "`lb_`m'', bootstrap standard error"
                    label variable `m'_lo "`lb_`m'', lower `level'% limit"
                    label variable `m'_hi "`lb_`m'', upper `level'% limit"
                    label variable `m'_prob "`lb_`m'', probability (%) above zero"
                }
                foreach m in pos neg {
                    local pl = cond("`m'" == "pos", "positive", "negative")
                    label variable tci_`m' "TCI of the `pl' system, estimate"
                    label variable tci_`m'_bias "TCI of the `pl' system, bootstrap bias"
                    label variable tci_`m'_bc "TCI of the `pl' system, bias-corrected"
                    label variable tci_`m'_se "TCI of the `pl' system, bootstrap standard error"
                    label variable tci_`m'_lo "TCI of the `pl' system, lower `level'% limit"
                    label variable tci_`m'_hi "TCI of the `pl' system, upper `level'% limit"
                }
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
                if "`stype'" == "tci" label variable prob "Probability (%) of a rise of at least `detect'% from the previous window"
                else label variable prob "Probability (%) that `slab' rose at least `detect'% from the previous window"
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
        if `osam' {
            quietly frame `frame': count if (sam_lo > 0 & sam_lo < .) | sam_hi < 0
            local nsam = r(N)
            di as txt "{p 0 0 2}Windows whose `level'% interval of the SAM excludes zero: " ///
                as res `nsam' as txt " of `nw'{p_end}"
            return scalar nsam = `nsam'
        }
        if `isdet' {
            local ofs = cond("`stype'" == "tci", "", " of the `slab'")
            di as txt "{p 0 0 2}Windows with a probability of at least `cutoff'% of a rise" ///
                "`ofs' of at least `detect'%: " as res `nflag' "{p_end}"
            if `nflag' > 0 & `nflag' <= 30 {
                frame `frame': list `=cond("`dv'" != "", "`dv'", "`tv'")' prob if flag == 1, noobs clean
            }
            return scalar nflag = `nflag'
        }
        return local frame "`frame'"
        // pairwise and NPDC bands: long frame, one row per window and ordered pair
        if `dopw' | `donp' {
            if "`pwframe'" == "" local pwframe `frame'_pw
            if inlist("`pwframe'", "`c(frame)'", "`fr'", "`frame'") {
                di as err "pwframe(`pwframe') may not be the current, results or bands frame"
                exit 498
            }
            mata: st_local("bexist", strofreal(st_frameexists("`pwframe'")))
            if `bexist' & "`replace'" == "" {
                di as err "frame {bf:`pwframe'} already exists; specify {cmd:replace} or another {cmd:pwframe()}"
                exit 110
            }
            capture frame drop `pwframe'
            frame create `pwframe'
            frame `pwframe' {
                mata: cn_bootpw(tokens("`names'"), `opw', `onp', "`tv'", cn__btv[cn__bw], ///
                    "`dv'", cn__bdv[cn__bw])
                format `tv' `tfmt'
                if "`dv'" != "" format `dv' `dfmt'
                label variable from "transmitting variable"
                label variable to "receiving variable"
                foreach blk in pw npdc {
                    capture confirm variable `blk'
                    if _rc continue
                    local bl = cond("`blk'" == "pw", "Spillover", "Net spillover")
                    label variable `blk' "`bl' from -from- to -to-, estimate"
                    label variable `blk'_bias "`bl', bootstrap bias"
                    label variable `blk'_bc "`bl', bias-corrected"
                    label variable `blk'_se "`bl', bootstrap standard error"
                    label variable `blk'_lo "`bl', lower `level'% limit"
                    label variable `blk'_hi "`bl', upper `level'% limit"
                }
                label data "connectedness: pairwise bootstrap bands (`reps' draws)"
                char _dta[cn_kind] "bootstrap_pairwise"
            }
            di as txt "Pairwise bands saved in frame " as res "`pwframe'" ///
                as txt " (`=`k'*(`k'-1)' ordered pairs per window)"
            return local pwframe "`pwframe'"
        }
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
    if `isev' | `isdet' return local stat = strtrim("`stype' `snames'")
    capture mata: mata drop cn__bsel
    capture mata: mata drop cn__BP cn__BB cn__BC cn__BS cn__BL cn__BU cn__BZ cn__BD cn__BN cn__BST cn__btv cn__bdv cn__bw
    capture mata: mata drop cn__bev cn__bl cn__bpos cn__P cn__i
end

// Ordering sensitivity of clustered identification (Buchwalter, Diebold and Yilmaz
// 2026b, Section 2.5): distribution of the net connectedness of each cluster (cluster
// table, eq. 2) over the orderings, its mean by position, and the paired
// first-versus-last difference (eq. 4), holding the relative order of the other
// clusters fixed; for rolling and TVP-VAR results, the series of the mean paired
// difference and of the dispersion over the orderings
program define Clorder, rclass
    syntax [, FRame(name) REPLACE DIGits(integer 2) ]

    if "`e(clusters)'" == "" {
        di as err "estat clorder requires results of {cmd:connectedness} with clusters()"
        exit 198
    }
    if "`e(clorder)'" == "simultaneous" {
        di as err "{p 0 0 2}clorder(simultaneous) does not depend on the cluster ordering; " ///
            "there is no ordering sensitivity to diagnose{p_end}"
        exit 198
    }
    if "`e(clorder)'" != "average" {
        di as err "estat clorder requires results with clorder(average)"
        exit 198
    }
    local nfix = (e(cl_ncontrols) == 1)
    local C = e(nclusters)
    local Cp = `C' - `nfix'
    if `Cp' < 2 {
        di as err "estat clorder requires at least two clusters to order"
        exit 198
    }
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }
    capture mata: mata describe cn_clpairs()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }
    local cl `e(clusters)'
    local rn
    forvalues c = 1/`C' {
        local rn `rn' `=strtoname("`: word `c' of `cl''")'
    }
    local pn
    forvalues c = `=`nfix'+1'/`C' {
        local pn `pn' `: word `c' of `rn''
    }
    local w2 = `digits' + 8
    local ptxt = cond(`nfix', " (positions among the clusters after the control cluster)", "")
    local fo `e(frame_clord)'
    local static = ("`fo'" != "")
    if `static' {
        mata: st_local("fexist", strofreal(st_frameexists("`fo'")))
        if !`fexist' {
            di as err "frame {bf:`fo'} not found; rerun {cmd:connectedness}"
            exit 111
        }
        frame `fo': local frc : char _dta[cn_clusters]
        if "`frc'" != "`cl'" {
            di as err "frame {bf:`fo'} does not hold the orderings of the current results; rerun {cmd:connectedness}"
            exit 498
        }
        tempname N P D A
        frame `fo' {
            mata: cn__opos = st_data(., "pos_*"); cn__onet = st_data(., "net_*")
        }
        mata: cn_clorddiag(cn__opos, cn__onet, `nfix', "`N'", "`P'", "`D'")
        matrix `A' = e(cl_table)
        matrix `A' = `A'[`=`C'+2', 1..`C']'
        matrix `N' = `A', `N'
        matrix rownames `N' = `rn'
        matrix colnames `N' = averaged mean sd min max
        matrix rownames `P' = `pn'
        local cp
        forvalues p = 1/`Cp' {
            local cp `cp' pos`p'
        }
        matrix colnames `P' = `cp'
        matrix rownames `D' = `pn'
        matrix colnames `D' = mean sd min max pairs
        di _n as txt "{p 0 0 2}Net connectedness of each cluster (cluster table, percent) over the " ///
            as res e(norderings) as txt " orderings: estimate averaged over orderings, mean, " ///
            "standard deviation, minimum and maximum{p_end}"
        matlist `N', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(14)
        di _n as txt "{p 0 0 2}Mean net connectedness by position in the ordering`ptxt'{p_end}"
        matlist `P', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(14)
        di _n as txt "{p 0 0 2}Paired first-versus-last difference in net connectedness, the cluster " ///
            "first minus the cluster last with the relative order of the other clusters held fixed " ///
            "(Buchwalter, Diebold and Yilmaz 2026b, eq. 4){p_end}"
        matlist `D', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(14)
        if "`frame'" != "" {
            if "`frame'" == "`c(frame)'" | "`frame'" == "`fo'" | "`frame'" == "`e(frame)'" {
                di as err "frame(`frame') may not be the current, results or orderings frame"
                exit 498
            }
            mata: st_local("bexist", strofreal(st_frameexists("`frame'")))
            if `bexist' & "`replace'" == "" {
                di as err "frame {bf:`frame'} already exists; specify {cmd:replace} or another {cmd:frame()}"
                exit 110
            }
            capture frame drop `frame'
            frame create `frame'
            frame `frame' {
                mata: cn_clordlong(cn__opos, cn__onet, `nfix', tokens("`cl'"))
                label variable ordering "cluster ordering"
                label variable cluster "cluster"
                label variable position "position in the ordering`=cond(`nfix', " (after the control cluster)", "")'"
                label variable net "net connectedness (cluster table, percent)"
                label variable dfl "first minus last, paired (rows with the cluster first)"
                label data "connectedness: net connectedness by cluster and ordering"
                char _dta[cn_kind] "clorder_long"
            }
            di _n as txt "One observation per ordering and cluster saved in frame " as res "`frame'"
            return local frame "`frame'"
        }
        capture mata: mata drop cn__opos cn__onet
        return matrix dfl = `D'
        return matrix bypos = `P'
        return matrix net = `N'
    }
    else {
        local fr `e(frame)'
        mata: st_local("fexist", strofreal(st_frameexists("`fr'")))
        if !`fexist' {
            di as err "results frame {bf:`fr'} not found; rerun {cmd:connectedness}"
            exit 111
        }
        capture frame `fr': confirm variable cl_dfl_1 cl_nsd_1
        if _rc {
            di as err "the results frame has no ordering diagnostics; rerun {cmd:connectedness}"
            exit 111
        }
        local tv `e(timevar)'
        local dv `e(datevar)'
        if "`dv'" == "" local dv `tv'
        tempname S
        frame `fr' {
            local dfmt : format `dv'
            mata: cn__ot = st_data(., "`dv'"); cn__od = st_data(., "cl_dfl_*"); cn__os = st_data(., "cl_nsd_*")
        }
        mata: cn__oi = J(1, `C', .); for (cn__c=1; cn__c<=`C'; cn__c++) cn__oi[cn__c] = (allof(cn__od[., cn__c], .) ? . : selectindex(cn__od[., cn__c] :== max(cn__od[., cn__c]))[1])
        mata: st_matrix("`S'", (cn_colmeanm(cn__od)', colmin(cn__od)', colmax(cn__od)', cn_colmeanm(cn__os)'))
        matrix rownames `S' = `rn'
        matrix colnames `S' = dfl_mean dfl_min dfl_max nsd_mean
        di _n as txt "{p 0 0 2}Ordering sensitivity over the " as res e(T0) as txt " periods: mean, " ///
            "minimum and maximum over time of the mean paired first-versus-last difference in net " ///
            "connectedness (cl_dfl_c), and mean over time of its standard deviation over the " ///
            as res e(norderings) as txt " orderings (cl_nsd_c)`ptxt'{p_end}"
        matlist `S', format(%`w2'.`digits'f) border(rows) lines(oneline) rowtitle(" ") twidth(14) nodotz
        forvalues c = `=`nfix'+1'/`C' {
            mata: st_local("ix", strofreal(cn__oi[`c']))
            if "`ix'" != "." {
                mata: st_local("dmax", strofreal(cn__ot[`ix']))
                local dmx : display `dfmt' `dmax'
                local nm : word `c' of `cl'
                di as txt "  largest difference for " as res "`nm'" as txt " in period " as res strtrim("`dmx'")
            }
        }
        di as txt "{p 0 0 2}The series are cl_dfl_c and cl_nsd_c in frame " as res "`fr'" ///
            as txt "; for example, {cmd:frame `fr': tsline cl_dfl_*}.{p_end}"
        capture mata: mata drop cn__ot cn__od cn__os cn__oi cn__c
        return matrix time = `S'
    }
    return scalar norderings = e(norderings)
    return scalar ncontrols = `nfix'
end

// k x k matrix bc with a star where the interval (lo, hi) excludes zero; diagonal
// shown as "."; columns in blocks that fit the line size
program define CnStarMat
    args bc lo hi digits
    local k = rowsof(`bc')
    local names : rownames `bc'
    local w = `digits' + 7
    local nb = max(1, floor((c(linesize) - 14)/`w'))
    forvalues c0 = 1(`nb')`k' {
        local c1 = min(`k', `c0' + `nb' - 1)
        local nc = `c1' - `c0' + 1
        di as txt "{hline 13}{c TT}{hline `=`w'*`nc''}"
        di as txt _col(14) "{c |}" _c
        forvalues j = `c0'/`c1' {
            local nm : word `j' of `names'
            di as txt %`=`w'-1's abbrev("`nm'", `w'-2) " " _c
        }
        di
        di as txt "{hline 13}{c +}{hline `=`w'*`nc''}"
        forvalues i = 1/`k' {
            local nm : word `i' of `names'
            di as txt %12s abbrev("`nm'", 12) " {c |}" _c
            forvalues j = `c0'/`c1' {
                if `i' == `j' {
                    di as txt %`=`w'-1's "." " " _c
                }
                else {
                    local s = cond(`lo'[`i', `j'] > 0 | `hi'[`i', `j'] < 0, "*", " ")
                    di as res %`=`w'-1'.`digits'f `bc'[`i', `j'] "`s'" _c
                }
            }
            di
        }
        di as txt "{hline 13}{c BT}{hline `=`w'*`nc''}"
    }
end
