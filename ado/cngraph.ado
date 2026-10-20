*! version 0.5.0  20oct2026
*! cngraph: graphs of dynamic connectedness measures after -connectedness-
program define cngraph
    version 16.1
    gettoken what 0 : 0, parse(" ,")
    local what = lower("`what'")
    local okm tci to from net npt npdc pci inf gtci tciext
    if !`: list what in okm' {
        di as err "cngraph: specify one of `okm'"
        exit 198
    }
    syntax [, FRame(name) COMPare(namelist) SELect(string) RAW TItle(string) ///
        NAME(string) SAVing(string) noDRAW COMBopts(string) BANd(integer 0) COMPonent(string) * ]

    local utitle `"`title'"'

    // ------------------------------------------------------------- base frame
    if "`frame'" == "" {
        if "`e(cmd)'" != "connectedness" error 301
        local frame `e(frame)'
    }
    mata: st_local("fexist", strofreal(st_frameexists("`frame'")))
    if !`fexist' {
        di as err "frame {bf:`frame'} not found"
        exit 111
    }
    frame `frame' {
        local names  : char _dta[cn_names]
        local k      : char _dta[cn_k]
        local tvar   : char _dta[cn_timevar]
        local dvar   : char _dta[cn_datevar]
        local corr   : char _dta[cn_corrected]
        local kind   : char _dta[cn_kind]
        local gtype  : char _dta[cn_gtype]
        local gnames : char _dta[cn_gnames]
        local appr   : char _dta[cn_approach]
        local nbands : char _dta[cn_nbands]
        local blabs  : char _dta[cn_bands]
        local comps  : char _dta[cn_comps]
    }
    if "`appr'" == "" local appr time
    if "`nbands'" == "" local nbands 0
    local suf
    if `band' != 0 {
        if "`appr'" != "frequency" {
            di as err "band() requires a frame estimated with approach(frequency)"
            exit 198
        }
        if `band' < 0 | `band' > `nbands' {
            di as err "band() must be between 1 and `nbands'"
            exit 198
        }
        local suf _f`band'
    }
    local component = lower("`component'")
    if "`component'" != "" {
        if "`comps'" != "c l" {
            di as err "component() requires a frame estimated with approach(r2) and lags() > 0"
            exit 198
        }
        if `band' {
            di as err "options band() and component() may not be combined"
            exit 184
        }
        if inlist("`component'", "c", "cont", "contemporaneous") local suf _c
        else if inlist("`component'", "l", "lag", "lagged") local suf _l
        else {
            di as err "component() must be contemporaneous or lagged"
            exit 198
        }
    }
    if "`appr'" == "joint" & inlist("`what'", "npt", "npdc", "pci", "inf") {
        di as err "`what' is not available after approach(joint)"
        exit 198
    }
    if "`names'" == "" {
        di as err "frame {bf:`frame'} was not created by {cmd:connectedness} or {cmd:estat group}"
        exit 198
    }
    local x = cond("`dvar'" != "", "`dvar'", "`tvar'")
    // daily dates: label the axis by year so long labels are not clipped
    local xlab
    frame `frame': local xfmt : format `x'
    if substr("`xfmt'", 1, 3) == "%td" local xlab xlabel(#6, format(%tdCCYY))
    local usec = ("`corr'" != "" & "`raw'" == "")
    if `"`name'"' != ""   local nameopt name(`name')
    if `"`saving'"' != "" local saveopt saving(`saving')

    // ------------------------------------------------------------- overlays
    local nov 0
    local lab0 ""
    if "`compare'" != "" {
        tempname hold
        _estimates hold `hold', copy restore
        foreach s of local compare {
            quietly estimates restore `s'
            if "`e(cmd)'" != "connectedness" {
                di as err "stored estimates `s' are not from {cmd:connectedness}"
                exit 322
            }
            local ++nov
            local of`nov' `e(frame)'
            local ol`nov' `s'
        }
        _estimates unhold `hold'
        local lab0 "current"
        if "`frame'" == "`e(frame)'" {
            local mn "VAR"
            if "`e(model)'" == "qvar" local mn "QVAR, tau `e(quantiles)'"
            if "`e(model)'" == "qvar" & "`e(quantile)'" != "" local mn "QVAR, tau `=e(quantile)'"
            if "`e(model)'" == "lad" local mn "LAD-VAR"
            if "`e(model)'" == "tvpvar" local lab0 "TVP-VAR"
            else if e(window) > 0 local lab0 "`mn', window `e(window)'"
            else local lab0 "`mn'"
        }
    }
    else if "`frame'" == "`e(frame)'" & "`e(asymmetric)'" != "" {
        local nov 2
        local of1 `e(frame_pos)'
        local ol1 positive
        local of2 `e(frame_neg)'
        local ol2 negative
        local lab0 "all"
    }

    // ------------------------------------------------------------- what to plot
    local sel 0
    if `"`select'"' != "" {
        local sel : list posof `"`select'"' in names
        if `sel' == 0 {
            di as err `"select(`select') is not one of: `names'"'
            exit 198
        }
    }
    local panels
    local ptitles
    if "`what'" == "tci" {
        local panels = cond(`usec', "ctci", "tci") + "`suf'"
        local single 1
        if `"`title'"' == "" local title = cond(`usec', "Corrected total connectedness index", ///
            "Total connectedness index")
        if "`appr'" == "joint" & `"`title'"' == "Total connectedness index" {
            local title "Joint total connectedness index"
        }
        // frequency results: the total and each band in one graph
        if "`appr'" == "frequency" & `band' == 0 & `nov' == 0 {
            local single 3
            local blegend `""Total""'
            forvalues b = 1/`nbands' {
                local panels `panels' `=cond(`usec', "ctci", "tci")'_f`b'
                local blegend `"`blegend' "`: word `b' of `blabs''""'
            }
        }
        if "`appr'" == "r2" & "`comps'" == "c l" & "`suf'" == "" & `nov' == 0 {
            local single 3
            local blegend `""Overall" "Contemporaneous" "Lagged""'
            local panels `panels' `=cond(`usec', "ctci", "tci")'_c `=cond(`usec', "ctci", "tci")'_l
        }
    }
    else if "`what'" == "tciext" {
        if "`gtype'" != "aggregated" {
            di as err "tciext is available only for frames from {cmd:estat group, type(aggregated)}"
            exit 198
        }
        local panels = cond(`usec', "ctciext", "tciext")
        local single 1
        if `"`title'"' == "" local title "External TCI (TCI_ext of ConnectednessApproach)"
    }
    else if "`what'" == "gtci" {
        if !inlist("`gtype'", "internal", "external") {
            di as err "gtci is available only for frames from {cmd:estat group, type(internal)} or type(external)"
            exit 198
        }
        local m : word count `gnames'
        forvalues g = 1/`m' {
            local panels `panels' `=cond(`usec', "cgtci", "gtci")'_`g'
        }
        local single 2
        if `"`title'"' == "" local title "Group total connectedness"
    }
    else if inlist("`what'", "to", "from", "net", "npt") {
        local single 0
        forvalues i = 1/`k' {
            if `sel' == 0 | `sel' == `i' {
                local panels `panels' `what'_`i'`suf'
                local ptitles `"`ptitles' `"`: word `i' of `names''"'"'
            }
        }
        if `"`title'"' == "" {
            if "`what'" == "to"   local title "Total directional connectedness TO others"
            if "`what'" == "from" local title "Total directional connectedness FROM others"
            if "`what'" == "net"  local title "Net total directional connectedness"
            if "`what'" == "npt"  local title "Net pairwise transmission"
        }
    }
    else {
        local single 0
        if "`kind'" == "group" & "`what'" != "npdc" {
            di as err "`what' is not available for group measures"
            exit 198
        }
        local sep = cond("`what'" == "npdc", " - ", ", ")
        forvalues i = 1/`=`k'-1' {
            forvalues j = `=`i'+1'/`k' {
                if `sel' == 0 | `sel' == `i' | `sel' == `j' {
                    local panels `panels' `what'_`i'_`j'`suf'
                    local ptitles `"`ptitles' `"`: word `i' of `names''`sep'`: word `j' of `names''"'"'
                }
            }
        }
        if `"`title'"' == "" {
            if "`what'" == "npdc" local title "Net pairwise directional connectedness"
            if "`what'" == "pci"  local title "Pairwise connectedness index"
            if "`what'" == "inf"  local title "Pairwise influence index"
        }
    }
    // in a block, so that the word of the band list is not looked up when band is 0
    if `band' & `"`utitle'"' == "" {
        local bl : word `band' of `blabs'
        local title `"`title', band `bl' periods"'
    }
    if inlist("`suf'", "_c", "_l") & `"`utitle'"' == "" {
        local title `"`title', `=cond("`suf'" == "_c", "contemporaneous", "lagged")'"'
    }
    local np : word count `panels'
    if `np' > 36 {
        di as err "`np' panels requested; use select() to focus on one variable"
        exit 198
    }

    // ------------------------------------------------------------- working frame
    local work `frame'
    if `nov' > 0 {
        tempname work
        frame copy `frame' `work'
        frame `work' {
            forvalues o = 1/`nov' {
                tempvar lnk
                quietly frlink 1:1 `tvar', frame(`of`o'') generate(`lnk')
                foreach v of local panels {
                    quietly frget __o`o'_`v' = `v', from(`lnk')
                }
            }
        }
    }
    local pats solid dash shortdash longdash dot dash_dot
    local note
    if `nov' > 0 {
        local note "`lab0': solid"
        forvalues o = 1/`nov' {
            local note "`note'; `ol`o'': `: word `=`o'+1' of `pats''"
        }
    }

    // ------------------------------------------------------------- draw
    local area = (`nov' == 0 & inlist("`what'", "to", "from", "net", "npdc"))
    if `single' {
        local plots
        local leg
        local c 0
        foreach v of local panels {
            local ++c
            local plots `plots' (line `v' `x')
            if `single' == 2 local leg `leg' `c' "`: word `c' of `gnames''"
            if `single' == 3 local leg `leg' `c' `"`: word `c' of `blegend''"'
        }
        forvalues o = 1/`nov' {
            local p : word `=`o'+1' of `pats'
            foreach v of local panels {
                local ++c
                local plots `plots' (line __o`o'_`v' `x', lpattern(`p'))
            }
        }
        if `nov' > 0 {
            local leg 1 "`lab0'"
            forvalues o = 1/`nov' {
                local leg `leg' `=`o'+1' "`ol`o''"
            }
        }
        if `"`leg'"' != "" local legopt legend(order(`leg') rows(1) position(6))
        else local legopt legend(off)
        frame `work': twoway `plots', ytitle("") xtitle("") title(`"`title'"') ///
            `xlab' `legopt' `nameopt' `saveopt' `draw' `options'
    }
    else {
        local glist
        local c 0
        foreach v of local panels {
            local ++c
            local pt : word `c' of `ptitles'
            if `area' {
                local plots (area `v' `x', base(0))
            }
            else {
                local plots (line `v' `x')
                forvalues o = 1/`nov' {
                    local p : word `=`o'+1' of `pats'
                    local plots `plots' (line __o`o'_`v' `x', lpattern(`p'))
                }
            }
            tempname g`c'
            frame `work': twoway `plots', ytitle("") xtitle("") legend(off) `xlab' ///
                title(`"`pt'"', size(medium)) nodraw name(`g`c'', replace) `options'
            local glist `glist' `g`c''
        }
        local noteopt
        if `"`note'"' != "" local noteopt note(`"`note'"')
        capture noisily graph combine `glist', title(`"`title'"') `noteopt' ///
            `nameopt' `saveopt' `draw' `combopts'
        local rc = _rc
        graph drop `glist'
        if `rc' {
            if `nov' > 0 frame drop `work'
            exit `rc'
        }
    }
    if `nov' > 0 frame drop `work'
end
