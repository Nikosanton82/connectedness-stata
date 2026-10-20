*! version 0.7.2  20oct2026
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
        NAME(string) SAVing(string) noDRAW COMBopts(string) BANd(integer 0) COMPonent(string) ///
        noSTYle * ]

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
    // style shared with cnnetwork: light background, muted axes and a fixed palette
    // (dark blue, mustard, teal, rust, grey, plum); nostyle keeps the graph scheme
    local sty = ("`style'" == "")
    local bg "252 251 248"
    local pal `""44 95 125" "201 154 46" "74 143 125" "170 78 61" "110 110 110" "123 97 155""'
    local axsty
    local regsty
    local tsty
    local ptsty size(medium)
    if `sty' {
        local axsty labsize(small) labcolor(gs6) grid glcolor(gs14) glwidth(vthin) glpattern(solid)
        local regsty graphregion(color(`"`bg'"') lcolor(`"`bg'"')) ///
            plotregion(color(`"`bg'"') lcolor(`"`bg'"') lstyle(none))
        local tsty , size(medlarge) color(black)
        local ptsty size(medium) color(gs3)
    }
    // daily, monthly and quarterly dates: label the axis by year so long labels are
    // not clipped
    local xfo
    frame `frame': local xfmt : format `x'
    if substr("`xfmt'", 1, 3) == "%td" local xfo format(%tdCCYY)
    if substr("`xfmt'", 1, 3) == "%tm" local xfo format(%tmCCYY)
    if substr("`xfmt'", 1, 3) == "%tq" local xfo format(%tqCCYY)
    // year ticks at the start of each labeled year: at most 8 labels in single graphs,
    // 6 with up to four panels and 4 with more panels
    local xlab
    local xlabp
    if "`xfo'" != "" {
        frame `frame': quietly summarize `x', meanonly
        local xmin = r(min)
        local xmax = r(max)
        local xty = substr("`xfo'", 10, 1)
        _cngraph_ticks, min(`xmin') max(`xmax') type(`xty') maxn(8)
        local xlab xlabel(`r(ticks)', `xfo' `axsty')
        _cngraph_ticks, min(`xmin') max(`xmax') type(`xty') maxn(6)
        local xlabp4 xlabel(`r(ticks)', `xfo' `axsty')
        _cngraph_ticks, min(`xmin') max(`xmax') type(`xty') maxn(4)
        local xlabp xlabel(`r(ticks)', `xfo' `axsty')
    }
    else if `sty' {
        local xlab xlabel(, `axsty')
        local xlabp `xlab'
        local xlabp4 `xlab'
    }
    if `sty' {
        local ysty ylabel(, angle(0) noticks `axsty') yscale(noline) xscale(lcolor(gs11))
        local xlab `xlab' `ysty'
        local xlabp `xlabp' `ysty'
        local xlabp4 `xlabp4' `ysty'
    }
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
            if "`e(model)'" == "lasso" local mn "LASSO-VAR"
            if "`e(model)'" == "ridge" local mn "ridge VAR"
            if "`e(model)'" == "elasticnet" local mn "elastic-net VAR"
            if "`e(model)'" == "tvpvar" local lab0 "TVP-VAR"
            else if "`e(model)'" == "dccgarch" local lab0 "DCC-GARCH"
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
    local dropw = (`nov' > 0)
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
    // net measures are drawn as two areas, positive (dark blue) and negative (mustard)
    local split = (`sty' & `nov' == 0 & inlist("`what'", "net", "npdc"))
    if `split' {
        tempname work
        local dropw 1
        frame `frame': frame put `x' `panels', into(`work')
        frame `work' {
            foreach v of local panels {
                quietly generate double __p_`v' = max(`v', 0) if `v' < .
                quietly generate double __n_`v' = min(`v', 0) if `v' < .
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
        local ns : word count `panels'
        foreach v of local panels {
            local ++c
            local lc
            if `sty' local lc lcolor(`"`: word `=mod(`c'-1, 6)+1' of `pal''"') lwidth(medthin)
            local plots `plots' (line `v' `x', `lc')
            if `single' == 2 local leg `leg' `c' "`: word `c' of `gnames''"
            if `single' == 3 local leg `leg' `c' `"`: word `c' of `blegend''"'
        }
        forvalues o = 1/`nov' {
            local p : word `=`o'+1' of `pats'
            local s 0
            foreach v of local panels {
                local ++c
                local ++s
                local lc
                if `sty' {
                    local ci = cond(`ns' == 1, `o' + 1, `s')
                    local lc lcolor(`"`: word `=mod(`ci'-1, 6)+1' of `pal''"') lwidth(medthin)
                }
                local plots `plots' (line __o`o'_`v' `x', lpattern(`p') `lc')
            }
        }
        if `nov' > 0 {
            local leg 1 "`lab0'"
            forvalues o = 1/`nov' {
                local leg `leg' `=`o'+1' "`ol`o''"
            }
        }
        local legsty
        if `sty' local legsty region(lcolor(none) fcolor(none)) size(small)
        if `"`leg'"' != "" local legopt legend(order(`leg') rows(1) position(6) `legsty')
        else local legopt legend(off)
        frame `work': twoway `plots', ytitle("") xtitle("") title(`"`title'"' `tsty') ///
            `xlab' `legopt' `regsty' `nameopt' `saveopt' `draw' `options'
    }
    else {
        local glist
        local c 0
        if `np' > 4 local xlab `xlabp'
        else local xlab `xlabp4'
        foreach v of local panels {
            local ++c
            local pt : word `c' of `ptitles'
            local c1 : word 1 of `pal'
            local c2 : word 2 of `pal'
            if `split' {
                local plots (area __p_`v' `x', base(0) fcolor(`"`c1'"') fintensity(100) ///
                    lcolor(`"`c1'"') lwidth(vthin)) (area __n_`v' `x', base(0) ///
                    fcolor(`"`c2'"') fintensity(100) lcolor(`"`c2'"') lwidth(vthin))
            }
            else if `area' {
                local ac
                if `sty' local ac fcolor(`"`c1'"') fintensity(100) lcolor(`"`c1'"') lwidth(vthin)
                local plots (area `v' `x', base(0) `ac')
            }
            else {
                local lc
                if `sty' local lc lcolor(`"`c1'"') lwidth(medthin)
                local plots (line `v' `x', `lc')
                forvalues o = 1/`nov' {
                    local p : word `=`o'+1' of `pats'
                    local lc
                    if `sty' local lc lcolor(`"`: word `=mod(`o', 6)+1' of `pal''"') lwidth(medthin)
                    local plots `plots' (line __o`o'_`v' `x', lpattern(`p') `lc')
                }
            }
            tempname g`c'
            frame `work': twoway `plots', ytitle("") xtitle("") legend(off) `xlab' ///
                title(`"`pt'"', `ptsty') `regsty' nodraw name(`g`c'', replace) `options'
            local glist `glist' `g`c''
        }
        local noteopt
        if `"`note'"' != "" {
            if `sty' local noteopt note(`"`note'"', size(vsmall) color(gs8))
            else local noteopt note(`"`note'"')
        }
        local cregsty
        if `sty' local cregsty graphregion(color(`"`bg'"') lcolor(`"`bg'"'))
        capture noisily graph combine `glist', title(`"`title'"' `tsty') `noteopt' ///
            `cregsty' `nameopt' `saveopt' `draw' `combopts'
        local rc = _rc
        graph drop `glist'
        if `rc' {
            if `dropw' frame drop `work'
            exit `rc'
        }
    }
    if `dropw' frame drop `work'
end

// Tick values at the start of years for a date axis: type d (daily), m (monthly) or
// q (quarterly); the step is the smallest of 1, 2, 5, 10, 20, 25 and 50 years that
// gives at most maxn labels
program define _cngraph_ticks, rclass
    syntax, min(real) max(real) type(string) maxn(integer)
    if "`type'" == "d" {
        local y0 = year(`min')
        local y1 = year(`max')
        if mdy(1, 1, `y0') < `min' local ++y0
    }
    else if "`type'" == "m" {
        local y0 = year(dofm(`min'))
        local y1 = year(dofm(`max'))
        if ym(`y0', 1) < `min' local ++y0
    }
    else {
        local y0 = year(dofq(`min'))
        local y1 = year(dofq(`max'))
        if yq(`y0', 1) < `min' local ++y0
    }
    if `y1' < `y0' {
        return local ticks "#4"
        exit
    }
    foreach st in 1 2 5 10 20 25 50 {
        local step `st'
        if floor(`y1'/`st') - ceil(`y0'/`st') + 1 <= `maxn' continue, break
    }
    local ticks
    forvalues y = `=ceil(`y0'/`step')*`step''(`step')`y1' {
        if "`type'" == "d" local ticks `ticks' `=mdy(1, 1, `y')'
        else if "`type'" == "m" local ticks `ticks' `=ym(`y', 1)'
        else local ticks `ticks' `=yq(`y', 1)'
    }
    if "`ticks'" == "" local ticks "#4"
    return local ticks "`ticks'"
end
