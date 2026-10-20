*! version 0.4.0  20oct2026
*! cnnetwork: network plot of pairwise connectedness after -connectedness-
program define cnnetwork
    version 16.1
    syntax [, Method(string) THReshold(real 0) AT(string) FRame(name) ///
        TItle(string) NAME(string) SAVing(string) noDRAW BANd(integer 0) * ]

    local method = lower("`method'")
    if "`method'" == "" local method npdc
    if !inlist("`method'", "npdc", "pci") {
        di as err "method() must be npdc or pci"
        exit 198
    }
    if `threshold' < 0 | `threshold' >= 1 {
        di as err "threshold() must lie in [0, 1)"
        exit 198
    }
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
        local names : char _dta[cn_names]
        local k     : char _dta[cn_k]
        local tvar  : char _dta[cn_timevar]
        local dvar  : char _dta[cn_datevar]
        local kind  : char _dta[cn_kind]
        local appr  : char _dta[cn_approach]
        local nbands : char _dta[cn_nbands]
        local blabs : char _dta[cn_bands]
    }
    if "`appr'" == "joint" {
        di as err "cnnetwork is not available after approach(joint); pairwise measures are not defined"
        exit 198
    }
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
    if "`names'" == "" {
        di as err "frame {bf:`frame'} was not created by {cmd:connectedness} or {cmd:estat group}"
        exit 198
    }
    if "`kind'" == "group" & "`method'" == "pci" {
        di as err "method(pci) is not available for group measures"
        exit 198
    }

    // one period (at()) or the average over all periods
    local r .
    local when "average over all periods"
    if `"`at'"' != "" {
        local atv = `at'
        frame `frame' {
            local r 0
            forvalues i = 1/`=_N' {
                if `tvar'[`i'] == `atv' local r `i'
                if "`dvar'" != "" {
                    if `dvar'[`i'] == `atv' local r `i'
                }
            }
            if `r' == 0 {
                di as err "at(`at') matches no period in frame `frame'"
                exit 498
            }
            local fmt : format `=cond("`dvar'" != "", "`dvar'", "`tvar'")'
            local when : display `fmt' `=cond("`dvar'" != "", "`dvar'", "`tvar'")'[`r']
            local when = strtrim("`when'")
        }
    }

    capture mata: mata describe cn_netdata()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }
    tempname N E nf
    frame `frame': mata: cn_netdata("`method'", `k', `r', `threshold', "`N'", "`E'", "`suf'")

    frame create `nf'
    frame `nf' {
        quietly svmat double `N', names(n)
        quietly svmat double `E', names(e)
        quietly generate str32 lbl = ""
        forvalues i = 1/`k' {
            quietly replace lbl = `"`: word `i' of `names''"' in `i'
        }
        local plots
        local nplot 0
        local ptype = cond("`method'" == "npdc", "pcarrow", "pcspike")
        forvalues c = 1/5 {
            quietly count if e6 == `c'
            if r(N) {
                local lw = 0.3*`c'
                local plots `plots' (`ptype' e2 e1 e4 e3 if e6 == `c', lwidth(`lw') ///
                    lcolor(gs9) mcolor(gs9) barbsize(`=0.6+0.3*`c'') msize(`=0.8+0.3*`c''))
                local ++nplot
            }
        }
        local thrtxt : display %4.2f `threshold'
        local leg
        quietly count if n6 == 1
        if r(N) {
            local plots `plots' (scatter n2 n1 [aw=n5] if n6 == 1, pstyle(p1) msymbol(O) msize(large))
            local ++nplot
            if "`method'" == "npdc" local leg `leg' `nplot' "Net transmitter"
        }
        quietly count if n6 == 0
        if r(N) {
            local plots `plots' (scatter n2 n1 [aw=n5] if n6 == 0, pstyle(p1) msymbol(O) mfcolor(white) mlwidth(medthick) msize(large))
            local ++nplot
            local leg `leg' `nplot' "Net receiver"
        }
        local plots `plots' (scatter n2 n1, msymbol(none) mlabel(lbl) mlabvposition(n3) ///
            mlabgap(*4) mlabcolor(black) mlabsize(medsmall))
        if `"`leg'"' != "" local legopt legend(order(`leg') rows(1))
        else local legopt legend(off)
        if `"`title'"' == "" {
            local title = cond("`method'" == "npdc", "Net pairwise directional connectedness", ///
                "Pairwise connectedness")
        }
        if `band' {
            local bl : word `band' of `blabs'
            local when "`when'; band `bl' periods"
        }
        if `"`name'"' != ""   local nameopt name(`name')
        if `"`saving'"' != "" local saveopt saving(`saving')
        twoway `plots', aspectratio(1) xscale(off range(-1.7 1.7)) yscale(off range(-1.4 1.4)) ///
            xlabel(none) ylabel(none) plotregion(style(none)) ///
            title(`"`title'"') subtitle("`when'") note("Edge width proportional to strength; threshold `thrtxt'") ///
            `legopt' `nameopt' `saveopt' `draw' `options'
    }
    frame drop `nf'
end
