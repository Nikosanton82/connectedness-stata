*! version 0.7.1  20oct2026
*! cnnetwork: network plot of pairwise connectedness after -connectedness-
program define cnnetwork
    version 16.1
    syntax [, Method(string) THReshold(real 0) AT(string) FRame(name) ///
        TItle(string) NAME(string) SAVing(string) noDRAW BANd(integer 0) COMPonent(string) ///
        STRaight CURVature(real 0.2) COLors(string asis) LABels(string) noNOTE * ]

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
    if `curvature' < -1 | `curvature' > 1 {
        di as err "curvature() must lie in [-1, 1]"
        exit 198
    }
    if "`straight'" != "" local curvature 0
    local labels = lower("`labels'")
    if "`labels'" == "" local labels auto
    if !inlist("`labels'", "auto", "inside", "outside") {
        di as err "labels() must be auto, inside or outside"
        exit 198
    }
    // node colors: net transmitters, net receivers
    local col1 "44 95 125"
    local col2 "201 154 46"
    if `"`colors'"' != "" {
        gettoken c1 rest : colors
        gettoken c2 rest : rest
        if `"`c2'"' == "" | strtrim(`"`rest'"') != "" {
            di as err `"colors() requires two colors, for example colors(navy orange) or colors("44 95 125" "201 154 46")"'
            exit 198
        }
        local col1 `"`c1'"'
        local col2 `"`c2'"'
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
        local comps : char _dta[cn_comps]
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
    frame `frame' {
        local dv = cond("`dvar'" != "", "`dvar'", "`tvar'")
        if "`dv'" != "" {
            capture confirm variable `dv'
            if !_rc {
                local fmt : format `dv'
                local first : display `fmt' `dv'[1]
                local last  : display `fmt' `dv'[_N]
                if _N == 1 local when "full-sample estimates"
                else local when = "average over " + strtrim("`first'") + " - " + strtrim("`last'")
            }
        }
    }
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
            local fmt : format `dv'
            local when : display `fmt' `dv'[`r']
            local when = strtrim("`when'")
        }
    }
    if `band' {
        local bl : word `band' of `blabs'
        local when "`when'; band `bl' periods"
    }
    if inlist("`suf'", "_c", "_l") local when "`when'; `=cond("`suf'" == "_c", "contemporaneous", "lagged")'"

    // labels inside the nodes for a few short names, otherwise outside the circle
    local maxlen 0
    foreach v of local names {
        local maxlen = max(`maxlen', udstrlen("`v'"))
    }
    if "`labels'" == "auto" local labels = cond(`k' <= 12 & `maxlen' <= 5, "inside", "outside")
    local lout = ("`labels'" == "outside")

    // sizes that shrink with the number of variables. Text sizes are in percent of
    // the graph height (Stata's small = 2.75), which is about the plot height
    // divided by 0.73; a character is about 0.6 times as wide as it is high.
    local sw   = cond(`k' <= 8, 1, cond(`k' <= 16, 0.7, cond(`k' <= 30, 0.5, 0.35)))
    local isz  = cond(`k' <= 8, "small", "vsmall")
    local osz  = cond(`k' <= 8, "medsmall", cond(`k' <= 16, "small", cond(`k' <= 30, "vsmall", "tiny")))
    local opct = cond(`k' <= 8, 3.05, cond(`k' <= 16, 2.75, cond(`k' <= 30, 2.1, 1.75)))
    local rmax = 1.2*min(0.15, 0.42*sin(c(pi)/`k'))
    // half widths of the plot: Lx leaves room for the labels at the sides, Ly for
    // one line of text above and below; the aspect ratio keeps the nodes round
    local leg  = ("`method'" == "npdc")
    local Lx   = 1 + `rmax' + 0.08
    local Ly   = `Lx'
    forvalues it = 1/3 {
        local gh  = (2*`Ly' + 0.3*`leg')/0.73
        local cw  = 0.62*`opct'*`gh'/100
        local th  = `opct'*`gh'/100
        local Lx  = 1 + `rmax' + 0.08 + `lout'*(`cw'*`maxlen' + 0.14)
        local Ly  = 1 + `rmax' + 0.08 + `lout'*(1.3*`th' + 0.03)
    }
    local gh   = (2*`Ly' + 0.3*`leg')/0.73
    local ybot = -`Ly' - 0.3*`leg'
    local asp  = (`Ly' - `ybot')/(2*`Lx')

    capture mata: mata describe cn_netplot()
    if _rc {
        quietly findfile cn_core.mata
        quietly do `"`r(fn)'"'
    }
    tempname N C A D nf
    frame `frame': mata: cn_netplot("`method'", `k', `r', `threshold', "`suf'", `curvature', ///
        `lout', `sw', `gh', "`N'", "`C'", "`A'", "`D'")

    frame create `nf'
    frame `nf' {
        quietly svmat double `N', names(n)
        quietly svmat double `C', names(c)
        quietly svmat double `A', names(a)
        quietly svmat double `D', names(d)
        if _N < 90 quietly set obs 90
        quietly generate str32 lbl = ""
        forvalues i = 1/`k' {
            quietly replace lbl = `"`: word `i' of `names''"' in `i'
        }
        // outside labels are centered on a point beyond the node, further out for
        // longer labels at the sides
        if `lout' {
            tempvar ang rho
            quietly generate double `ang' = atan2(n2, n1) in 1/`k'
            quietly generate double `rho' = 1 + n3 + 0.035 + ///
                (0.08 + 0.5*`cw'*udstrlen(lbl))*abs(cos(`ang')) + 0.6*`opct'*`gh'/100*abs(sin(`ang')) in 1/`k'
            quietly replace n5 = `rho'*cos(`ang') in 1/`k'
            quietly replace n6 = `rho'*sin(`ang') in 1/`k'
        }
        local plots
        // edges, weakest first; width and shade grow with strength
        forvalues q = 1/6 {
            local lw = `sw'*(0.15 + 0.3*`q')
            local g  = 12 - `q'
            quietly count if c3 == `q' & c1 < .
            if r(N) {
                local plots `plots' (line c2 c1 if c3 == `q', cmissing(n) lcolor(gs`g') lwidth(`lw'))
            }
            if "`method'" == "npdc" {
                quietly count if a5 == `q'
                if r(N) {
                    local ms = max(0.9, `sw'*(1.0 + 0.35*`q'))
                    local plots `plots' (pcarrow a2 a1 a4 a3 if a5 == `q', lcolor(gs`g') ///
                        mcolor(gs`g') lwidth(vthin) mlwidth(vthin) msize(`ms') barbsize(`ms'))
                }
            }
        }
        // nodes
        quietly count if d4 == 1 & d1 < .
        if r(N) local plots `plots' (rarea d2 d3 d1 if d4 == 1, cmissing(n) ///
            fcolor(`"`col1'"') fintensity(100) lcolor(`"`col1'"') lwidth(vthin))
        quietly count if d4 == 0 & d1 < .
        if r(N) local plots `plots' (rarea d2 d3 d1 if d4 == 0, cmissing(n) ///
            fcolor(`"`col2'"') fintensity(100) lcolor(`"`col2'"') lwidth(vthin))
        // labels
        if `lout' local plots `plots' (scatter n6 n5 in 1/`k', msymbol(none) mlabel(lbl) ///
            mlabposition(0) mlabcolor(gs3) mlabsize(`osz'))
        else local plots `plots' (scatter n2 n1 in 1/`k', msymbol(none) mlabel(lbl) ///
            mlabposition(0) mlabcolor(white) mlabsize(`isz'))
        // legend drawn in the plot, centered below the circle: a disc and a text for
        // each node type
        local legtext
        if `leg' {
            local y0  = -`Ly' - 0.16
            local rl  = 0.045*`gh'/3.5
            local lcw = 0.55*2.75*`gh'/100
            local w1  = 2*`rl' + 0.03 + 15*`lcw'
            local w2  = 2*`rl' + 0.03 + 12*`lcw'
            local sp  = 0.12*`gh'/3.5
            local xa  = -(`w1' + `sp' + `w2')/2 + `rl'
            local xb  = `xa' + `w1' + `sp'
            quietly generate double lgx = .
            quietly generate double lgu = .
            quietly generate double lgl = .
            quietly generate byte   lgt = .
            local row 0
            foreach t in 1 0 {
                local x0 = cond(`t', `xa', `xb')
                forvalues m = 0/40 {
                    local ++row
                    local ph = c(pi)*`m'/40
                    quietly replace lgx = `x0' - `rl'*cos(`ph') in `row'
                    quietly replace lgu = `y0' + `rl'*sin(`ph') in `row'
                    quietly replace lgl = `y0' - `rl'*sin(`ph') in `row'
                    quietly replace lgt = `t' in `row'
                }
                local ++row
                quietly replace lgt = `t' in `row'
            }
            local plots `plots' (rarea lgu lgl lgx if lgt == 1, cmissing(n) fcolor(`"`col1'"') ///
                fintensity(100) lcolor(`"`col1'"') lwidth(vthin)) ///
                (rarea lgu lgl lgx if lgt == 0, cmissing(n) fcolor(`"`col2'"') ///
                fintensity(100) lcolor(`"`col2'"') lwidth(vthin))
            local t1 = `xa' + `rl' + 0.03
            local t2 = `xb' + `rl' + 0.03
            local legtext text(`y0' `t1' "Net transmitter", place(e) size(small) color(gs4)) ///
                text(`y0' `t2' "Net receiver", place(e) size(small) color(gs4))
        }
        if `"`title'"' == "" {
            local title = cond("`method'" == "npdc", "Net pairwise directional connectedness", ///
                "Pairwise connectedness")
        }
        local thrtxt : display %4.2f `threshold'
        local notetxt "Edge width and shade proportional to strength"
        if `threshold' > 0 local notetxt "`notetxt'; edges at or below `=strtrim("`thrtxt'")' of the largest omitted"
        local noteopt note("`notetxt'", size(vsmall) color(gs8))
        if "`note'" != "" local noteopt
        if `"`name'"' != ""   local nameopt name(`name')
        if `"`saving'"' != "" local saveopt saving(`saving')
        local bg "252 251 248"
        twoway `plots', aspectratio(`asp') ///
            xscale(off range(-`Lx' `Lx')) yscale(off range(`ybot' `Ly')) xlabel(none) ylabel(none) ///
            legend(off) `legtext' ///
            graphregion(color(`"`bg'"') lcolor(`"`bg'"')) ///
            plotregion(color(`"`bg'"') lcolor(`"`bg'"') lstyle(none) margin(zero)) ///
            title(`"`title'"', size(medlarge) color(black)) ///
            subtitle(`"`when'"', size(small) color(gs7)) `noteopt' ///
            `nameopt' `saveopt' `draw' `options'
    }
    frame drop `nf'
end
