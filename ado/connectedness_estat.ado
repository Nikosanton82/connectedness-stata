*! version 0.7.0  20oct2026
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
