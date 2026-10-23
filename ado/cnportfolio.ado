*! version 0.8.0  23oct2026
*! cnportfolio: portfolio weights, hedge ratios and hedging effectiveness after -connectedness-
*! Port of the portfolio functions of the ConnectednessApproach R package (D. Gabauer, 1.0.4)
*! Author: Nikolaos Antonakakis.  License: GPL-3

program define cnportfolio, rclass
    version 16.1
    syntax varlist(numeric min=2) [if] [in], Method(string)          ///
        [ MATrix(string) FRame(string) TEST(string) CUMulative(string) ///
          noLONG SCale(real 0) SAVing(string) DIGits(integer 2) noTABle ]

    // ---------------------------------------------------------- method
    local method = lower(strtrim(`"`method'"'))
    if inlist("`method'", "mcp", "minimum", "minimumconnectedness")      local method mcp
    else if inlist("`method'", "mhp", "multivariate", "hedging")        local method mhp
    else if inlist("`method'", "rpp", "riskparity", "parity")            local method rpp
    else if inlist("`method'", "ewp", "equal", "equallyweighted")        local method ewp
    else if inlist("`method'", "bivariate", "bpp", "pairwise", "pairs")  local method bivariate
    else if inlist("`method'", "hedge", "hr", "hedgeratio", "ratio")     local method hedge
    else {
        di as err "method() must be {cmd:mcp}, {cmd:mhp}, {cmd:rpp}, {cmd:ewp}, " ///
            "{cmd:bivariate} or {cmd:hedge}"
        exit 198
    }
    if `digits' < 0 | `digits' > 8 {
        di as err "digits() must be between 0 and 8"
        exit 198
    }

    // ---------------------------------------------------------- results frame
    if `"`frame'"' == "" {
        if "`e(cmd)'" == "connectedness" local frame `e(frame)'
        else                             local frame cn_results
    }
    mata: st_local("fexist", strofreal(st_frameexists("`frame'")))
    if !`fexist' {
        di as err "results frame {bf:`frame'} not found; run {cmd:connectedness} first or specify {cmd:frame()}"
        exit 111
    }
    frame `frame' {
        local kind     : char _dta[cn_kind]
        local k        : char _dta[cn_k]
        local names    : char _dta[cn_names]
        local timevar  : char _dta[cn_timevar]
        local datevar  : char _dta[cn_datevar]
        local hascov   : char _dta[cn_cov]
        local approach : char _dta[cn_approach]
        local nres = _N
    }
    if "`kind'" != "main" {
        di as err "frame {bf:`frame'} does not hold the results of {cmd:connectedness}"
        exit 498
    }
    local nv : word count `varlist'
    if `nv' != `k' {
        di as err "`nv' return variables specified, but the results in frame `frame' are for `k' variables (`names')"
        exit 198
    }

    // ---------------------------------------------------------- matrix used
    local matrix = lower(strtrim(`"`matrix'"'))
    if "`method'" == "ewp" local matrix none
    else if "`matrix'" == "" local matrix = cond("`method'" == "mcp", "pci", "cov")
    else if inlist("`matrix'", "covariance", "variance") local matrix cov
    else if inlist("`matrix'", "correlation", "corr")    local matrix cor
    if !inlist("`matrix'", "pci", "cov", "cor", "none") {
        di as err "matrix() must be {cmd:pci}, {cmd:cov} or {cmd:cor}"
        exit 198
    }
    local hvars
    if inlist("`matrix'", "cov", "cor") {
        if "`hascov'" != "1" {
            di as err "the results frame holds no covariances;"
            if inlist("`approach'", "r2", "modelfree") ///
                di as err "approach(`approach') does not estimate a covariance matrix; use {cmd:matrix(pci)}"
            else di as err "rerun {cmd:connectedness} without {cmd:nocovariance}, or use {cmd:matrix(pci)}"
            exit 498
        }
        forvalues i = 1/`k' {
            forvalues j = `i'/`k' {
                local hvars `hvars' cov_`i'_`j'
            }
        }
    }
    else if "`matrix'" == "pci" {
        forvalues i = 1/`=`k'-1' {
            forvalues j = `=`i'+1'/`k' {
                local hvars `hvars' pci_`i'_`j'
            }
        }
        frame `frame': capture confirm numeric variable pci_1_2, exact
        if _rc {
            di as err "the results frame holds no pairwise connectedness index (PCI)"
            exit 498
        }
    }

    // ---------------------------------------------------------- options
    local test = lower(strtrim(`"`test'"'))
    if inlist("`test'", "", "fisher", "f")                  local test fisher
    else if inlist("`test'", "bartlett")                    local test bartlett
    else if inlist("`test'", "fk", "fligner", "fligner-killeen", "flignerkilleen") local test fk
    else if inlist("`test'", "levene")                      local test levene
    else if inlist("`test'", "bf", "brownforsythe", "brown-forsythe") local test bf
    else {
        di as err "test() must be {cmd:fisher}, {cmd:bartlett}, {cmd:fk}, {cmd:levene} or {cmd:bf}"
        exit 198
    }
    local cumulative = lower(strtrim(`"`cumulative'"'))
    if inlist("`cumulative'", "", "sum", "cumsum")          local cumulative sum
    else if inlist("`cumulative'", "prod", "cumprod", "product", "compound") local cumulative prod
    else {
        di as err "cumulative() must be {cmd:sum} or {cmd:prod}"
        exit 198
    }
    local islong = ("`long'" == "")
    if `scale' < 0 {
        di as err "scale() must be positive"
        exit 198
    }

    local sframe cn_portfolio
    local sreplace replace
    if `"`saving'"' != "" {
        gettoken sframe sopt : saving, parse(",")
        local sframe = strtrim("`sframe'")
        local sopt = strtrim(subinstr("`sopt'", ",", "", 1))
        if !inlist("`sopt'", "", "replace") {
            di as err "saving(): only option replace is allowed"
            exit 198
        }
        local sreplace `sopt'
        confirm name `sframe'
        if inlist("`sframe'", "`c(frame)'", "`frame'") {
            di as err "saving() must name a frame other than the current and the results frame"
            exit 498
        }
        mata: st_local("sexist", strofreal(st_frameexists("`sframe'")))
        if `sexist' & "`sreplace'" == "" {
            di as err "frame {bf:`sframe'} already exists; specify {cmd:saving(`sframe', replace)}"
            exit 110
        }
    }
    if inlist("`sframe'", "`c(frame)'", "`frame'") {
        di as err "saving() must name a frame other than the current and the results frame"
        exit 498
    }

    // ---------------------------------------------------------- align the data
    marksample touse
    capture confirm numeric variable `timevar', exact
    if _rc {
        di as err "the time variable `timevar' of the results is not in the current data"
        exit 111
    }
    local dkeep
    if "`datevar'" != "" {
        capture confirm numeric variable `datevar', exact
        if !_rc local dkeep `datevar'
    }
    tempname work
    frame put `timevar' `dkeep' `varlist' if `touse', into(`work')
    local static = (`nres' == 1)
    frame `work' {
        tempvar link
        local i 0
        foreach v of local varlist {
            local ++i
            rename `v' __r`i'
            local rvars `rvars' __r`i'
        }
        if !`static' & "`matrix'" != "none" {
            quietly frlink 1:1 `timevar', frame(`frame') generate(`link')
            quietly keep if !missing(`link')
            if _N == 0 {
                di as err "no period of the results frame matches the return data"
                exit 2000
            }
            quietly frget `hvars', from(`link')
            drop `link'
        }
        else if "`matrix'" != "none" {
            tempname hv
            foreach h of local hvars {
                frame `frame': scalar `hv' = `h'[1]
                quietly generate double `h' = `hv'
            }
        }
        else if !`static' {
            // equal weights: the periods of the results frame
            quietly frlink 1:1 `timevar', frame(`frame') generate(`link')
            quietly keep if !missing(`link')
            drop `link'
        }
        local miss 0
        foreach v of varlist `rvars' `hvars' {
            quietly count if missing(`v')
            local miss = `miss' + r(N)
        }
        if `miss' {
            di as err "missing values in the returns or in the matrices of the aligned periods"
            exit 416
        }
        sort `timevar'
        local N = _N
        if `N' < 3 {
            di as err "too few periods (`N')"
            exit 2001
        }

        // annualization of the Sharpe ratios: the periodicity of the dates, as
        // PerformanceAnalytics (daily 252, weekly 52, monthly 12, quarterly 4, yearly 1)
        local scnote
        if `scale' == 0 {
            local dv = cond("`dkeep'" != "", "`dkeep'", "`timevar'")
            local df : format `dv'
            if substr("`df'", 1, 3) == "%tm" | substr("`df'", 1, 3) == "%tM" local scale 12
            else if inlist(substr("`df'", 1, 3), "%tq", "%tQ")  local scale 4
            else if inlist(substr("`df'", 1, 3), "%tw", "%tW")  local scale 52
            else if inlist(substr("`df'", 1, 3), "%ty", "%tY")  local scale 1
            else if inlist(substr("`df'", 1, 3), "%th", "%tH")  local scale 2
            else if inlist(substr("`df'", 1, 3), "%td", "%tD") {
                tempvar dd
                quietly generate double `dd' = `dv'[_n] - `dv'[_n-1] in 2/l
                quietly summarize `dd', detail
                local md = r(p50)
                if `md' <= 1                  local scale 252
                else if `md' <= 7             local scale 52
                else if `md' <= 31            local scale 12
                else if `md' <= 92            local scale 4
                else                          local scale 1
                drop `dd'
            }
            else {
                local scale 1
                local scnote "(no date format found: returns and Sharpe ratios not annualized; see scale())"
            }
        }

        // ---------------------------------------------------------- compute
        capture mata: mata describe cn_pfmain()
        if _rc {
            quietly findfile cn_core.mata
            quietly do `"`r(fn)'"'
        }
        tempname tab
        local rowlab
        mata: cn_pfmain("`method'", "`rvars'", "`hvars'", "`matrix'", `islong', "`test'", ///
            "`cumulative'", `scale', "`names'", "`tab'")
        keep `timevar' `dkeep'
        mata: cn_pfpost()
        label data "cnportfolio results (`method')"
        char _dta[cn_kind]   "portfolio"
        char _dta[cn_method] "`method'"
        char _dta[cn_names]  "`names'"
        char _dta[cn_timevar] "`timevar'"
        char _dta[cn_datevar] "`dkeep'"
    }
    capture frame drop `sframe'
    frame rename `work' `sframe'
    frame `sframe': quietly tsset `timevar'

    // ---------------------------------------------------------- table
    local cn "Mean Std_Dev p5 p95 HE p_value Return Risk SR"
    matrix colnames `tab' = `cn'
    local rn
    foreach r of local rowlab {
        local rn `rn' `=subinstr("`r'", "/", "_", .)'
    }
    matrix rownames `tab' = `rn'

    local mtitle = cond("`method'" == "mcp", "Minimum connectedness portfolio",            ///
        cond("`method'" == "mhp", "Multivariate hedging portfolio",                        ///
        cond("`method'" == "rpp", "Risk-parity portfolio",                                 ///
        cond("`method'" == "ewp", "Equally weighted portfolio",                            ///
        cond("`method'" == "bivariate", "Bivariate portfolios (Kroner and Ng 1998)",        ///
        "Hedge ratios (Kroner and Sultan 1993)")))))
    local wtitle = cond("`method'" == "mhp", "Hedge weights",                              ///
        cond("`method'" == "hedge", "Hedge ratio", "Weight"))
    local mlab = cond("`matrix'" == "pci", "pairwise connectedness index",                 ///
        cond("`matrix'" == "cov", "covariance matrix",                                     ///
        cond("`matrix'" == "cor", "correlation matrix", "none")))
    local tlab = cond("`test'" == "fisher", "F test", cond("`test'" == "bartlett", "Bartlett", ///
        cond("`test'" == "fk", "Fligner-Killeen", cond("`test'" == "levene", "Levene", "Brown-Forsythe"))))

    if "`table'" == "" {
        di
        di as txt "`mtitle'"
        di as txt "Matrix: " as res "`mlab'" as txt "   Periods: " as res `N' ///
            as txt cond(`static', " (static estimates)", "") as txt "   Results: " as res "`frame'"
        local wl
        if inlist("`method'", "mcp", "bivariate") local wl = cond(`islong', "long only", "long and short")
        if "`method'" == "rpp" local wl "long only"
        di as txt "Variance test: " as res "`tlab'" as txt "   Annualization: " as res `scale' ///
            as txt "/year" as txt cond("`wl'" != "", "   Weights: ", "") as res "`wl'"
        if "`scnote'" != "" di as txt "`scnote'"
        local lw 12
        foreach r of local rowlab {
            local lw = max(`lw', udstrlen("`r'"))
        }
        local fmt %9.`digits'f
        local h1 = 30 - udstrlen("`wtitle'")
        di as txt "{hline `=`lw'+1'}{c TT}{hline 63}"
        di as txt %`lw's "" " {c |} {hline 3} `wtitle' {hline `h1'} {hline 3} Portfolio {hline 12}"
        di as txt %`lw's "" " {c |}" %9s "Mean" %9s "Std.Dev." %9s "5%" %9s "95%" ///
            %9s "HE" %9s "p-value" %9s "SR"
        di as txt "{hline `=`lw'+1'}{c +}{hline 63}"
        local i 0
        foreach r of local rowlab {
            local ++i
            di as txt %`lw's "`r'" " {c |}" as res                                         ///
                `fmt' `tab'[`i',1] `fmt' `tab'[`i',2] `fmt' `tab'[`i',3] `fmt' `tab'[`i',4] ///
                `fmt' `tab'[`i',5] %9.3f `tab'[`i',6] `fmt' `tab'[`i',9]
        }
        di as txt "{hline `=`lw'+1'}{c BT}{hline 63}"
        di as txt "{p 0 0 2}HE: hedging effectiveness, 1 - Var(portfolio)/Var(asset). p-value: " ///
            "`tlab' of equal variances of the asset and the portfolio. SR: annualized Sharpe "  ///
            "ratio of the portfolio; its annualized return and risk are in r(table). "          ///
            "Series in frame {bf:`sframe'}.{p_end}"
    }

    matrix colnames `tab'_p = fisher bartlett fk levene bf
    matrix rownames `tab'_p = `rn'
    return matrix pvalues = `tab'_p
    return matrix table = `tab'
    return scalar N     = `N'
    return scalar k     = `k'
    return scalar scale = `scale'
    return local method   "`method'"
    return local matrix   "`matrix'"
    return local test     "`test'"
    return local cumulative "`cumulative'"
    return local frame    "`frame'"
    return local saving   "`sframe'"
    return local rows     `"`rowlab'"'
    return local long = cond(`islong', "long", "nolong")
end
