* run_install.do: installs the package from this repository into a temporary
* ado directory (not PLUS), fetches the datasets with net get into a temporary
* data directory, and runs every replication do-file against the installed copy,
* as a user would. Usage, from the test folder:  do run_install.do
clear all
version 16.1
capture log close _all
capture mkdir logs
log using "logs/install_test.log", replace text name(inst)
quietly cd ..
local root "`c(pwd)'"
local tmp  "`c(tmpdir)'cn_install_test"
capture mkdir "`tmp'"
capture mkdir "`tmp'/data"
net set ado "`tmp'"
net install connectedness, from("`root'") replace
net set ado PLUS
adopath ++ "`tmp'"
which connectedness
which cnnetwork
which cnportfolio
quietly cd "`tmp'/data"
net get connectedness, from("`root'") replace
dir cn_*.dta
set graphics off
local rc 0
foreach f in 2009dy 2012dy 2018gg 2019ccm 2020acg 2021cg 2022adekoya 2022scg 2022cagt 2023cegg 2024cgp {
    display as text _n "{hline 78}" _n "Replication `f'" _n "{hline 78}"
    capture noisily do "`root'/replication/cn_rep_`f'.do"
    if _rc {
        local rc = _rc
        display as error "replication `f' failed with r(`rc')"
    }
    else display as text "replication `f': completed"
}
set graphics on
adopath - "`tmp'"
quietly cd "`root'/test"
display as text "Replication suite return code: " as result `rc'
log close inst
exit `rc'
