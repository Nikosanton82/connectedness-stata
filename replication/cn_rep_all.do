* Run every replication do-file and keep a plain-text log of each.
* Usage: from the folder that holds the datasets (after  net get connectedness ),
*        do "<path>/cn_rep_all.do"
* or, within the repository, cd to replication/ and  do cn_rep_all.do
version 16.1
local here "`c(pwd)'"
capture confirm file cn_dy2012.dta
if _rc {
    capture confirm file "../data/cn_dy2012.dta"
    if !_rc quietly cd "../data"
}
local repdir "`here'"
capture confirm file "`repdir'/cn_rep_2009dy.do"
if _rc local repdir "`c(pwd)'"
capture mkdir "`here'/logs"
foreach f in 2009dy 2012dy 2018gg 2018bk 2020acg 2021cg 2021lw 2021bgu 2021cgs 2022adekoya 2022scg 2023cegg 2023cgg 2023balli 2023gcs 2024gghm {
    capture log close rep
    log using "`here'/logs/cn_rep_`f'.log", replace text name(rep)
    do "`repdir'/cn_rep_`f'.do"
    log close rep
}
quietly cd "`here'"
