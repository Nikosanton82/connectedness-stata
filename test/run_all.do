* run_all.do: certification suite for -connectedness-
* Usage, from Stata:
*     cd "<repository>/test"
*     do run_all.do
* The log is written to test/logs/cert_connectedness.log

clear all
version 16.1
capture log close _all
capture mkdir logs
log using "logs/cert_connectedness.log", replace text name(cert)
display "Stata " c(stata_version) " (" c(flavor) "), " c(current_date) " " c(current_time)
adopath ++ "../ado"
discard
mata: mata clear
capture noisily do cert_connectedness.do
local rc = _rc
log close cert
exit `rc'
