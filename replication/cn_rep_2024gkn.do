* Replication: Greenwood-Nimmo, Kocenda and Nguyen (2024), detecting statistically
* significant changes in connectedness with a bootstrap. Daily volatilities of 19
* equity markets (Diebold and Yilmaz 2009; data from the authors' replication files,
* CC BY 4.0), rolling VAR(2) on windows of 200 days (window(202): 200 rows of data
* with two lags), horizon 10. For each event, the probability (percent) that the
* bias-corrected spillover index exceeds its level the day before, on the event day
* and 1, 5, 10 and 22 days later, as in Table 3 of the paper (generalized and
* orthogonalized FEVD; 1,000 draws per round). Event 12a (15 Nov 2003) fell on a
* Saturday; the paper's figures correspond to the previous trading day, 14 Nov 2003,
* which is used here (given a non-trading date, estat bootstrap uses the first window
* on or after it, here 17 Nov). Bootstrap probabilities differ from the paper's by
* Monte Carlo error: in Stata 19, 96 of the 100 probabilities of each table were
* within 5 percentage points of the paper's (mean difference 1.5 and 1.3), with the
* same 90 percent verdict for 99 and 100 of them. This do-file takes about 10 to 20
* minutes.
version 16.1
use cn_gkn2024, clear
tsset t
local ev 17oct1997 27oct1997 13jul1998 25aug1998 13jan1999 05jan2000 14apr2000 ///
    23jan2001 20mar2001 11sep2001 19jul2002 14nov2003 20nov2003 30jun2004 22feb2005 ///
    07jul2005 12may2006 19dec2006 08feb2007 09aug2007
connectedness US-TUR, lags(2) horizon(10) window(202) datevar(date) frame(gkn_gvd) replace notable
cngraph tci
estat bootstrap, event(`ev') leads(0 1 5 10 22) seed(1357) nodots
connectedness US-TUR, lags(2) horizon(10) window(202) orthogonalized datevar(date) ///
    frame(gkn_ovd) replace notable
estat bootstrap, event(`ev') leads(0 1 5 10 22) seed(1357) nodots
