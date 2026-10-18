* Replication: Gabauer and Gupta (2018), TVP-VAR connectedness decomposition
* Mirrors https://gabauerdavid.github.io/ConnectednessApproach/2018GabauerGupta
version 16.1
use cn_gg2018, clear
tsset month
connectedness US_FPU-JP_CPU, model(tvpvar) lags(1) horizon(10) window(200) ///
    kappa1(0.99) kappa2(0.99) prior(bayes) datevar(date) frame(gg2018) replace
local G `""US: US_FPU US_MPU US_TPU US_CPU" "JP: JP_FPU JP_MPU JP_TPU JP_CPU""'
estat group, type(internal) groups(`G') frame(gg2018_int) replace
cngraph gtci, frame(gg2018_int)
cngraph net, frame(gg2018_int)
cngraph npdc, frame(gg2018_int) select(US_FPU)
estat group, type(external) groups(`G') frame(gg2018_ext) replace
cngraph gtci, frame(gg2018_ext)
