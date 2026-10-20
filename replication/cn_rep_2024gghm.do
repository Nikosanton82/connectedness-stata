* Replication: Gabauer, Gupta, Marfatia and Miller (2024), LASSO, ridge and elastic-net
* VAR connectedness of U.S. state housing prices (data of the R package's gghm2022)
* The R package draws the cross-validation folds at random; set the seed for
* reproducible results.
version 16.1
use cn_gghm2022, clear
tsset t
foreach m in lasso ridge {
    connectedness AK-WY, model(`m') lags(1) horizon(10) seed(2024) nofevd notable ///
        datevar(date) frame(gghm_`m') replace
}
cnnetwork, threshold(0.5)
