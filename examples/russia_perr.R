#'Russia PERR estimation
#'
require(BVAR)
require(vars)
set.seed(1234567)

data_to_model <- read.csv("data/russian_data.csv")
data_to_model <- data_to_model[,-1]
data_to_model$neer_qoq <- data_to_model$neer_qoq*-1

bvar_model <- bvar(data_to_model, lags=1,
                   priors = bv_priors(hyper="full"),
                   mh = bv_mh(),
                   n_thin=2, n_burn=8000, n_draw=10000)
var_model <- VAR(data_to_model, p = 1, type = "const")

#Zero and sign matrix
horizon <- 2
nvars <- NCOL(data_to_model)
signs_ru <- array(NA, dim=c(nvars, nvars, horizon))
#supply shocks
signs_ru[1,1,1:2] <- 1 #self
signs_ru[2,1,1:2] <- -1 #inflation
signs_ru[5,1,1] <- 0 #oil
#demand shocks
signs_ru[2,2,1:2] <- 1 #self
signs_ru[1,2,1:2] <- 1 #gdp
signs_ru[3,2,1:2] <- 1 #intrate
signs_ru[4,2,1:2] <- -1 #exrate (versatile)
signs_ru[5,2,1] <- 0 #oil
#monetary shocks
signs_ru[3,3,1:2] <- 1 #self
signs_ru[1,3,1:2] <- -1 #gdp
signs_ru[2,3,1:2] <- -1 #inflation
signs_ru[4,3,1:2] <- -1 #exrate
signs_ru[5,3,1] <- 0 #oil
#exrate shocks
signs_ru[4,4,1:2] <- 1 #self
signs_ru[2,4,1:2] <- 1 #inflation
signs_ru[3,4,1:2] <- 1 #monetary
signs_ru[5,4,1] <- 0 #oil
#oil shocks
signs_ru[5,5,1:2] <- 1 #self
signs_ru[4,5,1:2] <- -1 #exrate

zerosign_restr_bvar <- zerosign_restr.bvar(bvar_model,
                                           restr_matrix = signs_ru,
                                           LR = FALSE, tries = 1000)
irf_zerosign_bvar <- irf.ZerosignR.result(zerosign_restr_bvar, horizon = 10)
plot.ZerosignR.irf(irf_zerosign_bvar, plot_ci = TRUE)
