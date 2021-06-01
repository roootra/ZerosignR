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
signs <- array(NA, dim=c(nvars, nvars, horizon))
#supply shocks
signs[1,1,1:2] <- 1 #self
signs[2,1,1:2] <- -1 #inflation
signs[5,1,1] <- 0 #oil
#demand shocks
signs[2,2,1:2] <- 1 #self
signs[1,2,1:2] <- 1 #gdp
signs[3,2,1:2] <- 1 #intrate
#signs[4,2,1:2] <- -1 #exrate (versatile)
signs[5,2,1] <- 0 #oil
#monetary shocks
signs[3,3,1:2] <- 1 #self
signs[1,3,1:2] <- -1 #gdp
signs[2,3,1:2] <- -1 #inflation
signs[4,3,1:2] <- -1 #exrate
signs[5,3,1] <- 0 #oil
#exrate shocks
signs[4,4,1:2] <- 1 #self
signs[2,4,1:2] <- 1 #inflation
signs[3,4,1:2] <- 1 #monetary
signs[5,4,1] <- 0 #oil
#oil shocks
signs[5,5,1:2] <- 1 #self
signs[4,5,1:2] <- -1 #exrate

zerosign_restr_bvar <- zerosign_restr.bvar(bvar_model,
                                           restr_matrix = signs,
                                           LR = FALSE, tries = 1000)
irf_zerosign_bvar <- irf.ZerosignR.result(zerosign_restr_bvar, horizon = 10)
