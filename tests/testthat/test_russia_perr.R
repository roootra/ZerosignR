#'Russia PERR estimation
#'
require(BVAR)
require(vars)
require(testthat)
set.seed(42)
data_to_model <- read.csv("russian_data.csv")
data_to_model <- data_to_model[,-1]
data_to_model$neer_qoq <- data_to_model$neer_qoq*-1

bvar_model <- bvar(data_to_model, lags=1,
                   priors = bv_priors(hyper="full"),
                   mh = bv_mh(),
                   n_thin=2, n_burn=100, n_draw=200)
var_model <- VAR(data_to_model, p = 1, type = "const")

#Zero and sign matrix
horizon <- 2
nvars <- NCOL(data_to_model)
signs <- array(NA, dim=c(nvars, nvars, horizon))
#supply shocks
signs[1,1,1:2] <- 1 #self
signs[5,1,1] <- 0 #oil
#demand shocks
signs[2,2,1:2] <- 1 #self
signs[5,2,1] <- 0 #oil
#monetary shocks
signs[3,3,1:2] <- 1 #self
signs[5,3,1] <- 0 #oil
#exrate shocks
signs[4,4,1:2] <- 1 #self
signs[5,4,1] <- 0 #oil
#oil shocks
signs[5,5,1:2] <- 1 #self
signs[4,5,1:2] <- -1 #exrate

zerosign_restr_bvar <- zerosign_restr.bvar(bvar_model,
                                           restr_matrix = signs,
                                           LR = FALSE, tries = 100)
irf_zerosign_bvar <- irf.ZerosignR.result(zerosign_restr_bvar, horizon = 10)

test_that("Check whether there are accepted models for a simple setup for bvar",{
  expect_gt(length(zerosign_restr_bvar$models), 0)
})
test_that("Check whether zero restictions are fulfilled for bvar",{
  expect_equal(irf_zerosign_bvar$median_irf[5,1:3], c(0,0,0))
})
test_that("Check whether normalizing sign restictions are fulfilled for bvar",{
  expect_equal(sum(sign(diag(irf_zerosign_bvar$median_irf[1:5,1:5]))), 5)
  expect_equal(sum(sign(diag(irf_zerosign_bvar$median_irf[6:10,1:5]))), 5)
})
test_that("Check whether other sign restictions are fulfilled for bvar",{
  expect_equal(sign(irf_zerosign_bvar$median_irf[4,5]), -1)
  expect_equal(sign(irf_zerosign_bvar$median_irf[9,5]), -1)
})

zerosign_restr_varest <- zerosign_restr.varest(var_model,
                                               restr_matrix = signs,
                                               LR = FALSE, tries = 1000)
irf_zerosign_varest <- irf.ZerosignR.result(zerosign_restr_varest, horizon = 10)
test_that("Check whether there are accepted models for a simple setup for varest",{
  expect_gt(length(zerosign_restr_varest$models), 0)
})
test_that("Check whether zero restictions are fulfilled for varest",{
  expect_equal(irf_zerosign_varest$median_irf[5,1:3], c(0,0,0))
})
test_that("Check whether normalizing sign restictions are fulfilled for varest",{
  expect_equal(sum(sign(diag(irf_zerosign_varest$median_irf[1:5,1:5]))), 5)
  expect_equal(sum(sign(diag(irf_zerosign_varest$median_irf[6:10,1:5]))), 5)
})
test_that("Check whether other sign restictions are fulfilled for varest",{
  expect_equal(sign(irf_zerosign_varest$median_irf[4,5]), -1)
  expect_equal(sign(irf_zerosign_varest$median_irf[9,5]), -1)
})

