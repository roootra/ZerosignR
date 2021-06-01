zerosign_restr.bvar <- function(bvar_model, restr_matrix, LR=FALSE, tries=300){
  #B, Sigma, p, n, draws, sr_sign_matrix, has_const=TRUE,
  B = bvar_model$beta
  Sigma = bvar_model$sigma
  p = bvar_model$meta$lags
  n = bvar_model$meta$M
  draws = bvar_model$meta$n_save
  has_const = bvar_model$meta$K - bvar_model$meta$M
  restr_matrix_stacked = NULL
  varnames <- bvar_model$variables
  SR = TRUE
  for(period in 1:dim(restr_matrix)[3]){
    restr_matrix_stacked = rbind(restr_matrix_stacked, restr_matrix[,,period])
  }
  #calculate horizon
  if(LR){
    horiz = dim(restr_matrix)[3]-2
  }
  else{
    horiz = dim(restr_matrix)[3]-1
  }
  cat("\n")
  cat("Restrictions horizon: ",
      ifelse(LR, paste0(ifelse(horiz < 0, "Only", paste0(horiz, " +")), " long-run"), horiz), "\n", sep="")
  satisfying_models = list()
  for(draw in 1:draws){
    cat("\rDraw ", draw, " of ", draws,". Accepted draws: ", length(satisfying_models), ".", sep="")
    bayesian_ortho_irfs = irf_ala_arias(B=B[draw,,], Sigma=Sigma[draw,,],
                                        p=p, n=n, horizon=horiz, LR=LR)
    succ_models_from_draw = zerosign_restr_ala_arias(irfs=bayesian_ortho_irfs, B=B[draw,,], Sigma=Sigma[draw,,],
                                                     zero_sign_matrix = restr_matrix_stacked,
                                                     tries = tries)
    if(length(succ_models_from_draw) != 0){
      satisfying_models = append(satisfying_models, succ_models_from_draw)
    }
  }
  cat("\n", sep="")
  if(length(satisfying_models) == 0){
    cat("No staisfying models are found!")
  }
  toreturn <- list("models" = satisfying_models)
  class(toreturn) <- "ZerosignR.result"
  toreturn$meta$varnames <- varnames
  toreturn$meta$modelclass <- "bvar"
  toreturn$meta$lags <- bvar_model$meta$lags
  toreturn$meta$nvars <- bvar_model$meta$M
  toreturn$meta$restr_matrix <- restr_matrix
  toreturn$meta$LR <- LR
  toreturn$meta$tries <- tries
  toreturn$meta$models_checked <- tries*draws
  toreturn$meta$accepted <- length(satisfying_models)
  return(toreturn)
}
