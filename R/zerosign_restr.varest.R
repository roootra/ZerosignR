zerosign_restr.varest <- function(varest_model, restr_matrix, LR=FALSE, tries=300){
  require(vars)
  if(!(varest_model$type %in% c("none", "const"))){
    stop(simpleError("VARs with trend are not supported."))
  }
  B = t(vars::Bcoef(varest_model)[,c(NCOL(vars::Bcoef(varest_model)),
                                     (1:(NCOL(vars::Bcoef(varest_model)) - 1)))])
  Sigma = summary(varest_model)$covres #how do I reach summary from vars namespace???
  p = varest_model$p
  n = varest_model$K
  has_const = (varest_model$type == "const")
  varnames <- colnames(varest_model$y)
  restr_matrix_stacked = NULL
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
  bayesian_ortho_irfs = irf_ala_arias(B=B, Sigma=Sigma,
                                      p=p, n=n, horizon=horiz, LR=LR)
  succ_models_from_draw = zerosign_restr_ala_arias(irfs=bayesian_ortho_irfs, B=B, Sigma=Sigma,
                                                   zero_sign_matrix = restr_matrix_stacked,
                                                   tries = tries)
  if(length(succ_models_from_draw) != 0){
    satisfying_models = append(satisfying_models, succ_models_from_draw)
  }
  cat("\rFrequentist model. Accepted models: ", length(succ_models_from_draw), ".", sep="")
  cat("\n", sep="")
  if(length(satisfying_models) == 0){
    cat("No staisfying models are found!")
  }
  toreturn <- list("models" = satisfying_models)
  class(toreturn) <- "ZerosignR.result"
  toreturn$meta$varnames <- varnames
  toreturn$meta$modelclass <- "varest"
  toreturn$meta$restr_matrix <- restr_matrix
  toreturn$meta$lags <- varest_model$p
  toreturn$meta$nvars <- varest_model$K
  toreturn$meta$LR <- LR
  toreturn$meta$tries <- tries
  toreturn$meta$models_checked <- tries
  toreturn$meta$accepted <- length(satisfying_models)
  return(toreturn)
}
