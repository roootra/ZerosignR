Q_mod <- function(Q){
  for(i in 1:length(diag(Q))){
    if(diag(Q)[i] < 0){
      Q[,i] <- -1*Q[,i]
    }
  }
  Q
}
