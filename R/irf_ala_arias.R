irf_ala_arias <- function(B, Sigma, p, n, horizon, LR=FALSE, Q=NULL){
  #' Calculate structural IRFs
  #'
  #' This function calculates orthogonal structural IRFs
  #' using linear algebra primitives. Returns a \eqn{n \times (n \cdot horizon + n \cdot 1_{LR})}{n x (n*horizon + 1*LR)} matrix,
  #' columns of which represent shocks, as ordered in model, and rows represent variables for each period of
  #' IRF, as ordered in model, written cyclically. First n rows represent contemporaneous relations of
  #' variables, next n rows represent same relations for h=1, and so on.
  #'
  #' @param B \eqn{B = A_+ \cdot (A_0)^(-1)}{B = A+ * A0^(-1)} --- matrix of reduced parameters in form \eqn{B = \left[c, B_1, \ldots, B_p \right]}{B = (c, B1, ..., Bp)}.
  #' If there are no constants in model, just add extra zero column to the beginning of B.
  #' @param Sigma Variance-covariance matrix of error term.
  #' @param p Order of model.
  #' @param n Number of variables/shocks, including unidentified ones.
  #' @param horizon Number of periods to calculate IRF for. -1 => no short-run IRFs.
  #' @param LR Boolean variable, whether to calculate long-run IRF (default=FALSE). Note if LR is TRUE, long-run IRFs.
  #' are added to the bottom of matrix returned (the last n rows are long-run IRFs).
  #' @param Q \emph{(optional)} Estimated orthonormal matrix used for structural IRF (default=NULL).
  #' @return A matrix of columnwise-stacked structural IRFs, columns = shocks (as ordered in model), rows = variables*horizon.
  require(expm)

  A0_inv = t(chol(Sigma)) #upper or lower?
  A0 = solve(A0_inv)
  A_plus = B %*% A0
  if(horizon >=0){
    #Finite-run IRFs
    F_mat = matrix(data=0, nrow = p*n, ncol = p*n)
    F_mat[,1:n] = t(B[-1,])
    if(p > 1){
      F_mat[1:((p-1)*n), (n+1):NCOL(F_mat)] = diag(nrow=(p-1)*n, ncol = (p-1)*n)
    }

    J_mat = matrix(data=0, nrow=p*n, ncol=n)
    J_mat[1:n, 1:n] = diag(nrow=n, ncol=n)

    irfs = matrix(nrow=(horizon+1)*n, ncol=n)
    for(time in 0:horizon){
      irfs[(1 + n*(time)):(n*(time+1)),] =
        A0_inv %*% t(J_mat) %*% (F_mat %^%(time)) %*% J_mat
      #delete A0_inv to get usual Phi
      #otherwise, get already Cholesky ortho IRFs
    }
  }
  else{
    irfs = matrix(nrow=0, ncol=n)
  }
  if(LR){
    #Long-run IRFs
    A_sum = 0
    for(lag in 1:p){
      A_sum = A_sum + t(A_plus)[1:n, ((2 + n*(lag-1)):(1 + n*lag))]
    }
    irfs_lr <- solve(t(A0) - A_sum)
    irfs <- rbind(irfs, irfs_lr)
  }
  if(!is.null(Q)){
    stopifnot("Orthonormal matrix Q should be square!"= NROW(Q) == NCOL(Q))
    stopifnot("Number of variables and dimension of matrix Q are incompatible."= NCOL(irfs)==NCOL(Q))
    irfs = irfs %*% Q
  }
  return(irfs)
}
