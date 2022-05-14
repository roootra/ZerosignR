zerosign_restr_ala_arias <- function(irfs, B, Sigma, zero_sign_matrix, tries,
                                     negative_Q=TRUE, perm_Q=FALSE){
  #' Perform simulations to identify shock under zero and sign restrictions — primitives
  #'
  #' @details
  #' This function simulates orthognormal random matrices \emph{Q} given number of times for
  #' given IRF matrix in order to identify structural shocks. It returns a list of
  #' \emph{Q} matrices and transformed IRFs that satisfy imposed restrictions.
  #'
  #' @param irfs IRF matrix of size \eqn{(nvars \times (nvars*horizons))}{nvars x (nvars*horizons)} ---
  #' columnwise-stacked IRFs by periods if there is more than one period of IRF.
  #' @param B \eqn{B = A_+ \cdot (A_0)^(-1)}{B = A+ * A0^(-1)} --- matrix of reduced parameters in form \eqn{B = \left[c, B_1, \ldots, B_p\right]}{B = (c, B1, ..., Bp)}.
  #' If there are no constants in model, just add extra zero column to the beginning of B.
  #' @param Sigma Variance-covariance matrix of error term.
  #' @param zero_sign_matrix Square matrix of zero and sign restrictions
  #' \eqn{(nvars \times nvars)}{nvars x nvars}.
  #' Columns = shocks, structural and residual ones. Put \emph{0} into \emph{(i,j)}-th cell of
  #' this matrix in order to set zero restriction to the \emph{i}-th variable response to
  #' the shock \emph{j}.
  #' @param tries Number of tries of random orthonormal matrix generation
  #' @param negative_Q (default = \code{TRUE}). Whether also check negative Q.
  #' @param perm_Q NOT IMPLEMENTED. Whether to account for permutation of matrix \emph{Q}
  #' columns to increase method's efficiency.
  #'
  #' @details
  #' Note that for matrix of zero and sign restrictions matrix
  #' user should order ariables in descending order of number of zero restrictions. If there are
  #' two or more variables with the same number of zero restrictions, their order
  #' can be arbitrary. There cannot be more than \emph{nvars - j} zero restrictions
  #' for \emph{j} th column (in total for all time periods).
  #'
  if(perm_Q){
    cat("Permutations of Q matrix are not implemented. Proceeding without them...\n")
  }
  nvars = NCOL(irfs)
  periods = NROW(irfs) / nvars
  S = na.fill(zero_sign_matrix, 0)
  Z = (zero_sign_matrix == 0)
  Z = na.fill(Z, 0)
  zero_restrictions_present = ifelse(any(Z == 1), TRUE, FALSE)
  satisfying_models = list()
  succ = 0
  fails = 0
  for(try in 1:tries){
    #each irf in period h (column-stacked) should be transformed
    if(zero_restrictions_present){
      #If there are sign restrictions, impose linear constraints
      #on Q and calculate it recursively, as given by Theorem 4 in (Arias et al., 2014)
      #Recursively-constructed Q-matrix is still orthonormal!
      Q = matrix(nrow=nvars, ncol=nvars)
      for(j in 1:nvars){
        Z_j = diag(mapply(Z[,j], FUN=as.numeric))
        Z_j = t(Z_j[,(colSums(Z_j) != 0)])
        n_zeros = sum(Z_j)
        if(j > 1){
          if(n_zeros > 0){
            R_j = rbind(Z_j %*% irfs, t(Q[,1:(j-1)]))
          }
          else{
            R_j = t(Q[,1:(j-1)])
          }
        }
        else{
          stopifnot("Please check that shocks are ordered in descending order of number of zero restrictions."= NCOL(Z_j) != 0)
          R_j = Z_j %*% irfs
        }
        N_jminus = Null(t(R_j))
        stopifnot("Error: check the number of restrictions for each shock: there cannot be more restrictions for one shock (in total for all IRF periods) than total number of shocks minus order of this shock."= NCOL(N_jminus) != 0)
        n_j = NCOL(N_jminus)
        y_j = rnorm(n_j)
        Q[,j] = N_jminus %*% y_j / norm(y_j, type="2")
      }
    }
    else{
      #If there are no zero restrictions, use Theorem 1
      X = matrix(rnorm(nvars^2), nrow=nvars)
      Q = qr.Q(qr(X))
    }
    #check for signs
    Q_minus = -1*Q
    flag_fail = FALSE
    #TBD: check and *-1 columns instead of whole matrix
    for(j in 1:nvars){
        S_j = diag(S[,j])
        S_j = t(S_j[,(colSums(S_j) != 0)])
        e_j = diag(nrow=nvars)[,j]
        if(any(S_j %*% irfs %*% Q[,j] < 0)){
          flag_fail = TRUE
          break
        }
    }
    if(flag_fail & negative_Q){ #Also check negative Q
      for(j in 1:nvars){
        S_j = diag(S[,j])
        S_j = t(S_j[,(colSums(S_j) != 0)])
        e_j = diag(nrow=nvars)[,j]
        if(any(S_j %*% irfs %*% Q_minus[,j] < 0)){
          flag_fail = TRUE
          break
        } else{flag_fail = FALSE}
      }
      if(!flag_fail){Q <- Q_minus}
    }
    if(flag_fail){
      fails = fails + 1
    } else{
      irfs_transformed = irfs %*% Q
      succ = succ + 1
      satisfying_models[[length(satisfying_models) + 1]] =
        list("Q" = Q, "B" = B, "Sigma" = Sigma, "Transformed_irfs" = irfs_transformed)
    }
  }
  return(satisfying_models)
}
