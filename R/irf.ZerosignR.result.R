irf.ZerosignR.result <- function(zerosignr, horizon=0, LR=FALSE, quantil=0.05){
  #'Calculate structural IRFs for zero and sign restrictions estimation result
  #'
  #'@description
  #'This function calculates arbitrary-horizon structural IRFs conditinal on given
  #'set of matrices Q --- orthonormal matrices that satisfy zero and sign restrictions.
  #'
  #'@param zerosignr \code{ZerosignR.result} object
  #'@param horizon Numeric. Horizon to calculate IRF for.
  #'0 = only contemporaneous IRF, 1 = contemporaneous + next period after shock IRF, etc.
  #'@param LR Boolean. Whether to calculate long-run IRF additionally to finite-run ones.
  #'@param quantil Numeric. Number in (0,1), calculate IRF quantile band.
  #'@return
  #'\code{ZerosignR.irf}, a list containing:
  #'\itemize{
  #'\item \code{struc_irfs}, 3-dim. array with IRFs. The 1-st dimension is for specific Q-matrix,
  #'the 2-nd one is \emph{nvars*(horizon+contemporaneous+LR)} --- row-stacked IRFs for different
  #'periods, starting from the earliest one, the 3-rd one is for shocks;
  #'\item \code{median_irf}, pointwise-median IRF, aggregated for Q-matrices;
  #'\item \code{quantile_irf}, IRF sample quantiles for band;
  #'\item \code{ctm_irf}, IRF from closest-to-median model, Euclidean distance;
  #'\item \eqn{\ldots}{...}, auxillary constants.
  #'}
  #'
  #'@note
  #'\code{struc_irfs} in returned value \code{ZerosignR.irf} is a row-stacked IRF with columns for shocks.
  #'Shocks' order is as in the restriction matrix from \code{ZerosignR.result} object.
  #'
  #'First \emph{nvars} rows are for variables' contemporaneous response (order is as in underlying model),
  #'Rows from \emph{nvars + 1} to \emph{2*nvars} are for variables' response in the period right after the shock, etc.
  #'If \code{LR} is \code{TRUE}, the last \emph{nvars} rows are long-run IRF.
  #'
  #'@details
  #'Arias et al. (2018) show that structural shock IRF for period \eqn{h}{h} can be calculated in
  #'matrix form as follows (see definitions for \eqn{J}{J} and \eqn{F}{F} in the original article):
  #'\deqn{L_h(A_0, A_+) = \left( A_0^{-1} J' F^h J \right)'.}{Lh(A0, A+) = (A0 * J' * F^h * J)'.}
  #'
  #'The special case is for long-run IRF, which is given as:
  #'\deqn{L_\infty(A_0, A_+) = \left( A^'_0 - \sum_{l=1}^{p} A^'_l \right)^{-1}.}{Linf(A0, A+) = (A0' - ∑ A_l')^(-1).}
  #'
  #'In order to obtain structural IRF conditional on orthonormal matrix Q, that satisfies
  #'zero and sign restrictions, one has just multiply IRF matrix by Q:
  #'\deqn{L_h(A_0, A_+)Q = L_h(A_0Q, A_+Q).}{Lh(A0, A+) * Q = Lh(A0 * Q, A+ * Q).}
  #'
  #'@seealso \code{\link{zerosign_restr}}
  #'@author Artur Zmanovskii. \email{anzmanovskii@gmail.com}
  #'@references{
  #'Arias, J.E. and Rubio-Ramirez, J. F. and Waggoner, D. F. (2018)
  #'Inference Based on Structural Vector Autoregressions Identifiied with Sign and Zero Restrictions:
  #'Theory and Applications. \emph{Econometrica}, 86, 2, 685-720, \url{https://doi.org/10.3982/ECTA14468}.
  #'}

  n = zerosignr$meta$nvars
  p = zerosignr$meta$lags
  accepted = zerosignr$meta$accepted
  varnames <- zerosignr$meta$varnames
  horizon_matrix <- ifelse(LR, horizon + 2, horizon + 1)
  struc_irfs = array(dim=c(accepted, horizon_matrix*n, n))
  for(i in 1:accepted){
    B_draw <- zerosignr$model[[i]]$B
    Sigma_draw <- zerosignr$model[[i]]$Sigma
    Q_est <- zerosignr$model[[i]]$Q
    irf_new <- irf_ala_arias(B = B_draw, Sigma = Sigma_draw, p = p, n = n,
                             horizon = horizon, LR = LR, Q = Q_est)
    struc_irfs[i,,] <- irf_new
  }
  #Median
  median_irfs <- apply(struc_irfs, c(2,3), median)
  #CTM
  dist_vec <- vector(length=accepted)
  for(i in 1:accepted){
    dist_vec[i] <- sum((median_irfs - struc_irfs[i,,])^2)
  }
  ctm_index <- which.min(dist_vec)
  ctm_irfs <- struc_irfs[ctm_index,,]
  #Quantiles
  quantile_irfs <- array(dim=c(2,horizon_matrix*n, n))
  quantile_irfs[1,,] <- apply(struc_irfs, c(2,3), quantile, probs=(quantil/2))
  quantile_irfs[2,,] <- apply(struc_irfs, c(2,3), quantile, probs=(1-quantil/2))
  transformed_irfs <- list("struc_irfs" = struc_irfs, "median_irf" = median_irfs,
                           "ctm_irf" = ctm_irfs, "quantile_irf" = quantile_irfs)
  transformed_irfs$accepted <- accepted
  transformed_irfs$n <- n
  transformed_irfs$p <- p
  transformed_irfs$modelclass <- zerosignr$meta$modelclass
  transformed_irfs$varnames <- varnames
  transformed_irfs$restr_matrix <- zerosignr$meta$restr_matrix
  transformed_irfs$quantile <- quantil
  class(transformed_irfs) <- "ZerosignR.irf"
  return(transformed_irfs)
  }

