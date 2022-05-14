plot.ZerosignR.irf <- function(zerosign_irf, plot_ci=TRUE, plot_ctm=FALSE){
  #'Plot structural IRFs for zero and sign restriction algorithm.
  #'
  #'@description This function plots IRFs for all the structural shocks.
  #'
  #'@param zerosign_irf \code{Zerosign.irf}, an object containing structural IRFs
  #'obtained using zero and sign algorithm.
  #'@param plot_ci Logical. Whether to plot confidence intervals.
  #'@param plot_ctm Logical. Whether to plot IRF from closest-to-median model instead of
  #'median IRF.
  #'
  #'@return \code{NULL}.
  #'
  #'@seealso \code{\link{zerosign_restr}, \link{irf.ZerosignR.result}}
  #'@author Artur Zmanovskii. \email{anzmanovskii@gmail.com}
  var_number <- zerosign_irf$n
  periods <- NROW(zerosign_irf$median_irf) %/%  var_number
  par(mfcol = c(var_number, var_number), mar = c(2.5, 2.25, 2, 0.5)) #mai = c(0.3, 0.4, 0.3, 0.2)
  if(plot_ctm){
    main_irf <- zerosign_irf$ctm_irf
  }
  else{
    main_irf <- zerosign_irf$median_irf
  }
  for(i_shock in 1:var_number){
    for(i_variable in 1:var_number){
      current_series <- rep(0.0, periods)
      irf_path_indices <- i_variable + var_number*(0:(periods-1))
      titl <- sprintf('%s → %s', zerosign_irf$shock_names[i_shock],
                                            zerosign_irf$var_names[i_variable])
      plot(main_irf[irf_path_indices, i_shock], main = titl, col = 'red', type = 'l',
           cex.main=1, xlab='Periods', ylab=NULL)
      if(plot_ci){
        polygon(c(1:periods, rev(1:periods)),
                  c(zerosign_irf$quantile_irf[1, irf_path_indices, i_shock],
                    rev(zerosign_irf$quantile_irf[2, irf_path_indices, i_shock])),
                      col = "grey90", border = NA)
        lines(main_irf[irf_path_indices, i_shock], col = 'red') # again for IRF to be visible
        lines(1:periods, zerosign_irf$quantile_irf[1, irf_path_indices, i_shock],
                                                      lty = 'dashed', col = 'blue')
        lines(1:periods, zerosign_irf$quantile_irf[2, irf_path_indices, i_shock],
                                                      lty = 'dashed', col = 'blue')
      }
      abline(h = 0, col = 'grey', lty = 'dashed')


    }
  }
}
