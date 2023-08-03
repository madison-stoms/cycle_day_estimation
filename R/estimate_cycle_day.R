
#' @param ynew vector of matrix (columns as hormones) of cross-sectional hormone values
#' @param hnames vector of hormone names to be used in the analysis
#' @param fpca_list list of univariate fpca esimates
#' @param pve cut off for percent variation explained

estimate_cycle_day = function(ynew, hnames, fpca_list, pve) {
  
  if(is.null(nrow(ynew))) { res = tibble(yin = ynew); names(res) = hnames }
  else{ res = as_tibble(ynew); names(res) = hnames }
  
  mle = c()
  for(i in 1:nrow(ynew)) { mle[i] = loglik(ynew[i,], hnames, fpca_list, pve)$thetahat }
  
  return(mutate(res, cycle_day = mle))
  
}
