

#' @param ynew vector of cross-sectional hormone values
#' @param hnames vector of hormone names to be used in the analysis
#' @param fpca_list list of univariate fpca esimates
#' @param pve cut off for percent variation explained

loglik = function(ynew, hnames, fpca_list, pve = 0.95) {
  
  
  if(length(hnames) != 1) { mfpca_list = mfpca_transform(hnames, fpca_list, pve)}
  if(length(hnames) == 1) { mfpca_list = NULL }
  
  full_ll = fpca_list[[1]]$argvals %>%
    map_dfr( ~ loglik_thetatilde(ynew, thetatilde = .x, hnames, fpca_list, pve))
  
  mle = full_ll %>% 
    filter(ll == max(ll)) %>%
    pull(thetatilde)
  
  if(!is.null(length(mle))) { mle = sample(mle,1) }
  
  return(list(thetahat = mle, ll = full_ll))
  
}


