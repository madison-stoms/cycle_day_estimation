

#' @param ynew vector of cross-sectional hormone values
#' @param thetatilde candidate theta
#' @param hnames vector of hormone names to be used in the analysis
#' @param fpca_list list of univariate FPCA esimates
#' @param pve percent variance explained by the included components


loglik_thetatilde = function(ynew, thetatilde, hnames, fpca_list, pve) {
  
  if(!(thetatilde %in% fpca_list[[1]]$argvals)) { stop("FPCA parameters not evaluated at thetatilde!") }
  
  # univariate case
  if(length(hnames) == 1) {
    
    # get fpca results
    res = fpca_list[[which(names(fpca_list) == hnames)]]
    if(is.null(res)) { stop(str_c("Hormone ", hnames, " is not included in fpca_list!")) }
    ind = which(thetatilde == res$argvals)
    
    # pull theta points
    phicur = matrix(res$eigenfunctions[ind, ], nrow = 1)
    mucur = res$mu[ind]
    
    # account for npc = 1
    if(length(res$eigenvalues) == 1) { emat = 1/res$eigenvalues }
    else { emat = diag(1/res$eigenvalues) }
    
    # calculate hypothetical score
    hscores = solve( t(phicur) %*% phicur + res$sigma2 * emat ) %*% t(phicur) %*% (ynew - mucur)
    
    # calculate normal densities
    p1 = dnorm(ynew, mean = mucur + phicur %*% hscores, sd = sqrt(res$sigma2))
    if(res$npc == 1) { p2 = as.numeric(dnorm(hscores, mean = 0, sd = sqrt(res$eigenvalues))) }
    if(res$npc != 1) { p2 = mvtnorm::dmvnorm(t(hscores), mean = rep(0,res$npc), sigma = diag(res$eigenvalues)) }
    
    ll = log(p1) + log(p2)
    
  }
  
  
  # multivariate case
  if(length(hnames) != 1) {
    
    mfpca_list = mfpca_transform(hnames, fpca_list, pve)
      
    Sig = c()
    for(l in 1:length(ynew)) {
      
      # get fpca results
      res = fpca_list[[which(names(fpca_list) == hnames[l])]]
      if(is.null(res)) { stop(str_c("Hormone ", hnames[l], " is not included in fpca_list!")) }
      ind = which(thetatilde == res$argvals)
      
      # pull theta points
      phicur = matrix(res$eigenfunctions[ind, ], nrow = 1)
      mucur = res$mu[ind]
      
      # account for PCnum = 1
      if(length(res$eigenvalues) == 1) { emat = 1/res$eigenvalues }
      else { emat = diag(1/res$eigenvalues) }
      
      # calculate hypothetical scores
      hscores = solve( t(phicur) %*% phicur + res$sigma2 * emat ) %*% t(phicur) %*% (ynew[l] - mucur)
      
      # concatenate scores
      Sig = c(Sig, hscores)
      
    }
    
    # transform scores
    npc = mfpca_list$npc; rho = c()
    for(m in 1:npc) { rho[m] = Sig %*% matrix(mfpca_list$cm[,m], ncol = 1) }
    
    # start likelihood loop
    ll = 0
    for(l in 1:length(ynew)) {
      
      # get fpca results
      res = fpca_list[[which(names(fpca_list) == hnames[l])]]
      ind = which(thetatilde == res$argvals)
      
      # get components
      psicur = matrix(mfpca_list$psi[[l]][ind, ], nrow = 1)
      mucur = res$mu[ind]
      
      # calculate likelihood probabilities
      p1 = dnorm(ynew[l], mean = mucur + psicur %*% rho, sd = sqrt(res$sigma2))
      
      if(npc == 1) { p2 = dnorm(rho, mean = 0, sd = sqrt(mfpca_list$vm)) }
      if(npc != 1) { p2 = mvtnorm::dmvnorm(t(rho), mean = rep(0,npc), sigma = diag(mfpca_list$vm)) }
      
      ll = log(p1) + log(p2) + ll
      
    }
  }
  
  return(tibble(thetatilde = thetatilde, ll = ll))
  
}
