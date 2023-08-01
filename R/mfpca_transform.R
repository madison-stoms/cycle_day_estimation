

#' @param hnames vector of hormone names to be used in the analysis
#' @param fpca_list list of univariate FPCA estimates
#' @param pve cut off for percent variation explained


mfpca_transform = function(hnames, fpca_list, pve = 0.95) {
  
  # step 1: Extact Univariate FPCA Params
  Phi = list(); Sigtemp = list()
  Mlist = list(); M = c()
  for(l in 1:length(hnames)) {
    
    index = which(names(fpca_list) == hnames[l])
    res = fpca_list[[index]]
    K = res$npc
    
    Phi[[l]] = res$eigenfunctions
    Sigtemp[[l]] = res$scores
    
    if(l == 1) { Mlist[[l]] = 1:K } 
    if(l != 1) { Mlist[[l]] = max(Mlist[[l-1]]) + 1:K }
    
    M = c(M, K)
    
  }
  Sig = do.call(cbind, Sigtemp)
  n = nrow(Sigtemp[[1]]); p = nrow(Phi[[1]])
  Mplus = sum(M)
  
  
  # Step 2: Construct Z
  Z = 1/(nrow(Sig)-1) * t(Sig) %*% Sig
  
  
  # Step 3: Perform Eigenanalysis on Z
  Zdecomp = eigen(Z)
  cm = Zdecomp$vectors
  vm = Zdecomp$values
  
  pvevec = vm / sum(vm)
  npc = min(which(cumsum(pvevec) > pve))
  
  
  # Step 4: Update Eigenfunctions and Scores
  psi = list(); rho = list()
  for(l in 1:length(hnames)) {
    psitemp = matrix(NA, nrow = p, ncol = npc)
    rho = matrix(NA, nrow = n, ncol = npc)
    for(m in 1:npc) {
      
      psitemp[,m] = Phi[[l]] %*% cm[Mlist[[l]],m] # eigenfns
      for(i in 1:n) { rho[i,m] = Sig[i,] %*% cm[,m] } # scores
      
    }
    psi[[l]] = psitemp
  }
  
  # name components
  names(psi) = hnames
  
  
  return(list(psi = psi, rho = rho, argvals = res$argvals, 
              vm = vm[1:npc], npc = npc, cm = cm))
  
}

