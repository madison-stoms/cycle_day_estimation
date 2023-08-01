

#' ARGUMENTS
#' `hnames` - vector of hormone names to be used in the analysis
#' `simlist` - list from simulate_data
#' `pve` - cut off for percent variation explained
#' `n` - number of simulation iterations


# function to find mse curve for one hormone
run_simulation = function(hnames, simlist, pve, numiter) {
  
  # reorder hnames
  rename = sort(as.numeric(str_remove(hnames, "hormone")))
  hnames = str_c("hormone", rename)
  hind = as.numeric(str_remove(hnames, "hormone"))
  
  # extract time stuff
  p = length(simlist$argvals); L = length(hnames); n = nrow(simlist$scores)
  
  # begin iteration
  sim_res = list(); iter = 1
  while(iter <= numiter) {
    
    # simulate scores
    scores = mvrnorm(1, rep(0, 3), Sigma = simlist$Sigma)
    
    # simulate new subject and construct fpca_list
    simsubj = matrix(NA, nrow = p, ncol = L); sim_fpca_list = list()
    for(m in 1:L) {
      
      l = hind[m]
      # simulate subject
      simsubj[,m] = abs(simlist$mu[,l] + simlist$eigenfunctions[,l]*scores[l] + rnorm(p, 0, sqrt(simlist$sigma2)))
      
      # fpca list
      sim_fpca_list[[m]] = list(argvals = simlist$argvals, mu = simlist$mu[,l],
                                eigenfunctions = matrix(simlist$eigenfunctions[,l], ncol = 1), scores = matrix(simlist$scores[,l], ncol = 1), 
                                sigma2 = simlist$sigma2, eigenvalues = simlist$eigenvalues[l], npc = 1)
      
    }
    names(sim_fpca_list) = hnames
    
    # extract ynew
    ind = sample(1:nrow(simsubj), 1)
    ynew = simsubj[ind,]; theta = simlist$argvals[ind]
    
    # calculate loglik
    ll = loglik(ynew, hnames, sim_fpca_list, pve)
    mle = ll$thetahat
    
    # calculate error
    
    daysoff = 28*sign(theta - mle) * min(abs(theta - mle), 1 - abs(theta - mle))  
    
    
    sim_res[[iter]] = tibble(iter = iter, hnames_set = paste(hnames, collapse = ", "), pve = pve,  
                     theta = theta, mle = ll$thetahat, daysoff = daysoff)
    
    iter = iter + 1
    
  }
  
  return(do.call(rbind, sim_res))
  
}
