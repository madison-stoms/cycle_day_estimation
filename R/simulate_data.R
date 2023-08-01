

#' @param n number of subjects in simulated data
#' @param sigma2 observation error variance
#' @param argvals grid on which estimated components will be evaluated

simulate_data = function(n, sigma2, argvals) {
  
  # create basis
  spline = bs(argvals, knots = 0:6/6)
  
  # create hormone 1 pc
  coef1 = c(0,0,0,0,0,1,0,0,0,0)
  pc1 = spline %*% coef1
  mu1 = pc1 + 0.08
  
  # create hormone 2 pc
  coef2 = c(0,0,0,0,1,0,0,0,0,0)
  pc2 = spline %*% coef2
  mu2 = pc2 + 0.08
  
  # create hormone 3 pc
  pc3 = matrix(rep(1, length(argvals)), ncol = 1)
  mu3 = pc3 - 0.4
  
  # create scores with correlation
  Sigma = matrix(c(0.1, 0.05, 0.05, 0.05, 0.1, 0.05, 0.05, 0.05, 0.04), nrow = 3)
  scores = mvrnorm(n, c(0,0,0), Sigma)
  
  # start data generating loop
  simdat = tibble()
  for(i in 1:n) {
  
  simdat = simdat %>%
    rbind(tibble(subj = i, argvals = argvals,
                 hormone1 = c(mu1 + pc1*scores[i,1] + rnorm(length(argvals), 0, sqrt(sigma2))),
                 hormone2 = c(mu2 + pc2*scores[i,2] + rnorm(length(argvals), 0, sqrt(sigma2))),
                 hormone3 = c(mu3 + pc3*scores[i,3] + rnorm(length(argvals), 0, sqrt(sigma2)))))
  
  }
  
  simdat = simdat %>%
    pivot_longer(cols = starts_with("h"), names_to = "hormone", values_to = "y") %>%
    mutate(hormone = factor(hormone, levels = c("hormone1", "hormone2", "hormone3"))) %>%
    mutate(y = ifelse(y < 0, 0, y))
  
  return(list(simdat = simdat, argvals = argvals, mu = cbind(mu1, mu2, mu3), scores = scores,
              eigenfunctions = cbind(pc1, pc2, pc3), 
              scores = scores, Sigma = Sigma, sigma2 = sigma2, eigenvalues = c(0.1, 0.1, 0.2)))

}
