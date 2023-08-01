

#' @param data data frame with four arguments: (1) subj: subject ID; (2) hormone: name of hormone; (3) argvals: observation times; (4) y: hormone levels
#' @param hnames vector of hormone names to be used in the analysis
#' @param pve cut off for percent variation explained
#' @param ntheta number of time points to test for subject
#' @param argvals new grid to run fpca on


loocv_analysis = function(data, hnames, pve, ntheta, argvals = NULL) {
  
  suppressWarnings({
  
  # prep data
  data = data %>%
    filter(hormone %in% hnames) %>%
    mutate(hormone = factor(hormone, levels = hnames))
  
  IDs = sort(unique(data$subj))
  
  if(is.null(argvals)) { argvals = sort(unique(data$argvals)) }
  
  # start iterating across IDs
  reslist = list()
  for(i in 1:length(IDs)) {
    
    # define datasets 
    train_dat = data %>% filter(subj != IDs[i]) 
    test_dat = data %>% filter(subj == IDs[i])
    
    # run FPCA on hormone set
    fpca_list = create_fpca_list(train_dat, hnames, argvals, pve)
    
    # choose common testing points (will be averaged across)
    testtemp = test_dat %>% group_by(argvals) %>% summarise(n = n()) %>% 
      filter(n == length(hnames)) %>% pull(argvals) 
    
    # uniformly sample across grid
    if(length(testtemp) < ntheta) { ntheta = length(testtemp) }
    sub1 = testtemp[which(testtemp < 0.5)]
    sub2 = testtemp[which(testtemp >= 0.5)]
    n1 = ifelse(length(sub1) < ceiling(ntheta/2), length(sub1), ceiling(ntheta/2))
    n2 = ifelse(length(sub2) < floor(ntheta/2), length(sub2), floor(ntheta/2))
    if((n1 + n2) != ntheta) { 
      if(n1 < n2) { n2 = ntheta - n1 }
      if(n2 < n1) { n1 = ntheta - n2 }
    }
    testtime = sort(c(sample(sub1, n1), sample(sub2, n2)))
    
    ## start loop for time points
    errorlist = list()
    for(j in 1:ntheta) {
      
      # define true theta and ynew
      theta = testtime[j]
  
      ynew = test_dat %>%
        filter(argvals == theta) %>%
        arrange(hormone) %>% pull(y)
      
      # calculate likelihood
      ll = loglik(ynew, hnames, fpca_list)
      mle = ll$thetahat
      
      # calculate error
      daysoff = 28*sign(theta - mle) * min(abs(theta - mle), 1 - abs(theta - mle)) 
      
      errorlist[[j]] = tibble(theta = theta, mle = mle, daysoff = daysoff)
    
    }
    reslist[[i]] = do.call(rbind, errorlist) %>% mutate(testID = IDs[i])
    
  }
  })
  return(do.call(rbind, reslist))
}
  