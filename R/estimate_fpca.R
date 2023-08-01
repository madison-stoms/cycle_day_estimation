
#' @param data a data frame with four arguments: (1) subj: subject ID; (2) hormone: name of hormone; (3) argvals: observation times; (4) y: hormone levels
#' @param hname hormone name 
#' @param newgrid grid on which estimated components will be evaluated
#' @param pve percent variance explained by the included components


estimate_fpca = function(data, hname, newgrid = NULL, pve = 0.998) {
  
  # extract hormone data
  temp = data %>% filter(hormone == hname) %>%
    dplyr::select(subj, argvals, y)
  
  # find common argvals
  if(is.null(newgrid)) { newgrid = sort(unique(temp$argvals)) }
  
  # run face.sparse
  mod = face.sparse(data = temp, argvals.new = newgrid, calculate.scores = TRUE, pve = pve)
  
  return(list(argvals = mod$argvals.new, mu = mod$mu.new, eigenfunctions = mod$eigenfunctions, 
              scores = mod$rand_eff$scores, sigma2 = mod$sigma2, eigenvalues = mod$eigenvalues, npc = mod$npc))
  
}
