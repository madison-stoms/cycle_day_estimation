
#' @param data a data frame with four arguments: (1) subj: subject ID; (2) hormone: name of hormone; (3) argvals: observation times; (4) y: hormone levels
#' @param hnames vector of hormone names to be used in the analysis
#' @param newgrid grid on which estimated components will be evaluated
#' @param pve percent variance explained by the included components


create_fpca_list = function(data, hnames, newgrid = NULL, pve = 0.998) {
  
  fpca_list = hnames %>%
    map( ~ estimate_fpca(data, hname = .x, newgrid, pve = pve))
  
  names(fpca_list) = hnames
  
  return(fpca_list)
  
} 