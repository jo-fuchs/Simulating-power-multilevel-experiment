simulate_individual_experiment <- function(nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat, dist_exp = "normal", dist_cell = "normal") {
  
  
  Readout_cell = simulate_multilevel_experiment(nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat, dist_exp, dist_cell) %>% 
    group_by(Experiment) %>% 
    mutate(NormC_Readout = (Readout/mean(Readout[Group=="Ctrl"])),    # normalize to control
           NormS_Readout = Readout - mean(Readout)                    # subtract experiment-mean
    ) %>% 
    ungroup() %>% 
    mutate(NormS_Readout = NormS_Readout + mean(Readout))             # rescale to original exp mean
  
return(Readout_cell)
  
}
