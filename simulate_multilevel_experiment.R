# 
# sdCell: cell to cell variability
# sdExp: Experiment to Experiment variability
# sdTreatment: Treatment stability (varying slope between experiments)

# all normally distributed
simulate_multilevel_experiment <- function(nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat, dist_exp, dist_cell) {
  
  a<- function(nExp, nCell){
    temp <- data.frame(Group = c("Ctrl", "Treat"))
    temp$Dummy <-c(0, 1)
    temp$Readout <- (mean1 + rnorm(1, mean =0, sd = sdExp)) + # culture to culture variability
      (mean2-mean1 +                                    # treatment effect
         rnorm(1, mean =0, sd = sdTreat))*(temp$Dummy) # treatment variability (-> differing at culture level)
    
    
    temp2 <- data.frame(Group = c(rep("Ctrl", nCell), rep("Treat", nCell)))
    temp2$Dummy <-c(rep(0, nCell), rep(1, nCell))
    
    temp2$Readout <- c(rep(temp$Readout[1], nCell), 
                       rep(temp$Readout[2], nCell)) + 
      rnorm(nCell*2, mean =0, sd = sdCell)
    return(temp2)
  }  
  
  true_readout_experiment = tibble(rerun(nExp, a(nExp, nCell))) %>% mutate(Experiment = row_number())
  true_readout_experiment <-unnest(true_readout_experiment,cols = c(`rerun(nExp, a(nExp, nCell))`))
  return(true_readout_experiment)
}

# 
# # all lognormal data
# simulate_multilevel_experiment <- function(nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat) {
# 
#   a<- function(nExp, nCell){
#     temp <- data.frame(Group = c("Ctrl", "PRG2"))
#
#     # convert SDs and means
#     # https://blogs.sas.com/content/iml/2014/06/04/simulate-lognormal-data-with-specified-mean-and-variance.html
#     log_sdExp = sqrt(log(1 + sdExp^2/mean1^2))
#     log_sdTreat = sqrt(log(1 + (sdTreat)^2 / mean2^2))
#     log_mean1 = log(mean1^2 / sqrt(sdExp^2 + mean1^2))
#     log_mean2 = log((mean1-mean2)^2 / sqrt((sdTreat)^2 + (mean1-mean2)^2))
# 
#     mean_exp_control <- exp(rnorm(1, mean = log_mean1, sd = log_sdExp ))
#     
#     temp$Readout <- c(mean_exp_control,
#                       mean_exp_control + exp(rnorm(1, mean = log_mean2, sd = log_sdTreat)))
# 
#     temp2 <- data.frame(Group = c(rep("Ctrl", nCell), rep("PRG2", nCell)))
# 
#     mean_cell = log(c(rep(temp$Readout[1], nCell),
#                       rep(temp$Readout[2], nCell)))
# 
#     # convert SDcell (not sure why exp but it makes sds more comparable to normal SD)
#     log_sdCell = sqrt(log(1 + sdCell^2 / mean_cell^2))
# 
#     temp2$Readout <- exp(rnorm(nCell*2, mean = mean_cell, sd = log_sdCell))
#     return(temp2)
#   }
# 
#   true_readout_experiment = tibble(rerun(nExp, a(nExp, nCell))) %>% mutate(Experiment = row_number())
#   true_readout_experiment <-unnest(true_readout_experiment,cols = c(`rerun(nExp, a(nExp, nCell))`))
#   return(true_readout_experiment)
# }

