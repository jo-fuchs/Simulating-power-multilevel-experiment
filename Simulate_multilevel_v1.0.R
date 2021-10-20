# just run all analyses on the same data

library(tidyverse)
library(future.apply)


source("simulate_multilevel_experiment.R")
source("simulate_individual_experiment.R")
source("simulate_simulations.R")
source("heatmapify_simulations.R")
source("plot_heatmap_tp.R")
source("plot_heatmap_diff.R")
source("plot_individual_simulation.R")




#### Define variables

simulations =  2^2
maxExp = 20
maxCell = 50
resExp = 2
resCell = 5
mean1 = 0.18
mean2 = 0.21
sdCell = 0.05
sdExp = 0.02
sdTreat = 0.02

one_experiment <- simulate_individual_experiment(nExp = 6, nCell = 15, mean1, mean2, sdCell, sdExp, sdTreat)
plot_individual_sim(one_experiment, Readout)
plot_individual_sim(one_experiment, NormC_Readout*mean1) 
plot_individual_sim(one_experiment, NormS_Readout)
plot_individual_sim(one_experiment, Readout)


## start simulating lots of combinations -> this takes time
plan(multisession, workers = 8)
results_matrix <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell, sdExp, sdTreat)


plot_heatmap_p(results_matrix, "matrix_raw", alpha = 0.05, direction = 1)
plot_heatmap_p(results_matrix, "matrix_normC", alpha = 0.05, direction = 1)
plot_heatmap_p(results_matrix, "matrix_normS", alpha = 0.05, direction = 1)
plot_heatmap_p(results_matrix, "matrix_error", alpha = 0.05, direction = 1)
plot_heatmap_p(results_matrix, "matrix_mixed", alpha = 0.05, direction = 1)
plot_heatmap_p(results_matrix, "matrix_pseudorep", alpha = 0.05, direction = 1)



plot_heatmap_diff(results_matrix, "matrix_raw", max_change = 100, only_significant = F)
plot_heatmap_diff(results_matrix, "matrix_normC", 100, only_significant = F)
plot_heatmap_diff(results_matrix, "matrix_normS", 100)
plot_heatmap_diff(results_matrix, "matrix_error", 100)
plot_heatmap_diff(results_matrix, "matrix_mixed", 100)
plot_heatmap_diff(results_matrix, "matrix_pseudorep", max_change = 100, only_significant = T)





# other scenarios (starting from the best-case high between exp, low treatment variability)
results_matrix_noDiff <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean1, sdCell, sdExp, sdTreat)
results_matrix_hiDiff <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2+0.03, sdCell, sdExp, sdTreat)
results_matrix_loCell <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell/2, sdExp, sdTreat)
results_matrix_loExp <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell, sdExp/2, sdTreat)
results_matrix_higTreatVar <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell, sdExp/2, sdTreat*2)

mean1 = 0.6
mean2 = 0.4
sdCell = 0.07
sdExp = 0.04
sdTreat = 0.02

results_matrix_mem <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell, sdExp, sdTreat)
results_matrix_mem_loCell <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell/2, sdExp, sdTreat)
results_matrix_mem_hiCell <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell*2, sdExp, sdTreat)
results_matrix_mem_loExp <- heatmapify_simulations(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell, sdExp, sdTreat)



