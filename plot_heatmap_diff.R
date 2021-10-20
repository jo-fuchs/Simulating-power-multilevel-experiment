# consider making mean1, mean2 not global

plot_heatmap_diff <- function(results_matrix, analysis_method, max_change, only_significant = F) {
   
  # difference of estimator and real difference -> rescale to percent change
  if(analysis_method == "matrix_normC") {
    # rescale to original scale, then to percent scale
    long_mat <- reshape2::melt(results_matrix[[analysis_method]]$diff)  %>% 
      mutate(norm_value = 100*(value*mean1/(mean2-mean1) - 1),
             norm_value = if_else(norm_value > max_change, max_change, norm_value),
             norm_value = if_else(norm_value < max_change * -1, max_change* -1, norm_value),
             significant = abs(reshape2::melt(results_matrix[[analysis_method]]$p_val)$value) < 0.05)
  } else {
    long_mat <- reshape2::melt(results_matrix[[analysis_method]]$diff)  %>% 
      mutate(norm_value = 100*((value/(mean2-mean1) - 1)),
             norm_value = if_else(norm_value > max_change, max_change, norm_value),
             norm_value = if_else(norm_value < max_change * -1, max_change* -1, norm_value),
             significant = abs(reshape2::melt(results_matrix[[analysis_method]]$p_val)$value) < 0.05)
  }
  
  
  
  ggplot(long_mat, aes(x=Var1, y=Var2, fill = norm_value)) + 
    # only plot those that are significant
    {if(only_significant) geom_raster(aes(alpha = significant))} +
    {if(!only_significant) geom_raster()} +
    labs(x = "Analysed cells per experiment",
         y = "Analysed experiments",
         fill = "% Estimate error",
         title = glue::glue("Simulating experiments with large cell-to-cell variability,\n{results_matrix[[analysis_method]]$name}"),
         subtitle = "blue = overestimating difference, yellow = equivalence, red = underestimating",
         caption = paste("Mean1 = ", mean1, ", Mean2 = ", mean2, ", sdCell = ", sdCell, ", sdExp = ", sdExp, ", sdTreat = ", sdTreat)) + 
    
    scale_fill_distiller(type = "div", limits = c(-max_change, max_change), 
                         palette = "RdYlBu", direction = 1)+
    
    scale_x_continuous(breaks = seq(ceiling(sqrt(simulations)/2), maxCell/resCell*sqrt(simulations), sqrt(simulations)),
                       labels = seq(resCell, maxCell, resCell)) +
    scale_y_continuous(breaks = seq(ceiling(sqrt(simulations)/2), maxExp/resExp*sqrt(simulations), sqrt(simulations)),
                       labels = seq(resExp, maxExp, resExp)) +
    
    geom_vline(xintercept=seq(0.5, maxCell/resCell*sqrt(simulations), sqrt(simulations)), color = "grey") +
    geom_hline(yintercept=seq(0.5, maxExp/resExp*sqrt(simulations), sqrt(simulations)), color = "grey") + 
    
    theme_minimal() + theme(panel.grid = element_blank())
}

