

# plotting p-value
plot_heatmap_p <- function(results_matrix, analysis_method, alpha, direction) {
  
  # set axis labels for P-value
  scale_p_lab <- c("-0.001","-0.01", "-0.05", "-0.2", "1", "0.2", "0.05", "0.01", "0.001")
  scale_p_val <- c(-3, -2, -1.3, -0.68 , 0, 0.68, 1.3, 2, 3)
  
  
  long_mat <- reshape2::melt(results_matrix[[analysis_method]]$p_val) %>% 
    mutate(log_value = log10(abs(1/value))*sign(value),
           log_value = if_else(log_value > log10(1/alpha), log10(1/alpha), log_value),
           log_value = if_else(log_value < log10(1/alpha)*-1, log10(1/alpha)*-1, log_value))
  
  ggplot(long_mat, aes(x=Var1, y=Var2, fill = log_value*direction)) + geom_raster() +
    labs(x = "Analysed cells per experiment",
         y = "Analysed experiments",
         fill = "p-value",
         title = glue::glue("Simulating experiments with large cell-to-cell variability,\n{results_matrix[[analysis_method]]$name}"),
         subtitle = "dark-green = significant, white = not significant, pink = significant wrong direction",
         caption = paste("Mean1 = ", mean1, ", Mean2 = ", mean2, ", sdCell = ", sdCell, ", sdExp = ", sdExp, ", sdTreat = ", sdTreat)) + 
    
    scale_fill_distiller(type = "div", 
                         limits = c(-log10(1/alpha),log10(1/alpha)),
                         breaks = scale_p_val,
                         labels = scale_p_lab,
                         palette = "PiYG", direction = 1)+
    
    scale_x_continuous(breaks = seq(ceiling(sqrt(simulations)/2), maxCell/resCell*sqrt(simulations), sqrt(simulations)),
                       labels = seq(resCell, maxCell, resCell)) +
    scale_y_continuous(breaks = seq(ceiling(sqrt(simulations)/2), maxExp/resExp*sqrt(simulations), sqrt(simulations)),
                       labels = seq(resExp, maxExp, resExp)) +
    
    geom_vline(xintercept=seq(0.5, maxCell/resCell*sqrt(simulations), sqrt(simulations)), color = "white") +
    geom_hline(yintercept=seq(0.5, maxExp/resExp*sqrt(simulations), sqrt(simulations)), color = "white") + 
    theme_minimal() + theme(panel.grid = element_blank())
}



# 
# # plotting t-Statistic (not correct for low n)
# plot_heatmap_t <- function(results_matrix, analysis_method, minmax, direction) {
#   
#   long_mat <- reshape2::melt(results_matrix[[analysis_method]]$t_val) # difference of estimator and real difference
#   
#   ggplot(long_mat, aes(x=Var1, y=Var2, fill = value*direction)) + geom_raster() +
#     labs(x = "Analysed cells per experiment",
#          y = "Analysed experiments",
#          fill = "T-statistic",
#          title = glue::glue("Simulating experiments with large cell-to-cell variability,\n{results_matrix[[analysis_method]]$name}"),
#          subtitle = "dark-green = significant, white = equivalence, pink = wrong direction",
#          caption = paste("Mean1 = ", mean1, ", Mean2 = ", mean2, ", sdCell = ", sdCell, ", sdExp = ", sdExp, ", sdTreat = ", sdTreat)) + 
#     
#     scale_fill_distiller(type = "div", limits = minmax, palette = "PiYG", direction = 1, na.value = "#238b45")+
#     scale_x_continuous(breaks = seq(ceiling(sqrt(simulations)/2), maxCell/resCell*sqrt(simulations), sqrt(simulations)),
#                        labels = seq(resCell, maxCell, resCell)) +
#     scale_y_continuous(breaks = seq(ceiling(sqrt(simulations)/2), maxExp/resExp*sqrt(simulations), sqrt(simulations)),
#                        labels = seq(resExp, maxExp, resExp)) +
#     geom_vline(xintercept=seq(0.5, maxCell/resCell*sqrt(simulations), sqrt(simulations)), color = "white") +
#     geom_hline(yintercept=seq(0.5, maxExp/resExp*sqrt(simulations), sqrt(simulations)), color = "white") + 
#     theme_minimal() + theme(panel.grid = element_blank())
# }
# 
# 
