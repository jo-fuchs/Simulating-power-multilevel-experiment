
plot_individual_sim <- function(data, readout) { 
  
  Readout_exp = data %>%
    group_by(Experiment, Group) %>% 
    summarise(Readout = mean(Readout),
              NormC_Readout = mean(NormC_Readout),
              NormS_Readout = mean(NormS_Readout),
              .groups = "rowwise") %>% 
    ungroup()
  
  sdCell_measured <- data %>% summarise(round(sd({{readout}}), 2))
  sdExp_measured <- Readout_exp %>% summarise(round(sd({{readout}}),2))
  
  plot1 <- ggplot(data, aes(x = Group, y = {{readout}}, color = factor(Experiment))) + 
    ggbeeswarm::geom_quasirandom() +
    facet_wrap(~Experiment, ncol = 6) +            
    
    stat_summary(fun = mean, geom = "crossbar", color = "black", fill = NA, size = 0.5) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.5, size = 0.5) +
    
    labs(title = "Individual cultures") +
    coord_cartesian(ylim = c(0, max(data$Readout))) +
    ggsci::scale_color_npg()+
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none", 
          strip.text = element_blank(), axis.title.x = element_blank())
  
  plot2 <- ggplot(Readout_exp, aes(x = Group, y = {{readout}}, color = factor(Experiment))) + 
    ggbeeswarm::geom_quasirandom(data = data, alpha = 0.3, shape = 16) +
    geom_point(data = Readout_exp, size = 4)+
    
    stat_summary(fun = mean, geom = "crossbar", color = "black", fill = NA, size = 0.5) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", color = "black", width = 0.5, size = 0.5) +
    
    labs(color = "Culture",
         title = "Total experiment") +
    coord_cartesian(ylim = c(0, max(data$Readout))) +
    ggsci::scale_color_npg()+
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())
  
  figures <- cowplot::plot_grid(plot1, plot2, ncol = 2, rel_widths = c(2,1), 
                                align = "h", axis = "tb")
  
  subtitleplot = ggplot() +
    labs(caption = glue::glue("measured sdCell = {sdCell_measured} measured, sdExp = {sdExp_measured}\nsdCell = {sdCell}, sdExp = {sdExp}, sdTreat = {sdTreat}, mean1 = {mean1}, mean2 = {mean2}"))
  
  cowplot::plot_grid(figures, subtitleplot, nrow = 2, rel_heights = c(10,1))
}

