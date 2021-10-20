
# repeating the simulation multiple times and analysing by t.tests instead of lm
simulate_simulations <- function(iterations, nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat){
  analyse_one_exp <- function(nExp, nCell) {  
    temp2 = tibble(model = c("Raw data", 
                             "Normalized to control", 
                             "Normalized to mean of experiment",
                             "ANOVA with Error(Experiment)",
                             "Nested mixed effects model",
                             "Pseudoreplicated") , 
                   t_val = NA, 
                   p_val = NA,
                   difference = NA)
    
    Readout_cell = simulate_multilevel_experiment(nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat) %>% 
      group_by(Experiment) %>% 
      mutate(NormC_Readout = (Readout/mean(Readout[Group=="Ctrl"])),    # normalize to control
             NormS_Readout = Readout - mean(Readout)                    # subtract experiment-mean
      ) %>% 
      ungroup() %>% 
      mutate(NormS_Readout = NormS_Readout + mean(Readout))             # rescale to original exp mean
    
    # Summarize per culture
    Readout_exp = Readout_cell %>%
      group_by(Experiment, Group) %>% 
      summarise(Readout = mean(Readout),
                NormC_Readout = mean(NormC_Readout),
                NormS_Readout = mean(NormS_Readout),
                .groups = "rowwise")
    
    
    # Analysis strategies
    model_exp <- t.test(data = Readout_exp, Readout ~  Group)                 # analyse raw
    model_exp_normC <- t.test(data = Readout_exp, NormC_Readout ~  Group)     # analyse normalized to control
    model_exp_normS <- t.test(data = Readout_exp, NormS_Readout ~  Group)     # analyse normalized by mean of exp
    model_exp_error <- aov(data = Readout_exp, Readout ~  Group + Error(Experiment))  # analyse experiment as error-term
    model_exp_lme4 <- lmerTest::lmer(data = Readout_cell, Readout ~ Group + (1 + Group|Experiment)) # analyse nested mixed model
    model_pseudorep <- t.test(data = Readout_cell, Readout ~  Group)          # analyse with n = cells
    
    # Extracting the relevant info per strategy
    temp2$difference[[1]] <- broom::tidy(model_exp)$estimate2 - broom::tidy(model_exp)$estimate1  # required for sign
    temp2$t_val[[1]] <- abs(broom::tidy(model_exp)$statistic) * sign(temp2$difference[[1]]) 
    temp2$p_val[[1]] <- broom::tidy(model_exp)$p.value * sign(temp2$difference[[1]])     # to have the sign 
    
    
    temp2$difference[[2]] <- broom::tidy(model_exp_normC)$estimate2 - broom::tidy(model_exp_normC)$estimate1  # required for sign
    temp2$t_val[[2]] <- abs(broom::tidy(model_exp_normC)$statistic)  * sign(temp2$difference[[1]])
    temp2$p_val[[2]] <- broom::tidy(model_exp_normC)$p.value * sign(temp2$difference[[2]])     
    
    temp2$difference[[3]] <- broom::tidy(model_exp_normS)$estimate2 - broom::tidy(model_exp_normS)$estimate1  # required for sign
    temp2$t_val[[3]] <- abs(broom::tidy(model_exp_normS)$statistic)  * sign(temp2$difference[[1]])
    temp2$p_val[[3]] <- broom::tidy(model_exp_normS)$p.value * sign(temp2$difference[[3]])     
    
    temp2$difference[[4]] <- coef(model_exp_error)$Within
    temp2$t_val[[4]] <- broom::tidy(model_exp_error)[broom::tidy(model_exp_error)$term == "Group",]$statistic * sign(temp2$difference[[4]])
    temp2$p_val[[4]] <- broom::tidy(model_exp_error)[broom::tidy(model_exp_error)$term == "Group",]$p.value * sign(temp2$difference[[4]])
    
    temp2$difference[[5]] <- coef(summary(model_exp_lme4))[, "Estimate"][2]
    temp2$t_val[[5]] <- coef(summary(model_exp_lme4))[, "t value"][2]
    temp2$p_val[[5]] <- coef(summary(model_exp_lme4))[, "Pr(>|t|)"][2] * sign(temp2$difference[[5]])
    
    temp2$difference[[6]] <- broom::tidy(model_pseudorep)$estimate2 - broom::tidy(model_pseudorep)$estimate1  # required for sign
    temp2$t_val[[6]] <- abs(broom::tidy(model_pseudorep)$statistic) * sign(temp2$difference[[5]]) 
    temp2$p_val[[6]] <- broom::tidy(model_pseudorep)$p.value * sign(temp2$difference[[5]])     # to have the sign 
    
    
    return(temp2)
  }
  
  # repeat 
  tests <- tibble(future_replicate(iterations, analyse_one_exp(nExp, nCell), simplify = F))
  tests <- unnest(tests, cols = c(`future_replicate(iterations, analyse_one_exp(nExp, nCell), simplify = F)`))
  return(tests)
}





# repeating the simulation multiple times and analysing (linear models)
# simulate_simulations_lm <- function(iterations, nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat){
#   a <- function(nExp, nCell) {  
#     temp2 = tibble(model = c("Raw data", 
#                              "Normalized to control", 
#                              "Normalized to mean of experiment",
#                              "ANOVA with Error(Experiment)",
#                              "Nested mixed effects model") , 
#                    t_val = NA, 
#                    difference = NA)
#     
#     Readout_cell = simulate_multilevel_experiment(nExp, nCell, mean1, mean2, sdCell, sdExp, sdTreat) %>% 
#       group_by(Experiment) %>% 
#       mutate(NormC_Readout = (Readout/mean(Readout[Group=="Ctrl"])),    # normalize to control
#              NormS_Readout = Readout - mean(Readout)                    # subtract experiment-mean
#       ) %>% 
#       ungroup() %>% 
#       mutate(NormS_Readout = NormS_Readout + mean(Readout))             # rescale to original exp mean
#     
#     
#     Readout_exp = Readout_cell %>%
#       group_by(Experiment, Group) %>% 
#       summarise(Readout = mean(Readout),
#                 NormC_Readout = mean(NormC_Readout),
#                 NormS_Readout = mean(NormS_Readout),
#                 .groups = "rowwise")
#     
#     model_exp <- lm(data = Readout_exp, Readout ~  Group)                 # analyse raw
#     model_exp_normC <- lm(data = Readout_exp, NormC_Readout ~  Group)     # analyse normalized to control
#     model_exp_normS <- lm(data = Readout_exp, NormS_Readout ~  Group)     # analyse normalized by mean of exp
#     model_exp_error <- aov(data = Readout_exp, Readout ~  Group + Error(Experiment))  # analyse experiment as error-term
#     model_exp_lme4 <- lmerTest::lmer(data = Readout_cell, Readout ~ Group + (1 + Group|Experiment)) # analyse nested mixed model
#     
#     temp2$t_val[[1]] <- broom::tidy(model_exp)[broom::tidy(model_exp)$term == "GroupPRG2",]$statistic
#     temp2$difference[[1]] <- broom::tidy(model_exp)[broom::tidy(model_exp)$term == "GroupPRG2",]$estimate
#     
#     temp2$t_val[[2]] <- broom::tidy(model_exp_normC)[broom::tidy(model_exp_normC)$term == "GroupPRG2",]$statistic
#     temp2$difference[[2]] <- broom::tidy(model_exp_normC)[broom::tidy(model_exp_normC)$term == "GroupPRG2",]$estimate
#     
#     temp2$t_val[[3]] <- broom::tidy(model_exp_normS)[broom::tidy(model_exp_normS)$term == "GroupPRG2",]$statistic
#     temp2$difference[[3]] <- broom::tidy(model_exp_normS)[broom::tidy(model_exp_normS)$term == "GroupPRG2",]$estimate
#     
#     temp2$t_val[[4]] <- broom::tidy(model_exp_error)[broom::tidy(model_exp_error)$term == "Group",]$statistic
#     temp2$difference[[4]] <- coef(model_exp_error)$Within
#     
#     temp2$t_val[[5]] <- coef(summary(model_exp_lme4))[, "t value"][2]
#     temp2$difference[[5]] <- coef(summary(model_exp_lme4))[, "Estimate"][2]
#     
#     
#     return(temp2)
#   }
#   
#   tests <- tibble(future_replicate(iterations, a(nExp, nCell), simplify = F))
#   tests <- unnest(tests, cols = c(`future_replicate(iterations, a(nExp, nCell), simplify = F)`))
#   return(tests)
# }