
heatmapify_simulations <- function(maxExp, maxCell, resExp, resCell, simulations, mean1, mean2, sdCell, sdExp, sdTreat) {
  
  # pre-assign lists of matrices
  matrix_raw   = list(name = NA,
                      t_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      p_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      diff = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)))
  matrix_normC = list(name = NA,
                      t_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      p_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      diff = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)))
  matrix_normS = list(name = NA,
                      t_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      p_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      diff = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)))
  matrix_error = list(name = NA,
                      t_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      p_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      diff = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)))
  matrix_mixed = list(name = NA,
                      t_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      p_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      diff = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)))  

  matrix_pseudorep = list(name = NA,
                      t_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      p_val = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)),
                      diff = matrix(ncol = maxExp/resExp*sqrt(simulations), nrow = maxCell/resCell*sqrt(simulations)))  
  
  matrix_store_results <- function(matrix_input, model_name, simulation_results_tb, row, col, simulations) {
    tests2 <- simulation_results_tb %>% filter(model == model_name)
    
    matrix_input$name = model_name
    matrix_input$t_val[((row-1)*sqrt(simulations)+1):(row*sqrt(simulations)), 
                       ((col-1)*sqrt(simulations)+1):(col*sqrt(simulations))] = matrix(tests2$t_val, ncol = sqrt(simulations), nrow = sqrt(simulations))
    matrix_input$diff[((row-1)*sqrt(simulations)+1):(row*sqrt(simulations)), 
                      ((col-1)*sqrt(simulations)+1):(col*sqrt(simulations))] = matrix(tests2$difference, ncol = sqrt(simulations), nrow = sqrt(simulations))
    matrix_input$p_val[((row-1)*sqrt(simulations)+1):(row*sqrt(simulations)), 
                       ((col-1)*sqrt(simulations)+1):(col*sqrt(simulations))] = matrix(tests2$p_val, ncol = sqrt(simulations), nrow = sqrt(simulations))
    
    return(matrix_input)
  }
  
  model_names <-c("Raw data", 
                  "Normalized to control", 
                  "Normalized to mean of experiment",
                  "ANOVA with Error(Experiment)",
                  "Nested mixed effects model",
                  "Pseudoreplicated")
  
  for (col in 1:(maxExp/resExp)) {
    print(paste("Experiments: ",col*resExp))
    for (row in 1:(maxCell/resCell)) {
      print(row*resCell)
      tests <- simulate_simulations(simulations, col*resExp, row*resCell, mean1, mean2, sdCell, sdExp, sdTreat)
      
      matrix_raw = matrix_store_results(matrix_raw, model_names[1], tests, row, col, simulations)
      matrix_normC = matrix_store_results(matrix_normC, model_names[2], tests, row, col, simulations)
      matrix_normS = matrix_store_results(matrix_normS, model_names[3], tests, row, col, simulations)
      matrix_error = matrix_store_results(matrix_error, model_names[4], tests, row, col, simulations)
      matrix_mixed = matrix_store_results(matrix_mixed, model_names[5], tests, row, col, simulations)
      matrix_pseudorep = matrix_store_results(matrix_pseudorep, model_names[6], tests, row, col, simulations)
      
    }
    
  }
  result = list(matrix_raw   = matrix_raw, 
                matrix_normC = matrix_normC , 
                matrix_normS = matrix_normS, 
                matrix_error = matrix_error, 
                matrix_mixed = matrix_mixed,
                matrix_pseudorep = matrix_pseudorep)
  return(result)
}
