## Bayes factors of Spearman's rank correlation coefficient

library('neatStats')


input_dir <- 'E://dfi_experiment_figures//Paper_figures//iAF//iAF_betw_betabinom//R'


# Helper function
get_neat_bfs <- function(var1, var2) {
  # Compute parametric and non-parametric Bayes factors for
  # correlation between var1 and var2
  
  print('Parametric correlation')
  bf_parametric = corr_neat(
    var1,
    var2,
    nonparametric = FALSE,
    bf_added = TRUE,
    direction = NULL,
  )
  
  print('Non parametric correlation')
  bf_non_parametric = corr_neat(
    var1,
    var2,
    nonparametric = TRUE,
    bf_added = TRUE,
    direction = NULL,
  )
  
  return(list(bf_parametric, bf_non_parametric))
}


# Loop through files (tasks and conditions)
tasks <- c('2ifc', 'yesno', 'yn_threshold')
for (itask in c(1, 2, 3)) {
  
  print(sprintf('================= %s =================', tasks[itask]))
  
  for (icond in c(1, 2, 3)) {
    
    print(sprintf('Condition = %i', icond))
    
    filename <- sprintf('%s_%i.csv', tasks[itask], icond)
    data <- read.csv(paste(input_dir, filename, sep="//"))
  
    outp <- get_neat_bfs(data[['Var1']], data[['Var2']])
    
    print(sprintf('Parametric BF = %f, Non-parametric BF = %f', 
                  outp[[1]][3], outp[[2]][3]))
  }
}
