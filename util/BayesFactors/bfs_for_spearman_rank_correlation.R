## Bayes factors of Spearman's rank correlation coefficient

library('neatStats')


input_dir <- 'E://dfi_experiment_figures//Paper_figures//iAF//iAF_betw_betabinom//R'


# Try for 1Sound condition bias consistency effect
x <- c(0.0406, 0.6680, -0.0554, -0.3090, -0.2059, -0.3038, 0.4857, -0.1877,
       -0.0746, 0.1938, 0.3333, 0.0540, 0.1402, -0.3878, -0.0006, -0.1540, 0.0484,
       -0.0200, 0.1005, -0.0860)

y <- c(0.1351, 0.7735, -0.1262, -0.0206, 0.1062, -0.1090, 0.1708, -0.0562, -0.2365, 
       -0.0582, -0.0933, 0.2269, -0.2196, -0.0508, 0.0297, -0.1268, -0.0978,
       0.0555, -0.0859, 0.0111)

get_neat_bfs(x, y)



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
