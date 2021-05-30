## Bayes factors of Spearman's rank correlation coefficient

library('neatStats')


# ------------------------------------------------------
#
#         Between subject frequency analysis
#
# ------------------------------------------------------

# Select input directory

# iAF ~ threshold; eyes-open
input_dir <- 'E://dfi_experiment_figures//Paper_figures//iAF//iAF_betw_betabinom//R'

# iAF ~ threshold; eyes-closed
input_dir <- 'E://dfi_experiment_figures//Paper_figures//iAF//iAF_betw_betabinom//R//eyes-closed'

# iAF ~ threshold; source space (lcmv)
input_dir <- 'E://dfi_experiment_figures//Paper_figures//iAF//iAF_betw_betabinom//R//lcmv'


# Output directory will be the same as input directory
output_dir <- input_dir


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
for (itask in 1:3) {
  
  print(sprintf('================= %s =================', tasks[itask]))
  
  bfs_mat = matrix(data=NA, nrow=10, ncol=3)
  for (irep in 1:10) {
    
    bfs_for_matlab = c()
    for (icond in 1:3) {
      
      print(sprintf('Condition = %i', icond))
      
      # input
      filename <- sprintf('%s_%i.csv', tasks[itask], icond)
      data <- read.csv(paste(input_dir, filename, sep="//"))
    
      # compute bfs
      outp <- get_neat_bfs(data[['Var1']], data[['Var2']])
      
      print(sprintf('Parametric BF = %f, Non-parametric BF = %f', 
                    outp[[1]][3], outp[[2]][3]))
      
      # output 1 iteration
      bfs_for_matlab <- c(bfs_for_matlab, outp[[2]][[3]])
    }
    
    # output all iterations
    bfs_mat[irep, ] <- bfs_for_matlab
  }
  
  # save output for Matlab
  df <- data.frame(colMeans(bfs_mat))
  filename_out <- sprintf('%s_non_param_bfs.csv', tasks[itask])
  write.csv(df, paste(output_dir, filename_out, sep="//"))
}



# ------------------------------------------------------
#
#         Within subject frequency analysis
#           Effect consistency over tasks
#
# ------------------------------------------------------

# todo...
