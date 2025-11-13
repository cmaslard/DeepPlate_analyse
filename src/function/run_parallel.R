run_parallel <- function(group_list, function_to_run,varlist, n_cores = detectCores() - 1,...) {
  os <- Sys.info()[["sysname"]]
  
  
  cat_col(paste0("- Running enrichment with ", n_cores, " cores \n"), "green")
  
  if (os == "Linux" || os == "Darwin") {
    cat_col(paste0("- function run_parallel for ",os,"\n"), "magenta")
    # macOS or Linux: use mclapply()
    all_results <- mclapply(group_list, function_to_run, mc.cores = n_cores, ...)
  } else if (os == "Windows") {
    cat_col(paste0("- function run_parallel for ",os,"\n"), "magenta")
    # Create a cluster for Windows
    cl <- makeCluster(n_cores)
    
    # Export packages and environment to cluster nodes
    #Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = varlist,
                  envir = environment())
    # clusterExport(cl, varlist = "loaded_packages", envir = environment())
    
    # Load packages in cluster nodes
    # clusterEvalQ(cl, {
    #   lapply(loaded_packages, require, character.only = TRUE)
    # })
    
    # Export all necessary objects to cluster nodes
    clusterExport(cl, varlist = ls(envir = environment()), envir = environment())
    
    # Run the function in parallel
    all_results <- parLapply(cl, group_list, function_to_run, ...)
    
    # Stop the cluster
    stopCluster(cl)
  } else {
    stop("Unsupported operating system")
  }
  
  # Remove NULL results
  
  # Return combined results
  
  return(all_results)
}
