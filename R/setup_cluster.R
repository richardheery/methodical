#' Setup a cluster for use with doParallel
#'
#' @param ncores Number of cores. 
#' @param packages A character vector with the names of packages to load.
#' @param outfile Outfile for the cluster. Default is the NULL device.  
#' @return A SOCKcluster if ncores is greater than 1 or NULL otherwise.
setup_cluster = function(ncores, packages = NULL, outfile = NULL){
  
  # Define %do% and %dopar%
  assign("%do%", foreach::`%do%`, envir = parent.frame())
  assign("%dopar%", foreach::`%dopar%`, envir = parent.frame())
  
  # Create cluster, register it and load specified packages if ncores greater than 1
  if(ncores > 1){
    if(!is.null(outfile)){
      cl = parallel::makeCluster(ncores, outfile = outfile)
    } else {
      cl = parallel::makeCluster(ncores)
    }
    doParallel::registerDoParallel(cl, ncores)
    if(!is.null(packages)){
      lapply(packages, function(x)
        parallel::clusterCall(cl,
          function(x){invisible(suppressMessages(library(x, character.only = TRUE)))}, x))
    }
    # Return the cluster
    return(cl)
  } else {
    foreach::registerDoSEQ()
    return(NULL)
  }
  
}