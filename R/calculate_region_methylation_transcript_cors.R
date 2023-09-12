#' Calculate the correlation values between the methylation of genomic regions and the expression of associated transcripts
#'
#' @param genomic_regions A GRanges object. 
#' @param genomic_region_names Names for genomic_regions. If not provided, attempts to use genomic_regions$name. 
#' @param genomic_region_transcripts Names of transcripts associated with each region in genomic_regions. 
#' If not provided, attempts to use genomic_regions$transcript_name. All transcripts must be present in transcript_table.
#' @param meth_rse A RangedSummarizedExperiment with methylation values for CpG sites which will be used to calculate methylation values for genomic_regions. 
#' There must be at least 3 samples in common between meth_rse and transcript_table
#' @param transcript_table A table with the expression values for different transcripts in different samples. 
#' Row names should give be the transcript name and column names should be the name of samples. 
#' @param genomic_region_methylation Optional preprovided table with methylation values for genomic_regions 
#' such as created using summarize_region_methylation(). Table will be created if it is not provided which will increase running time.
#' Row names should match genomic_region_names and column names should match those of transcript_table 
#' @param samples_subset Optional sample names used to subset meth_rse and transcript_table. 
#' Provided samples must be found in both meth_rse and transcript_table.
#' Default is to use all samples in meth_rse and transcript_table.
#' @param cor_method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor(). Default is "spearman".
#' @param p_adjust_method Method used to adjust p-values. Same as the methods from p.adjust.methods. Default is Benjamini-Hochberg.
#' @param region_methylation_summary_function A function that summarized column values. Default is colMeans.
#' @param ncores Number of cores to use for parallel processing. Default is 1.
#' @return A data.frame with the correlation values between the methylation of genomic regions and expression of transcripts associated with them
#' @export
calculate_region_methylation_transcript_cors = function(genomic_regions, genomic_region_names = NULL, 
  genomic_region_transcripts = NULL, meth_rse, transcript_table, genomic_region_methylation = NULL,
  samples_subset = NULL, cor_method = "spearman", p_adjust_method = "BH", region_methylation_summary_function = colMeans, ncores = 1){
  
  # Check that samples_subset are in meth_rse and transcript_table
  if(!is.null(samples_subset)){
    if(any(!samples_subset %in% colnames(meth_rse))){
      stop("Some samples_subset are not in meth_rse")
    } else {meth_rse = meth_rse[, samples_subset]}
    if(any(!samples_subset %in% colnames(transcript_table))){
      stop("Some samples_subset are not in transcript_table")
    } else {transcript_table = dplyr::select(transcript_table, dplyr::all_of(samples_subset))}
  }

  # Check that names of meth_rse and transcript_table match
  if(!all(colnames(meth_rse) == names(transcript_table))){stop(
    "Sample names in meth_rse and transcript_table do not match")
  }
  
  # Check that there are at least three samples and give a warning if there are less than 20 samples
  n_samples = ncol(meth_rse) 
  if(n_samples < 3){stop("There are not enough samples to calculate correlations")}
  if(n_samples < 20){
    warning(paste("There are only", n_samples, "samples. It is recommended to have at least 20 samples to calculate correlations"))
  }
  
  # If genomic_region_names not provided, set to genomic_regions$name
  if(is.null(genomic_region_names)){
    genomic_region_names = genomic_regions$name
  }
  
  # If genomic_region_transcripts not provided, set to genomic_regions$transcript_name
  if(is.null(genomic_region_transcripts)){
    genomic_region_transcripts = genomic_regions$transcript_name
  }
  
  # Check that all genomic_region_transcripts are in transcript_table
  if(!any(genomic_region_transcripts %in% row.names(transcript_table))){
    stop("All genomic_region_transcripts must be present in transcript table")
  }
  
  # Check that there are no duplicates for genomic_region_names and that it has the same length as genomic_regions
  if(length(genomic_region_names) != length(genomic_regions)){stop("genomic_region_names must have the same length as genomic_regions")}
  if(anyDuplicated(genomic_region_names)){stop("genomic_region_names cannot contain duplicates")}
  
  # If genomic_region_methylation provided, check that it has the correct row and column names. 
  if(!is.null(genomic_region_methylation)){
    if(any(row.names(genomic_region_methylation) != genomic_region_names)){
      stop("Row names of genomic_region_methylation do not match genomic_region_names")
    }
    if(any(colnames(genomic_region_methylation) != colnames(transcript_table))){
      if(is.null(samples_subset)){
        stop("Column names of genomic_region_methylation do not match those of transcript table")
      } else {
        stop("Column names of genomic_region_methylation do not match samples_subset")
      }
    }
  }
  
  # Print correlation mthod being used
  message(paste("Using", match.arg(cor_method, c("pearson", "kendall", "spearman")), "correlation method"))
  
  # Create genomic_region_methylation if it is not provided
  if(is.null(genomic_region_methylation)){
  
    # Summarize methylation values for genomic_regions
    cat("Summarizing region methylation\n")
    genomic_region_methylation = methodical::summarize_region_methylation(
      meth_rse = meth_rse, genomic_regions = genomic_regions, 
      genomic_regions_names = genomic_region_names, 
      summary_function = region_methylation_summary_function, n_chunks_parallel = ncores
    )
    
    # Set region_name of genomic_region_methylation as row.names
    genomic_region_methylation = tibble::column_to_rownames(genomic_region_methylation, "region_name")
    
  }
  
  # Calculate correlations between region methylation values and transcript expression
  cat("Calculating transcript correlations")
  
  # Create a data.frame matching genomic_region_names and transcript_names
  feature_matches_df = 
    data.frame(genomic_region_names = genomic_region_names, transcript_name = genomic_region_transcripts)
  
  # Create a list of the transcripts associated with each genomic region
  feature_matches = split(feature_matches_df$genomic_region_names, feature_matches_df$transcript_name)
  
  # Create cluster if ncores greater than 1
  if(ncores > 1){
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, ncores)
    `%dopar%` = foreach::`%dopar%`
    `%do%` = foreach::`%do%`
  } else {
    `%dopar%` = foreach::`%do%`
    `%do%` = foreach::`%do%`
  }
  
  # For each transcript, calculate the correlation between its expression and the methylation of regions associated with it
  methylation_transcript_correlations = foreach::foreach(transcript = names(feature_matches)) %dopar% {
    
    methodical:::rapid_cor_test(
      table1 = t(genomic_region_methylation[feature_matches[[transcript]], ]), 
      table2 = setNames(data.frame(unlist(transcript_table[transcript, ])), transcript),
      table1_name = "genomic_region_name", table2_name = "transcript_name",
      cor_method = cor_method, p_adjust_method = "none")
  }
  
  # Stop the cluster if present
  if(ncores > 1){parallel::stopCluster(cl)}
  
  # Combine results into a single table
  methylation_transcript_correlations = dplyr::bind_rows(methylation_transcript_correlations)
  
  # Adjust p-values
  methylation_transcript_correlations$q_val = p.adjust(methylation_transcript_correlations$p_val, method = p_adjust_method)
  
  # Put methylation_transcript_correlations in order of genomic_region_names
  methylation_transcript_correlations = methylation_transcript_correlations[match(
    genomic_region_names, methylation_transcript_correlations$genomic_region_name), ]
  row.names(methylation_transcript_correlations) = NULL
  
  # Return results
  return(methylation_transcript_correlations)
  
}