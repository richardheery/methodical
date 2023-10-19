#' Calculate correlation between expression of transcripts and methylation of sites surrounding their TSS
#'
#' @param meth_rse A RangedSummarizedExperiment for methylation sites. 
#' @param assay_number The assay from meth_rse to extract values from. Default is the first assay. 
#' @param transcript_expression_table A table with the expression values for transcripts, where row names are transcript names and columns sample names. 
#' There should be a row corresponding to each transcript associated with each range in tss_gr. 
#' Names of samples must match those in meth_rse unless samples_subset provided.
#' @param samples_subset Sample names used to subset meth_rse and transcript_expression_table. Provided samples must be found in both meth_rse and transcript_expression_table.
#' Default is to use all samples in meth_rse and transcript_expression_table.
#' @param tss_gr A GRanges object with the locations of transcription start sites (TSS). Each region should have a width of 1. 
#' names(tss_gr) should give the name of the transcript associated with the TSS, which must be present in transcript_expression_table.
#' @param expand_upstream Number of bases to add upstream of each TSS. Must be numeric vector of length 1 or equal to the length of tss_gr. Default is 5000.
#' @param expand_downstream Number of bases to add downstream of each TSS. Must be numeric vector of length 1 or equal to the length of tss_gr. Default is 5000.
#' @param cor_method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor(). Default is "pearson".
#' @param add_distance_to_region A logical value indicating whether to add the distance of methylation sites to the TSS. Default value is TRUE.
#' Setting to FALSE will roughly half the total running time.
#' @param max_sites_per_chunk The approximate maximum number of methylation sites to try to load into memory at once. 
#' The actual number loaded may vary depending on the number of methylation sites overlapping each region, 
#' but so long as the size of any individual regions is not enormous (>= several MB), it should vary only very slightly. 
#' Some experimentation may be needed to choose an optimal value as low values will result in increased running time, 
#' while high values will result in a large memory footprint without much improvement in running time. 
#' Default is floor(62500000/ncol(meth_rse)), resulting in each chunk requiring approximately 500 MB of RAM. 
#' @param ncores Number of cores to use. Default is 1.
#' @param results_dir An optional path to a directory to save results as RDS files. There will be one RDS file for each transcript. 
#' If not provided, returns the results as a list. 
#' @return If results_dir is NULL, a list of data.frames with the correlation of methylation sites surrounding a specified 
#' genomic region with a given feature, p-values and adjusted q-values for the correlations.
#' Distance of the methylation sites upstream or downstream to the center of the region is also provided.
#' If results_dir is provided, instead returns a list with the paths to the RDS files with the results. 
#' @export
#' @examples 
#' 
#' # Load TUBB6 TSS GRanges, RangedSummarizedExperiment with methylation values for CpGs around TUBB6 TSS and TUBB6 transcript counts
#' data(tubb6_tss, package = "methodical")
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#' data(tubb6_transcript_counts, package = "methodical")
#' 
#' # Calculate correlation values between methylation values and transcript values for TUBB6
#' tubb6_cpg_meth_transcript_cors <- methodical::calculate_meth_site_transcript_cors(meth_rse = tubb6_meth_rse, 
#'   transcript_expression_table = tubb6_transcript_counts, tss_gr = tubb6_tss, expand_upstream = 5000, expand_downstream = 5000)
#'   
calculate_meth_site_transcript_cors <- function(meth_rse, assay_number = 1, transcript_expression_table, samples_subset = NULL, tss_gr, expand_upstream = 5000,
  expand_downstream = 5000, cor_method = "pearson", add_distance_to_region = TRUE, max_sites_per_chunk = NULL, ncores = 1, results_dir = NULL){
  
  # Check that all regions in tss_gr have a length of 1 and that they have names
  if(any(width(tss_gr) != 1)){stop("All regions in tss_gr must have a width of 1")}
  if(is.null(names(tss_gr))){stop("names(tss_gr) should be the names of the transcripts associated with each TSS and cannot be NULL")}
     
  # Check that there are no duplicated transcript names in tss_gr or transcript_expression_table
  if(anyDuplicated(names(tss_gr)) | anyDuplicated(row.names(transcript_expression_table))){
    stop("There cannot be duplicated transcript names in names(tss_gr) or row.names(transcript_expression_table)")
  }
  
  # Check that all transcripts associated with tss_gr are present in transcript_expression_table
  # and subset transcript_expression_table for these transcripts if so
  if(!any(names(tss_gr) %in% row.names(transcript_expression_table))){
    stop("There are transcripts in tss_gr that are not present in transcript_expression_table")
  } else {
    transcript_expression_table <- transcript_expression_table[names(tss_gr), ]
  }
  
  # Check that samples_subset are in meth_rse and transcript_expression_table
  if(!is.null(samples_subset)){
    if(any(!samples_subset %in% colnames(meth_rse))){
      stop("Some samples_subset are not in meth_rse")
    } else {meth_rse <- meth_rse[, samples_subset]}
    if(any(!samples_subset %in% names(transcript_expression_table))){
      stop("Some samples_subset are not in transcript_expression_table")
    } else {transcript_expression_table <- dplyr::select(transcript_expression_table, dplyr::all_of(samples_subset))}
  }

  # Check that names of meth_rse and transcript_expression_table match
  if(!all(colnames(meth_rse) == names(transcript_expression_table))){stop(
    "Sample names in meth_rse and transcript_expression_table do not match")
  }
  
  # Check that there are at least three samples and give a warning if there are less than 20 samples
  n_samples <- ncol(meth_rse) 
  if(n_samples < 3){stop("There are not enough samples to calculate correlations")}
  if(n_samples < 20){
    warning(paste("There are only", n_samples, "samples. It is recommended to have at least 20 samples to calculate correlations"))
  }
  
  # Create results_dir is it doesn't exist
  if(!is.null(results_dir)){
    if(!dir.exists(results_dir)){
      dir.create(results_dir)
    } else {
      stop(paste("Directory called", results_dir, "already exists"))
    }
  }
  
  # Print correlation method being used
  message(paste("Using", match.arg(cor_method, c("pearson", "kendall", "spearman")), "correlation method"))

  # Expand tss_gr and subset meth_rse for these regions
  tss_gr_expanded <- GenomicRanges::promoters(tss_gr, expand_upstream, expand_downstream + 1)
  
  # Split tss_gr_expanded into chunks based on the number of methylation sites that they cover
  genomic_region_bins <- .chunk_regions(meth_rse = meth_rse, genomic_regions = tss_gr_expanded, 
    max_sites_per_chunk = max_sites_per_chunk, ncores = 1)
  
  # Create cluster if ncores greater than 1
  cl <- setup_cluster(ncores = ncores, packages = c("methodical", "HDF5Array"), outfile = NULL)
  if(ncores > 1){on.exit(parallel::stopCluster(cl))}
  
  # For each sequence get methylation of the associated regions
  all_correlations <- foreach::foreach(chunk = seq_along(genomic_region_bins)) %do% {
    
    # Print name of sequence which is being summarized
    message(paste("Calculating correlations for chunk", chunk, "of", length(genomic_region_bins)))

    # Get regions from chunk and associated TSS
    regions_for_chunk <- genomic_region_bins[[chunk]]
    tss_for_chunk <- tss_gr[names(regions_for_chunk)]
    
    # Get transcript values for chunk transcripts
    transcript_expression_table_chunk <- transcript_expression_table[names(regions_for_chunk), ]
    
    # Subset meth_rse_for_chunk for regions overlapping regions_for_chunk
    meth_rse_for_chunk <- subsetByOverlaps(meth_rse, regions_for_chunk)
    invisible(gc()) 
    
    # Find the overlaps of regions_for_chunk and meth_rse_for_chunk
    overlaps_df <- data.frame(findOverlaps(regions_for_chunk, meth_rse_for_chunk))
    
    # Add region names to overlaps_df
    overlaps_df$genomic_region_name <- names(regions_for_chunk)[overlaps_df$queryHits]
    
    # Create a list matching region names to rows of meth_rse_for_chunk
    tss_region_indices_list <- split(overlaps_df$subjectHits, overlaps_df$genomic_region_name)
    rm(overlaps_df); gc()
    
    # Put tss_region_indices_list in same order as transcript_expression_table_chunk
    tss_region_indices_list <- tss_region_indices_list[row.names(transcript_expression_table_chunk)]
    
    # Read all values from specified assay of meth_rse_for_chunk into memory and run the garbage collection
    meth_values_chunk <- as.matrix(SummarizedExperiment::assay(meth_rse_for_chunk, i = assay_number))
    gc()
    
    # Add meth site names as row.names
    row.names(meth_values_chunk) <- as.character(SummarizedExperiment::rowRanges(meth_rse_for_chunk))
    
    # Create lists with the methylation values, transcript values and TSS for each transcript
    tss_meth_values <- lapply(tss_region_indices_list, function(x) meth_values_chunk[x, , drop = FALSE])
    transcript_values <- lapply(names(tss_region_indices_list), function(x) transcript_expression_table_chunk[x, , drop = FALSE])
    tss_for_chunk <- split(tss_for_chunk, names(tss_for_chunk))[names(tss_meth_values)] 
    
    # Remove meth_values_chunk and tss_region_indices_list
    rm(meth_values_chunk, tss_region_indices_list)
    
    # Calculate correlations. 
    chunk_correlations <- foreach::foreach(
      meth_table = tss_meth_values, 
      transcript_table = transcript_values,
      transcript_tss = tss_for_chunk, transcript_name = names(tss_meth_values)) %dopar% {
      
      # Transpose meth_table
      meth_table <- t(meth_table)
      
      # Transpose transcript_table and name it with transcript name
      transcript_table <- setNames(data.frame(t(transcript_table)), transcript_name)
        
      tryCatch({
        transcript_meth_site_cors <- methodical::rapid_cor_test(
          table1 = meth_table, table2 = transcript_table, 
          table1_name = "meth_site", table2_name = "transcript_name", 
          cor_method = cor_method, p_adjust_method = "none")
  
          # Add meth site distance to region if specified
          if(add_distance_to_region){transcript_meth_site_cors$distance_to_tss <- 
            methodical::stranded_distance(query_gr = GenomicRanges::GRanges(transcript_meth_site_cors$meth_site), 
              subject_gr = transcript_tss)
          }
        
        # Add transcript TSS as an attribute to correlation data.frame
        attributes(transcript_meth_site_cors)$tss_range <- transcript_tss
        
        # Save results to a RDS file in results_dir if it is provided
        if(!is.null(results_dir)){
          
          # Set name of RDS file as the transcript_name with .rds appended
          rds_file <- paste0(results_dir, "/", transcript_name, ".rds")
          
          # Write the results to the specified RDS file and return the RDS filepath
          saveRDS(transcript_meth_site_cors, rds_file)
          transcript_meth_site_cors <- rds_file
        }
        
        transcript_meth_site_cors
        
      }, error = function(x) data.frame())
    }
  
    # Add names of transcript to list and return
    chunk_correlations <- setNames(chunk_correlations, names(tss_meth_values))
      
    chunk_correlations
  }

  # Run garbage collection one final time and return result
  invisible(gc())
  return(unlist(all_correlations, recursive = FALSE))
  
}
