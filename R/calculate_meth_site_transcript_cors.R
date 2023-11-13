#' Create an iterator function for use with bpiterate
#'
#' @param meth_values_chunk A table with methylation values for current chunk
#' @param tss_region_indices_list A list with the indices for methylation sites associated with each TSS.
#' @param transcript_values A list with expression values for transcripts.
#' @param tss_for_chunk A list of GRanges with the TSS for the current chunk.
#' @param cor_method Correlation method to use. 
#' @param add_distance_to_region Logical value indicating whether to add distance to TSS.
#' @param results_dir Location of results directory. 
#' @return An iterator function which returns a list with the parameters necessary for .tss_correlations. 
.tss_iterator = function(meth_values_chunk, tss_region_indices_list, transcript_values, tss_for_chunk, 
  cor_method, add_distance_to_region, results_dir){
  
  n <- length(tss_for_chunk)
  i <- 0L
  function() {
    if(i >= n)
      NULL
    else {
      i <<- i + 1L
      list(
        meth_table = meth_values_chunk[tss_region_indices_list[[i]], , drop = FALSE], 
        transcript_table = transcript_values[[i]],
        transcript_tss = tss_for_chunk[[i]], 
        transcript_name = names(tss_for_chunk)[i],
        cor_method = cor_method, add_distance_to_region = add_distance_to_region, results_dir = results_dir
        )
    }
  }
}

#' Calculate meth site-transcript correlations for given TSS
#'
#' @param correlation_objects A list with a table of methylation values, 
#' expression values for transcripts, a GRanges for the TSS and the name of the transcript. 
#' @return A data.frame with the correlation values
.tss_correlations = function(correlation_objects){
      
  meth_table <- correlation_objects[["meth_table"]]
  transcript_table <- correlation_objects[["transcript_table"]]
  transcript_tss <- correlation_objects[["transcript_tss"]]
  transcript_name <- correlation_objects[["transcript_name"]]
  cor_method <- correlation_objects[["cor_method"]]
  add_distance_to_region <- correlation_objects[["add_distance_to_region"]]
  results_dir <- correlation_objects[["results_dir"]]
      
    # Transpose meth_table
    meth_table <- t(meth_table)
    
    # Transpose transcript_table and name it with transcript name
    transcript_table <- setNames(data.frame(t(transcript_table)), transcript_name)
        
    transcript_meth_site_cors <- tryCatch({
      transcript_meth_site_cors <- methodical::rapidCorTest(
        table1 = meth_table, table2 = transcript_table, 
        table1_name = "meth_site", table2_name = "transcript_name", 
        cor_method = cor_method, p_adjust_method = "none")
  
        # Add meth site distance to region if specified
        if(add_distance_to_region){transcript_meth_site_cors$distance_to_tss <- 
          methodical::strandedDistance(query_gr = GenomicRanges::GRanges(transcript_meth_site_cors$meth_site), 
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
        
    
    # Return transcript_meth_site_cors or else an emoty data.frame if there was an error
    transcript_meth_site_cors}, error = function(error) data.frame())
    
}

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
#' @param BPPARAM A BiocParallelParam object for parallel processing. Defaults to `BiocParallel::bpparam()`. 
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
#' tubb6_cpg_meth_transcript_cors <- methodical::calculateMethSiteTranscriptCors(meth_rse = tubb6_meth_rse, 
#'   transcript_expression_table = tubb6_transcript_counts, tss_gr = tubb6_tss, expand_upstream = 5000, expand_downstream = 5000)
#'   
calculateMethSiteTranscriptCors <- function(meth_rse, assay_number = 1, transcript_expression_table, 
  samples_subset = NULL, tss_gr, expand_upstream = 5000, expand_downstream = 5000, cor_method = "pearson", 
  add_distance_to_region = TRUE, max_sites_per_chunk = NULL, BPPARAM = BiocParallel::bpparam(), results_dir = NULL){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), is(assay_number, "numeric"),
    is(transcript_expression_table, "data.frame") | is(transcript_expression_table, "matrix"),
    is(samples_subset, "character") | is.null(samples_subset), 
    is(tss_gr, "GRanges"), is(expand_upstream, "numeric"), is(expand_downstream, "numeric"),
    is(cor_method, "character"), is(add_distance_to_region, "logical"),
    is(max_sites_per_chunk, "numeric") & max_sites_per_chunk >= 1 | is.null(max_sites_per_chunk), 
    is(BPPARAM, "BiocParallelParam"), is(results_dir, "character") | is.null(results_dir))
  cor_method = match.arg(cor_method, c("pearson", "kendall", "spearman"))
  
  # Check that all regions in tss_gr have a length of 1 and that they have names
  if(any(width(tss_gr) != 1)){stop("All regions in tss_gr must have a width of 1")}
  if(is.null(names(tss_gr))){stop("names(tss_gr) should be the names of the transcripts associated with each TSS and cannot be NULL")}
     
  # Check that there are no duplicated transcript names in tss_gr or transcript_expression_table
  if(anyDuplicated(names(tss_gr)) | anyDuplicated(row.names(transcript_expression_table))){
    stop("There cannot be duplicated transcript names in names(tss_gr) or row.names(transcript_expression_table)")
  }
  
  # Identify transcripts in common between tss_gr and transcript_expression_table.
  # Throw an error if there are no common transcripts and subset tss_gr and transcript_expression_table
  common_transcripts = intersect(names(tss_gr), row.names(transcript_expression_table))
  if(length(common_transcripts) == 0){
    stop("There are no common transcripts/genes between names(tss_gr) and row.names(transcript_expression_table)")
  } else {
    message(paste("There are", length(common_transcripts), "genes/transcripts in common between tss_gr and transcript_expression_table"))
  }
  transcript_expression_table <- transcript_expression_table[common_transcripts, ]
  tss_gr <- tss_gr[common_transcripts]
  
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
    message(paste("There are only", n_samples, "samples. It is recommended to have at least 20 samples to calculate correlations"))
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
  genomic_region_bins <- methodical:::.chunk_regions(meth_rse = meth_rse, genomic_regions = tss_gr_expanded, 
    max_sites_per_chunk = max_sites_per_chunk, ncores = 1)
  
  # For each sequence get methylation of the associated regions
  `%do%` <- foreach::`%do%`
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
    transcript_values <- lapply(names(tss_region_indices_list), function(x) transcript_expression_table_chunk[x, , drop = FALSE])
    tss_for_chunk <- split(tss_for_chunk, names(tss_for_chunk))[names(tss_region_indices_list)] 
    
    # Create an iterator function for TSS sites
    tss_iter = methodical:::.tss_iterator(meth_values_chunk, tss_region_indices_list, transcript_values, tss_for_chunk, 
      cor_method, add_distance_to_region, results_dir)
    
    # Calculate correlations for all TSS in chunk. 
    suppressWarnings(chunk_correlations <- BiocParallel::bpiterate(ITER = tss_iter, FUN = .tss_correlations, BPPARAM = BPPARAM))
    
    # Add names of transcript to list and return
    chunk_correlations <- setNames(chunk_correlations, names(tss_for_chunk))
      
    chunk_correlations
  }

  # Run garbage collection one final time and return result
  invisible(gc())
  all_correlations <- unlist(all_correlations, recursive = FALSE)
  return(all_correlations)
  
}
