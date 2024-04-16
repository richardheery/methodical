#' Calculate the correlation values between the methylation of genomic regions and the expression of associated transcripts
#'
#' @param meth_rse A RangedSummarizedExperiment with methylation values for CpG sites which will be used to calculate methylation values for genomic_regions. 
#' There must be at least 3 samples in common between meth_rse and transcript_expression_table.
#' @param assay_number The assay from meth_rse to extract values from. Default is the first assay. 
#' @param transcript_expression_table A table with the expression values for different transcripts in different samples. 
#' Row names should give be the transcript name and column names should be the name of samples. 
#' @param samples_subset Optional sample names used to subset meth_rse and transcript_expression_table. 
#' Provided samples must be found in both meth_rse and transcript_expression_table.
#' Default is to use all samples in meth_rse and transcript_expression_table.
#' @param genomic_regions A GRanges object. 
#' @param genomic_region_names Names for genomic_regions. If not provided, attempts to use names(genomic_regions). 
#' @param genomic_region_transcripts Names of transcripts associated with each region in genomic_regions. 
#' If not provided, attempts to use genomic_regions$transcript_id. All transcripts must be present in transcript_expression_table.
#' @param genomic_region_methylation Optional preprovided table with methylation values for genomic_regions 
#' such as created using summarizeRegionMethylation(). Table will be created if it is not provided which will increase running time.
#' Row names should match genomic_region_names and column names should match those of transcript_expression_table 
#' @param cor_method A character string indicating which correlation coefficient is to be computed. 
#' One of either "pearson" or "spearman" or their abbreviations. 
#' @param p_adjust_method Method used to adjust p-values. Same as the methods from p.adjust.methods. Default is Benjamini-Hochberg.
#' @param region_methylation_summary_function A function that summarizes column values. Default is colMeans.
#' @param BPPARAM A BiocParallelParam object for parallel processing. Defaults to `BiocParallel::bpparam()`. 
#' @param ... Additional arguments to be passed to summary_function. 
#' @return A data.frame with the correlation values between the methylation of genomic regions and expression of transcripts associated with them
#' @export
#' @examples 
#' 
#' # Load TUBB6 TMRs, RangedSummarizedExperiment with methylation values for CpGs around TUBB6 TSS and TUBB6 transcript counts
#' data(tubb6_tmrs, package = "methodical")
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#' data(tubb6_transcript_counts, package = "methodical")
#' 
#' # Calculate correlation values between TMRs identified for TUBB6 and transcript expression
#' tubb6_tmrs_transcript_cors <- methodical::calculateRegionMethylationTranscriptCors(
#'   meth_rse = tubb6_meth_rse, transcript_expression_table = tubb6_transcript_counts,
#'   genomic_regions = tubb6_tmrs, genomic_region_names = tubb6_tmrs$tmr_name)
#' tubb6_tmrs_transcript_cors
#'  
calculateRegionMethylationTranscriptCors <- function(meth_rse, assay_number = 1, transcript_expression_table, samples_subset = NULL, 
  genomic_regions, genomic_region_names = NULL, genomic_region_transcripts = NULL, genomic_region_methylation = NULL,
  cor_method = "pearson", p_adjust_method = "BH", region_methylation_summary_function = colMeans, BPPARAM = BiocParallel::bpparam(), ...){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), is(assay_number, "numeric"),
    is(transcript_expression_table, "data.frame") | is(transcript_expression_table, "matrix"),
    is(samples_subset, "character") | is.null(samples_subset), 
    is(genomic_regions, "GRanges"), is(genomic_region_names, "character") | is.null(genomic_region_names),
    is(genomic_region_transcripts, "character") | is.null(genomic_region_transcripts), 
    is(genomic_region_methylation, "data.frame") | is(genomic_region_methylation, "matrix") | is.null(genomic_region_methylation),
    is(p_adjust_method, "character") & p_adjust_method %in% p.adjust.methods,
    is(region_methylation_summary_function, "function"), is(BPPARAM, "BiocParallelParam"))
    match.arg(cor_method, choices = c("pearson", "spearman"))
  
  # Check that assay_number is not greater than the number of assays in meth_rse
  if(assay_number > length(SummarizedExperiment::assays(meth_rse))){
    stop(paste0("Provided value for assay_number (", assay_number, ") is greater than the number of assays in meth_rse"))
  }
  
  # Check that samples_subset are in meth_rse and transcript_expression_table
  if(!is.null(samples_subset)){
    if(any(!samples_subset %in% colnames(meth_rse))){
      stop("Some samples in samples_subset are not in meth_rse")
    } else {meth_rse <- meth_rse[, samples_subset]}
    if(any(!samples_subset %in% colnames(transcript_expression_table))){
      stop("Some samples in samples_subset are not in transcript_expression_table")
    } else {transcript_expression_table <- transcript_expression_table[, samples_subset, drop = FALSE]}
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
  
  # If genomic_region_names not provided, set to names(genomic_regions)
  if(is.null(genomic_region_names)){
    genomic_region_names <- names(genomic_regions)
  }
  
  # Add names to genomic_regions if they are not already present and also check that no names are duplicated. 
  if(is.null(genomic_region_names)){
    message("No names for provided regions so naming them region_1, region_2, etc.")
    genomic_region_names <- paste0("region_", seq_along(genomic_regions))
    names(genomic_regions) <- genomic_region_names
  } else {
    if(length(genomic_region_names) != length(genomic_regions)){
      stop("genomic_region_names must be the same length as genomic_regions")
    } 
    if(anyDuplicated(genomic_region_names)){
      stop("genomic_region_names cannot contain duplicates")
    } else {
      names(genomic_regions) <- genomic_region_names
    }
  }
  
  # If genomic_regions$transcript_id is NULL, set to genomic_region_transcripts
  if(!is.null(genomic_region_transcripts)){
    if(length(genomic_region_transcripts) != length(genomic_regions)){
      stop("If provided, length of genomic_region_transcripts must equal length of genomic_regions")
    } else {genomic_regions$transcript_id <- genomic_region_transcripts}
  } else if(is.null(genomic_regions$transcript_id)){
     stop("genomic_regions$transcript_id and genomic_region_transcripts cannot both be NULL")
  }
    
  # Add names and transcript IDs to genomic_regions
  names(genomic_regions) <- genomic_region_names
  
  # Identify transcripts in common between tss_gr and transcript_expression_table.
  # Throw an error if there are no common transcripts and subset tss_gr and transcript_expression_table
  common_transcripts <- intersect(genomic_regions$transcript_id, row.names(transcript_expression_table))
  if(length(common_transcripts) == 0){
    stop("There are no common transcripts/genes between genomic_region_transcripts and row.names(transcript_expression_table)")
  } else {
    message(paste("There are", length(common_transcripts), "genes/transcripts in common between genomic_region_transcripts and transcript_expression_table"))
  }
  transcript_expression_table <- transcript_expression_table[common_transcripts, , drop = FALSE]
  genomic_regions <- genomic_regions[genomic_regions$transcript_id %in% common_transcripts]
  
  # Check that there are no duplicates for genomic_region_names and that it has the same length as genomic_regions
  # if(length(genomic_region_names) != length(genomic_regions)){stop("genomic_region_names must have the same length as genomic_regions")}
  if(anyDuplicated(genomic_region_names)){stop("genomic_region_names cannot contain duplicates")}
  
  # If genomic_region_methylation provided, check that it has the correct row and column names. 
  if(!is.null(genomic_region_methylation)){
    if(any(row.names(genomic_region_methylation) != genomic_region_names)){
      stop("Row names of genomic_region_methylation do not match genomic_region_names")
    }
    if(any(colnames(genomic_region_methylation) != colnames(transcript_expression_table))){
      if(is.null(samples_subset)){
        stop("Column names of genomic_region_methylation do not match those of transcript table")
      } else {
        stop("Column names of genomic_region_methylation do not match samples_subset")
      }
    }
  }
  
  # Print correlation method being used
  message(paste("Using", match.arg(cor_method, c("pearson", "kendall", "spearman")), "correlation method"))
  
  # Create genomic_region_methylation if it is not provided
  if(is.null(genomic_region_methylation)){
  
    # Summarize methylation values for genomic_regions
    genomic_region_methylation <- summarizeRegionMethylation(
      meth_rse = meth_rse, assay_number = assay_number, genomic_regions = genomic_regions, 
      genomic_region_names = names(genomic_regions), 
      summary_function = region_methylation_summary_function, BPPARAM = BPPARAM, ...
    )
    
  }
  
  # Calculate correlations between region methylation values and transcript expression
  message("Calculating transcript correlations")
  
  # Create a data.frame matching genomic_region_names and transcript_names
  feature_matches_df <- 
    data.frame(genomic_region_names = names(genomic_regions), transcript_id = genomic_regions$transcript_id)
  
  # Create a list of the transcripts associated with each genomic region
  feature_matches <- split(feature_matches_df$genomic_region_names, feature_matches_df$transcript_id)
  
  # For each transcript, calculate the correlation between its expression and the methylation of regions associated with it
  methylation_transcript_correlations <- suppressMessages(BiocParallel::bplapply(X = names(feature_matches), 
    FUN = function(X) {rapidCorTest(
      table1 = t(genomic_region_methylation[feature_matches[[X]], ]),
      table2 = setNames(data.frame(unlist(transcript_expression_table[X, ])), X),
      table1_name = "genomic_region_name", table2_name = "transcript_name"
      )}, BPPARAM = BPPARAM))
    
  # Combine results into a single table
  methylation_transcript_correlations <- dplyr::bind_rows(methylation_transcript_correlations)
  
  # Adjust p-values
  methylation_transcript_correlations$q_val <- p.adjust(methylation_transcript_correlations$p_val, method = p_adjust_method)
  
  # Put methylation_transcript_correlations in order of genomic_region_names
  methylation_transcript_correlations <- methylation_transcript_correlations[match(
    names(genomic_regions), methylation_transcript_correlations$genomic_region_name), ]
  row.names(methylation_transcript_correlations) <- NULL
  
  # Return results
  return(methylation_transcript_correlations)
  
}
