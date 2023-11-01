#' Summarize methylation values for regions in a chunk
#'
#' @param chunk_regions Chunk with genomic regions of interest. 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param assay_number The assay from meth_rse to extract values from.
#' @param summary_function summary_function A function that summarizes column values.
#' @param na.rm A logical value indicating whether to remove NA values when calculating summaries.
#' @param ... Additional arguments to be passed to summary_function. 
#' @return A function which returns a list with the 
.summarize_chunk_methylation = function(chunk_regions, meth_rse, assay_number, summary_function, na.rm, ...){
  
  # Subset meth_rse_for_chunk for regions overlapping chunk_regions
  meth_rse_for_chunk <- subsetByOverlaps(meth_rse, chunk_regions)
  invisible(gc()) 
  
  # Find the overlaps of chunk_regions and meth_rse_for_chunk
  overlaps_df <- data.frame(findOverlaps(chunk_regions, meth_rse_for_chunk))
  
  # Add region names to overlaps_df
  overlaps_df$genomic_region_name <- names(chunk_regions)[overlaps_df$queryHits]
  
  # Create a list matching region names to rows of meth_rse_for_chunk
  region_names_to_rows_list <- split(overlaps_df$subjectHits, overlaps_df$genomic_region_name)
  
  # Read all values from specified assay of meth_rse_for_chunk into memory and run the garbage collection
  meth_values <- as.matrix(SummarizedExperiment::assay(meth_rse_for_chunk, i = assay_number))
  gc()
  
  # Summarize methylation values 
  meth_summary <- lapply(region_names_to_rows_list, function(x) 
    summary_function(meth_values[x, , drop = FALSE], na.rm = na.rm, ...))
  rm(meth_values); gc()
  
  # Combine meth_summary into a single table
  meth_summary <- do.call("rbind", meth_summary)
  
  # Convert meth_summary to a data.frame
  meth_summary <- data.frame(meth_summary)
  gc()
  meth_summary

}

#' Summarize methylation of genomic regions
#' 
#' @param meth_rse A RangedSummarizedExperiment with methylation values.
#' @param assay_number The assay from meth_rse to extract values from. Default is the first assay. 
#' @param genomic_regions GRanges object with regions to summarize methylation values for. 
#' @param keep_metadata_cols A logical value indicating whether to add the metadata columns of genomic_regions to the output. Default is FALSE.
#' @param genomic_region_names A vector of names to give genomic_regions in the output table. There cannot be any duplicated names. 
#' Default is to attempt to use `names(genomic_regions)` if they are present or to name them region_1, region_2, etc otherwise.
#' @param max_sites_per_chunk The approximate maximum number of methylation sites to try to load into memory at once. 
#' The actual number loaded may vary depending on the number of methylation sites overlapping each region, 
#' but so long as the size of any individual regions is not enormous (>= several MB), it should vary only very slightly. 
#' Some experimentation may be needed to choose an optimal value as low values will result in increased running time, 
#' while high values will result in a large memory footprint without much improvement in running time. 
#' Default is floor(62500000/ncol(meth_rse)), resulting in each chunk requiring approximately 500 MB of RAM. 
#' @param summary_function A function that summarizes column values. Default is base::colMeans. 
#' @param na.rm A logical value indicating whether to remove NA values when calculating summaries. Default value is TRUE. 
#' @param BPPARAM A BiocParallelParam object. Defaults to `BiocParallel::bpparam()`. 
#' @param ... Additional arguments to be passed to summary_function. 
#' @return A data.table with the summary of methylation of each region in genomic_regions for each sample.
#' @export
#' @examples 
#' 
#' # Load sample RangedSummarizedExperiment with CpG methylation data
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#' 
#' # Create a sample GRanges
#' test_gr <- GRanges(c("chr18:12303400-12303500", "chr18:12303600-12303750", "chr18:12304000-12306000"))
#' names(test_gr) <- paste("region", 1:3, sep = "_")
#' 
#' # Calculate mean methylation values for chr1 CpG islands in meth_h5 
#' test_gr_methylation <- methodical::summarizeRegionMethylation(tubb6_meth_rse, genomic_regions = test_gr,
#'   genomic_region_names = names(test_gr))
#' 
summarizeRegionMethylation <- function(meth_rse, assay_number = 1, genomic_regions, genomic_region_names = NULL, 
  keep_metadata_cols = FALSE, max_sites_per_chunk = NULL, summary_function = base::colMeans, na.rm = TRUE, BPPARAM = BiocParallel::bpparam(), ...){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), is(assay_number, "numeric"),
    is(genomic_regions, "GRanges"), is(genomic_region_names, "character") | is.null(genomic_region_names),
    is(keep_metadata_cols, "logical"), 
    (is(max_sites_per_chunk, "numeric") & max_sites_per_chunk >= 1) | is.null(max_sites_per_chunk),
    is(summary_function, "function"), is(na.rm, "logical"), 
    is(BPPARAM, "BiocParallelParam"))
  
  # If genomic_region_names if NULL, attempt to use names of genomic_regions
  genomic_region_names = names(genomic_regions)
    
  # Add names to genomic_regions if they are not already present and also check that no names are duplicated. 
  if(is.null(genomic_region_names)){
    message("No names for provided regions so naming them region_1, region_2, etc.")
    genomic_region_names <- paste0("region_", 1:length(genomic_regions))
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
  
  # Check that there are some seqlevels from genomic_regions present in meth_rse and 
  # give an error if not and a warning if some are missing 
  if(all(!seqlevels(genomic_regions) %in% seqlevels(meth_rse))){
    stop("There are no seqlevels in common between genomic_regions and meth_rse")
  } 
  if(any(!seqlevels(genomic_regions) %in% seqlevels(meth_rse))){
    message("There are some seqlevels from genomic_regions missing from meth_rse")
  } 
  
  # Split genomic regions into chunks based on the number of methylation sites that they cover
  genomic_region_bins <- .chunk_regions(meth_rse = meth_rse, genomic_regions = genomic_regions, 
    max_sites_per_chunk = max_sites_per_chunk, ncores = BiocParallel::bpnworkers(BPPARAM))
  
  # For each sequence get methylation of the associated regions
  message("Summarizing genomic region methylation")
  region_methylation <- BiocParallel::bpmapply(FUN = .summarize_chunk_methylation, 
    chunk_regions = genomic_region_bins, MoreArgs = list(meth_rse = meth_rse, 
      assay_number = assay_number, summary_function = summary_function, na.rm = na.rm, ...), 
    BPPARAM = BPPARAM, SIMPLIFY = FALSE)

  # Combine data.frames for each chunk
  region_methylation <- dplyr::bind_rows(region_methylation)
  gc()
  
  # Turn rownames into a column and convert the result to a data.table
  region_methylation <- data.table::data.table(tibble::rownames_to_column(region_methylation, "region_name"))
  
  # Create a data.table with the genomic_region_names in the correct order
  genomic_region_names_df <- data.table::data.table(region_name = genomic_region_names)
  
  # Put rows in same order as regions in genomic_regions. Adds rows with NA values for regions which didn't overlap any methylation sites. 
  region_methylation <- data.table::merge.data.table(genomic_region_names_df, region_methylation, 
    by = "region_name", all.x = TRUE, sort = FALSE)

  # Add metadata from genomic_regions if specified
  if(keep_metadata_cols){
    region_methylation <- cbind(region_methylation, data.frame(mcols(genomic_regions)))
  }
  
  # Run the garbage collection and return region_methylation
  gc()
  return(region_methylation)

}
