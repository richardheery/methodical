#' Summarize methylation values for CpG regions surrounding genomic regions. MAYBE LIMIT THE FUNCTIONS THAT CAN BE USED WITH SUMMARY_FUNCTION. 
#' MAYBE CAN SPEED UP SAMPLE_GROUP2 DIFFERENCE.
#'
#' @param methrix_object A methrix object. Names of samples must match those in transcript_expression_table unless samples_group1 is provided.
#' @param samples_group1 Sample names to subset methrix object and transcript_expression_table. Provided samples must be found in both methrix_object and transcript_expression_table.
#' Default is to use all samples in methrix_object and transcript_expression_table. 
#' @param samples_group2 Sample names to subset methrix object and transcript_expression_table. Provided samples must be found in both methrix_object and transcript_expression_table.
#' Default is to use all samples in methrix_object and transcript_expression_table. 
#' @param genomic_regions Genomic regions to calculate mean CpG methylation for. There should be one region for each row in transcript_expression_table. 
#' Must have a metadata column called transcript_id that is identical to the row.names of transcript_expression_table.
#' @param expand_upstream Number of bases to add upstream of each region in genomic regions. Must be numeric vector of length 1 or equal to the length of genomic_regions. Default is 0.
#' @param expand_downstream Number of bases to add downstream of each region in genomic regions. Must be numeric vector of length 1 or equal to the length of genomic_regions. Default is 0.
#' @param summary_function A summary function to apply to the methylation values of each CpG site. Default is rowMeans().
#' @param summary_column_name The name to give the column with the output of the summary function. Default is "methylation_summary".
#' @param add_distance_to_region A logical value indicating whether to add the distance of CpG sites to the center of the associated genomic region. Default value is TRUE. 
#' Setting to FALSE will roughly half the total running time. 
#' @param ncores Number of cores to use. Default is 1.
#' @return A list of data.frames with a summary of the methylation values of CpG sites surrounding a specified genomic region.
#' Distance of the CpG upstream or downstream to the center of the region is also provided.
#' @export
cpg_methylation_summary = function(methrix_object, samples_group1 = NULL, samples_group2 = NULL, genomic_regions, expand_upstream = 0, 
  expand_downstream = 0, summary_function = rowMeans, summary_column_name = "methylation_summary", add_distance_to_region = T, ncores = 1){
  
  # Check that if samples_group2 is supplied, so is samples_group1
  if(!is.null(samples_group2)){
    if(is.null(samples_group1)){stop("If samples_group2 is provided so must samples_group1")}
  }
  
  # Check that there are no samples in common between samples_group1 and samples_group2
  if(any(samples_group1 %in% samples_group2)){stop("There cannot be samples in common between samples_group1 and samples_group2")}
  
  # Check that samples_group1 and samples_group2 are in methrix object. If samples_group1 not provided, it is set to be all samples in methrix_object.
  if(!is.null(samples_group1)){
    if(any(!c(samples_group1, samples_group2) %in% row.names(colData(methrix_object)))){
      stop("Some provided samples are not in methrix_object")
    } else {methrix_object = suppressMessages(subset_methrix(methrix_object, samples = c(samples_group1, samples_group2)))}
  } else {samples_group1 = row.names(colData(methrix_object))}
  
  # Get sequences common to both genomic_regions and methrix_object. 
  common_sequences = gtools::mixedsort(intersect(unique(rowData(methrix_object)$chr), unique(seqnames(genomic_regions))))
  
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
  
  # Find the center of the original genomic regions
  genomic_regions_centers = resize(genomic_regions, width = 1, fix = "center")
  
  # Expand genomic_regions and subset methrix object for these regions
  genomic_regions_expanded = methodical::expand_granges(genomic_regions, expand_upstream, expand_downstream)
  
  # Get methylation values associated with chrom. Takes 25 seconds for chromosome setup. Takes 6 minutes for chr1 with 1 cores. Takes 3 minutes with 3 cores. 
  chrom_results = foreach::foreach(chrom = common_sequences) %do% {
    system2("echo", paste("Starting regions on", chrom))
    genomic_regions_expanded_chrom =  genomic_regions_expanded[seqnames(genomic_regions_expanded) == chrom]
    methrix_chrom = suppressMessages(subset_methrix(methrix_object, regions = reduce(genomic_regions_expanded_chrom, ignore.strand = T)))
    
    # Extract methylation values for the regions located on the chromosome and then run the garbage collection
    cpg_values_chrom = suppressMessages(as.matrix(methodical::extract_cpg_values_from_methrix(methrix = methrix_chrom, genomic_regions = genomic_regions_expanded_chrom)))
    invisible(gc())
    
    # Find names of CpG sites overlapping each region
    cpg_overlaps = data.frame(findOverlaps(GRanges(row.names(cpg_values_chrom)), genomic_regions_expanded_chrom))
    cpg_overlaps$transcript_id = genomic_regions_expanded_chrom$transcript_id[cpg_overlaps$subjectHits]
    cpg_overlaps = split(cpg_overlaps$queryHits, cpg_overlaps$transcript_id)
    
    # Calculate summaries. Took 6 minutes with 1 core. 
    all_summaries = foreach::foreach(transcript = names(cpg_overlaps), .packages = c("methodical")) %dopar% {
      cpg_summary = setNames(data.frame(summary_function(cpg_values_chrom[cpg_overlaps[[transcript]], samples_group1], na.rm = T)), summary_column_name)
      row.names(cpg_summary) = row.names(cpg_values_chrom[cpg_overlaps[[transcript]], samples_group1])
      
      if(!is.null(samples_group2)){cpg_summary = cpg_summary - setNames(data.frame(summary_function(cpg_values_chrom[cpg_overlaps[[transcript]], samples_group2], na.rm = T)), summary_column_name)}
      
      # Add CpG distance to region if specified
      if(add_distance_to_region){cpg_summary$distance_to_region = methodical::cpg_distances_to_region(
        query_region = genomic_regions_centers[genomic_regions$transcript_id == transcript], cpg_names = row.names(cpg_summary))}
      
      cpg_summary
    }
    
    # Add names of transcript to list and return
    all_summaries = setNames(all_summaries, names(cpg_overlaps))
    all_summaries
  }
  
  # Stop cluster, un garbage collection one final time and return result
  if(ncores > 1){stopCluster(cl)}
  invisible(gc())
  return(unlist(chrom_results, recursive = F))
}