#' Calculate correlation between methylation of CpGs surrounding genomic regions and transcription of associated features
#' 
#' AT THE MOMENT, IT RELIES ON transcript_expression_table AND genomic_regions HAVING MATCHING TRANSCRIPT IDS. MAYBE CHANGE THIS. 
#' MAYBE ADD OPTION NOT TO ADD DISTANCE OR SIGNIFICANCE
#' ADD PARAMETER WITH DESIDRED CHROMOSOME ORDER MAYBE
#' ADD FUNCTION TO COMBINE RESULTS INTO A SINGLE TABLE
#' Sort by distance? Maybe not if it increases running time a lot
#' index using numbers rather than characters
#'
#' @param methrix_object A methrix object. Names of samples must match those in transcript_expression_table unless subset_samples provided.
#' @param transcript_expression_table A table with the expression values for transcripts, where rows correspond to transcripts and columns to samples. 
#' There should be one row for each region in genomic_regions. Names of samples must match those in methrix_object unless subset_samples provided.
#' @param subset_samples Sample names to subset methrix object and transcript_expression_table. Provided samples must be found in both methrix_object and transcript_expression_table.
#' Default is to use all samples in methrix_object and transcript_expression_table. 
#' @param genomic_regions Genomic regions to calculate mean CpG methylation for. There should be one region for each row in transcript_expression_table. 
#' Must have a metadata column called transcript_id that is identical to the row.names of transcript_expression_table.
#' @param expand_upstream Number of bases to add upstream of each region in genomic regions. Must be numeric vector of length 1 or equal to the length of genomic_regions. Default is 0.
#' @param expand_downstream Number of bases to add downstream of each region in genomic regions. Must be numeric vector of length 1 or equal to the length of genomic_regions. Default is 0.
#' @param cor_method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor(). Default is "pearson".
#' @param n_permutations Number of permutations for permuation test. Default is 10,000.
#' @param add_distance_to_region A logical value indicating whether to add the distance of CpG sites to the center of the associated genomic region. Default value is TRUE. 
#' Setting to FALSE will roughly half the total running time. 
#' @param ncores Number of cores to use. Default is 1.
#' @param random_seed Integer to use as a random seed. Default is 123.
#' @return A list of data.frames with the correlation of CpG sites surrounding a specified genomic region with a given feature, p-values and adjusted q-values for the correlations.
#' Distance of the CpG upstream or downstream to the center of the region is also provided.
#' @export
transcription_cpg_correlations = function(methrix_object, transcript_expression_table, subset_samples = NULL, genomic_regions, expand_upstream = 0, 
  expand_downstream = 0, cor_method = "pearson", n_permutations = 10000, add_distance_to_region = T, ncores = 1, random_seed = 123){
  
  # Check that cor_method is wither pearson or spearman
  match.arg(cor_method, choices = c("pearson", "spearman"), several.ok = F)
  
  # Check that GRanges object and transcript_expression_table have the correct dimensions
  if(length(genomic_regions) != nrow(transcript_expression_table)){stop("Length of genomic_regions should equal the number of rows in transcript_expression_table")}
  
  # Check that transcript_id of genomic_regions corresponds to row.names of transcript_expression_table
  if(!all(genomic_regions$transcript_id == row.names(transcript_expression_table))){stop("genomic_regions$transcript_id must match row.names(transcript_expression_table)")}
  
  # Check that subset_samples are in methrix object and transcript_expression_table
  if(!is.null(subset_samples)){
    if(any(!subset_samples %in% row.names(colData(methrix_object)))){
      stop("Some subset_samples are not in methrix_object")
    } else {methrix_object = suppressMessages(subset_methrix(methrix_object, samples = subset_samples))}
    if(any(!subset_samples %in% names(transcript_expression_table))){
      stop("Some subset_samples are not in transcript_expression_table")
    } else {transcript_expression_table = dplyr::select(transcript_expression_table, dplyr::all_of(subset_samples))}
  }
  
  # Check that names of methrix_object and transcript_expression_table match
  if(!all(row.names(colData(methrix_object)) == names(transcript_expression_table))){stop("Sample names in methrix_object and transcript_expression_table do not match")}
  
  # Get sequences common to both genomic_regions and transcript_expression_table. !! MAYBE GET INTERSECT OF SEQUENCES FOR TRANSCRIPTS, METHRIX AND GENOMIC_RANGES
  common_sequences = gtools::mixedsort(intersect(unique(rowData(methrix_object)$chr), unique(seqnames(genomic_regions))))
  
  # Create cluster if ncores greater than 1
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, ncores)
  `%do%` = foreach::`%do%`
  `%dopar%` = foreach::`%dopar%`
  `%dorng%` = doRNG::`%dorng%`
  
  # Find the center of the original genomic regions
  genomic_regions_centers = resize(genomic_regions, width = 1, fix = "center")
  
  # Expand genomic_regions and subset methrix object for these regions
  genomic_regions_expanded = methodical::expand_granges(genomic_regions, expand_upstream, expand_downstream)
  
  # Get methylation values associated with chrom. Takes 25 seconds for chromosome setup. Takes 8 minutes for chr1 with 3 cores. 
  chrom_results = foreach::foreach(chrom = common_sequences) %do% {
    system2("echo", paste("Starting regions on", chrom))
    genomic_regions_expanded_chrom =  genomic_regions_expanded[seqnames(genomic_regions_expanded) == chrom]
    transcript_expression_table_chrom = transcript_expression_table[genomic_regions_expanded_chrom$transcript_id, ]
    methrix_chrom = suppressMessages(subset_methrix(methrix_object, regions = reduce(genomic_regions_expanded_chrom, ignore.strand = T)))
    
    # Extract methylation values for the regions located on the chromosome and then run the garbage collection
    cpg_values_chrom = suppressMessages(methodical::extract_cpg_values_from_methrix(methrix = methrix_chrom, genomic_regions = genomic_regions_expanded_chrom))
    invisible(gc())
    
    # Find names of CpG sites overlapping each region
    cpg_overlaps = data.frame(findOverlaps(GRanges(row.names(cpg_values_chrom)), genomic_regions_expanded_chrom))
    cpg_overlaps$transcript_id = genomic_regions_expanded_chrom$transcript_id[cpg_overlaps$subjectHits]
    cpg_overlaps = split(cpg_overlaps$queryHits, cpg_overlaps$transcript_id)
    
    # Calculate correlations. Takes 5 minutes for chromosome 1 with 3 cores
    system.time({all_correlations = foreach::foreach(transcript = names(cpg_overlaps), .packages = c("methodical"), .options.RNG = random_seed) %dorng% {
      
      # Create table where columns are CpGs and rows are samples
      current_cpg_table = t(cpg_values_chrom[cpg_overlaps[[transcript]], ])
      
      # Create a table with one column for the transcript expression values
      current_transcript_table = t(transcript_expression_table_chrom[transcript, ])
      
      # If cor_method is Spearman, change CpG values and transcript values to ranks
      if(cor_method == "spearman"){
        current_cpg_table = apply(current_cpg_table, 2, function(x) rank(x, na.last = "keep"))
        current_transcript_table = apply(current_transcript_table, 2, function(x) rank(x, na.last = "keep"))
      }
      
      # Create a table with the transcript values shuffled 
      shuffled_transcript_table = matrix(unlist(lapply(1:n_permutations, function(x) sample(current_transcript_table))), ncol = n_permutations)
      
      # Calculate actual correlations between CpG values and the transcript
      transcript_cpg_correlations = setNames(data.frame(WGCNA::cor(current_cpg_table, current_transcript_table, method = "pearson", use = "p")), "correlation")
      
      # Calculate random correlations 
      random_cors = WGCNA::cor(current_cpg_table, shuffled_transcript_table, method = "pearson", use = "p")
      
      # HAVING >= VERY IMPORTANT!
      transcript_cpg_correlations$p_value = (rowSums(abs(random_cors) >= (abs(transcript_cpg_correlations[["correlation"]]) -0.001)) + 1)/(n_permutations + 1)
      
      # Add CpG distance to region if specified
      if(add_distance_to_region){transcript_cpg_correlations$distance_to_region = methodical::cpg_distances_to_region(
        query_region = genomic_regions_centers[genomic_regions$transcript_id == transcript], cpg_names = row.names(transcript_cpg_correlations))}
      
      transcript_cpg_correlations
    }})
    
    # Add names of transcript to list and return
    all_correlations = setNames(all_correlations, names(cpg_overlaps))
    all_correlations
  }
  
  # Stop the cluster if present
  parallel::stopCluster(cl)
  
  # Run garbage collection one final time and return result
  invisible(gc())
  return(unlist(chrom_results, recursive = F))
}

#' Created a data.frame combining the results from the output of transcription_cpg_correlations()
#' 
#' @param transcription_cpg_correlations_list A list created by transcription_cpg_correlations()
#' @param correct_p_values Logical value indicating whether p-values in the final table should be corrected, Default is TRUE. 
#' @param p_adjust_method Method used to adjust p-values if correct_p_values is specified. Same as the methods from p.adjust.methods. Default is Benjamini-Hochberg. 
#' @export
combine_transcription_cpg_correlations = function(transcription_cpg_correlations_list, correct_p_values = T, p_adjust_method = "BH"){
  
  # Combine results into a table and turn rownames into a column called cpg_name
  transcription_cpg_correlations_df = tibble::rownames_to_column(dplyr::bind_rows(transcription_cpg_correlations_list, .id = "transcript_id"), "cpg_name")
  
  # Remove ...1 etc added to duplicate cpg names. 
  transcription_cpg_correlations_df$cpg_name = gsub("\\..*", "", transcription_cpg_correlations_df$cpg_name)
  
  # Correct p-values if specified
  if(correct_p_values){transcription_cpg_correlations_df$q_value = p.adjust(transcription_cpg_correlations_df$p_value, method = p_adjust_method)}
  
  return(transcription_cpg_correlations_df)
}