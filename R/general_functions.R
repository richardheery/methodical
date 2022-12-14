# NEED TO SPEED UP sample_cpgs

#' Add dummy coverage values columns to bedGraphs so that they can be used with methrix::read_bedgraphs()
#' 
#' Each CpG site is given a coverage value of 50. 
#'
#' @param bedgraphs A list of paths to bedGraph files 
#' @param output_directory A directory in which to save the modified bedGraph files. Directory is created if it doesn't exist. Default name is "dummy_coverage_bedgraphs". 
#' Can overwrite the original bedGraphs if output_directory is the directory in which they are located. 
#' @param ncores Number of cores to use. Default is 1.
#' @tmpdir Directory to store temporary files in. Default is tempdir().
#' @return Invisibly returns TRUE if all bedGraphs converted without errors. 
#' @export
add_dummy_coverage_to_bedgraphs = function(bedgraphs, output_directory = "dummy_coverage_bedgraphs", ncores = 1, tmpdir = tempdir()){
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_directory)){dir.create(output_directory)}
  
  # Create temporary directory 
  temporary_directory = tempfile(pattern = "dummy_coverage_bedgraphs_temp_", tmpdir = tmpdir)
  dir.create(temporary_directory)
  
  # Set up cluster if ncores > 1
  if(ncores > 1){
    cl = parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl, ncores)
    `%dopar%` = foreach::`%dopar%`
  } else {`%dopar%` = foreach::`%do%`}
  
  # Loop through each bedgraph and add dummy coverage using awk
  suppressWarnings({foreach::foreach(bg = bedgraphs) %dopar% {
    
    # Make a temporary file
    temp_file = tempfile(tmpdir = temporary_directory)
    
    # Detect if file is compressed using gzip
    if(tools::file_ext(bg) == "gz"){
      compress_command = "gzip -c"
    } else {compress_command = "cat"}
    
    # Add dummy coverage value of 50 and write to temporary file, move it to output_directory and then make sure the temporary file is deleted. 
    system(sprintf("zcat -f %s | awk -v OFS='\t' '{print $0, 50}' |  %s > %s", bg, compress_command, temp_file))
    file.copy(temp_file, paste(output_directory, basename(bg), sep = "/"), overwrite = T)
    unlink(temp_file)
  }})
  
  # Stop cluster
  if(ncores > 1){
    parallel::stopCluster(cl)
  }
  
  unlink(temporary_directory, recursive = T)
  
  return(invisible(TRUE))
}

#' Expand GRanges. 
#'
#' Expand ranges in a GRanges object upstream and downstream by specified numbers of bases, taking strand into account.
#' If GRanges do not have a strand, they are treated like they are on the "+" strand. 
#'
#' @param genomic_regions A GRanges object
#' @param upstream Number of bases to add upstream of each region in genomic_regions. Must be numeric vector of length 1 or else equal to the length of genomic_regions. Default is 0.
#' @param downstream Number of bases to add downstream of each region in genomic_regions. Must be numeric vector of length 1 or else equal to the length of genomic_regions. Default is 0.
#' @return None
#' @export
expand_granges = function(genomic_regions, upstream = 0, downstream = 0) {
  
  # Check that upstream and downstream are positive integers 
  if(!all(upstream == abs(floor(upstream))) | !length(upstream) %in% c(1, length(genomic_regions))){
    stop("upstream should be a vector of positive intergers of length 1 or the length of genomic_regions")}
  if(!all(downstream == abs(floor(downstream))) | !length(downstream) %in% c(1, length(genomic_regions))){
    stop("downstream should be a vector of positive intergers of length 1 or the length of genomic_regions")}
  
  # Check for each range if it's on the negative or positive strand
  strand_is_minus = as.character(GenomicRanges::strand(genomic_regions)) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  
  # Adjust ranges based on whether they are on the positive or negative strand
  GenomicRanges::start(genomic_regions)[on_plus] = GenomicRanges::start(genomic_regions)[on_plus] - upstream
  GenomicRanges::start(genomic_regions)[on_minus] = GenomicRanges::start(genomic_regions)[on_minus] - downstream
  GenomicRanges::end(genomic_regions)[on_plus] = GenomicRanges::end(genomic_regions)[on_plus] + downstream
  GenomicRanges::end(genomic_regions)[on_minus] = GenomicRanges::end(genomic_regions)[on_minus] + upstream
  
  # Remove any out-of-bounds regions
  genomic_regions = GenomicRanges::trim(genomic_regions)
  
  return(genomic_regions)
} 

# cpg_granges_from_methrix = function(methrix_object){
#   cpg_df = rowData(methrix_object)
#   cpg_granges = GenomicRanges::GRanges(paste0(cpg_df$chr, ":", cpg_df$start)))
#   gc()
#   return(cpg_granges)
# }

#' Get methylation values of CpG sites overlapping a GRanges object from a methrix object
#'
#' @param methrix A methrix object
#' @param genomic_regions Genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a GenomicRanges object. Default is to extract all CpG sites in the methrix.
#' @param samples Sample names to subset by. Default is to get values for all samples.
#' @return A data.frame with the methylation values of overlapping CpG sites with CpG names as first column
#' @export
extract_cpg_values_from_methrix = function(methrix, genomic_regions = NULL, samples = NULL){
  
  # Subset methrix for genomic_regions if provided
  if(!is.null(genomic_regions)){
    methrix = methrix::subset_methrix(methrix, GenomicRanges::reduce(genomic_regions, ignore.strand = T))
  }
  
  # Subset methrix by samples 
  methrix = methrix::subset_methrix(methrix, samples = samples)
  
  # Extract CpG values for specified regions adding location of CpG sites
  cpg_values = methodical:::get_matrix_updated(methrix, type = "M", add_loci = T)
  
  # Add CpG names as row names and remove the loci columns
  cpg_values = tibble::column_to_rownames(dplyr::select(dplyr::mutate(cpg_values, cpg_name = paste(chr, start, sep = ":")), -chr, -start, -strand), "cpg_name") 
    
  return(cpg_values)
}

#' Find distances between start of a region and CpG sites
#'
#' Find distances between a query region and CpG sites, giving a negative distance if CpGs are upstream of the query region and positive values if they are downstream.
#' Adds 1 to the distance found using CpGs GenomicRanges::distance() so that CpGs touching a query_region have a distance of 1 and 
#' CpGs located within the region are assigned a distance of 0 bases from the regions.
#'
#' @param query_region A GRanges object with a single range for one query region.
#' @param cpg_names Names of CpG sites in the form sequence:start where start is 1-based e.g. chr1:10469 for the first CpG site on chr1 from hg38.  
#' @return None
#' @export
cpg_distances_to_region = function(query_region, cpg_names){
  
  # Check that query_region contains only a single range
  if(length(query_region) != 1 | class(query_region) != "GRanges"){stop("query_region should be a GRanges object with only one range")}
  
  # Get GRanges for CpG names
  cpg_granges = GenomicRanges::GRanges(cpg_names)

  # Get distances between query_region and cpg_granges. 
  distance_to_regions = start(cpg_granges) - start(query_region) 
  
  # Check if the cpg_ranges are upstream or downstream of the query_region 
  if(as.character(strand(query_region)) == "-"){distance_to_regions = distance_to_regions * -1}
  
  return(distance_to_regions)
}

#' Liftover CpGs sites from one genome build to another. 
#'
#' @param cpg_values A data.frame with methylation values for CpG values in different samples such as returned by extract_cpg_values_from_methrix().
#' @param liftover_chain A liftover Chain object specific to the genome build of the CpG sites in cpg_values and the desired target genome build. 
#' @return A data.frame with CpG methylation values for those CpGs which could be uniquely lifted over to the genome build indicated by liftover_chain. 
#' CpG sites which map to multiple sites in the target genome build or which are mapped to the same site as another CpG site are excluded. 
#' @export
liftover_cpgs = function(cpg_values, liftover_chain){
  
  # Create a GRanges object for the original CpG sites
  original_cpgs_gr = GenomicRanges::GRanges(row.names(cpg_values), original_position = row.names(cpg_values))
  
  # Perform the liftover using the chain file
  liftover_result = unlist(rtracklayer::liftOver(original_cpgs_gr, liftover_chain))
  liftover_result$liftover_position = paste(seqnames(liftover_result), GenomicRanges::start(liftover_result), sep = ":")
  
  # Create a data.frame with original CpG sites and their new positions after liftover 
  liftover_cpg_df = data.frame(original_position = liftover_result$original_position, liftover_position = liftover_result$liftover_position)
  
  # Remove liftovers which involve multiple-mappings in either direction
  liftover_cpg_df = dplyr::filter(liftover_cpg_df, !original_position %in% names(which(table(liftover_cpg_df$original_position) > 1)), 
    !liftover_position %in% names(which(table(liftover_cpg_df$liftover_position) > 1)))
  
  # Subset cpg_values for those which could be lifted over unambiguously and rename CpG sites according to their new position
  liftover_cpg_values = cpg_values[liftover_cpg_df$original_position, ]
  row.names(liftover_cpg_values) = liftover_cpg_df$liftover_position
  
  return(liftover_cpg_values)
}

#' Calculate correlation between methylation of CpG sites and a feature 
#' 
#' @param cpg_table A data.frame with methylation values for CpG sites where CpG names are rownames and column names are samples. 
#' @param feature A numeric vector for a feature of interest with sample names as observation names.
#' @param method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor().
#' @return A data.frame with the correlation values for all CpGs in table2 with feature
#' @export
cpg_feature_cor = function(cpg_table, feature, method = "pearson"){

    # Check that allowed values provided for correlation method 
    match.arg(arg = method, choices = c("pearson", "p", "spearman", "s"), several.ok = F)

    # Find sample names common to both tables and subset tables so that they contain these samples
    common_names = intersect(names(feature), colnames(cpg_table))
    if(length(common_names) == 0){stop("There are no samples in common between feature and cpg_table")}
    if(length(common_names) < length(names(feature))){warning("There are samples in the cpg_table which are not in feature. These are being ignored.")}
    if(length(common_names) < length(colnames(cpg_table))){warning("There are samples in feature which are not in cpg_table. These are being ignored.")}
    feature = feature[common_names]
    cpg_table = cpg_table[, common_names]

    # Transpose cpg_table so that features are now columns
    cpg_table = data.frame(t(cpg_table[, , drop = F]))
    
    # Calculate the correlation results for all CpGs
    cor_test_results = lapply(cpg_table, function(x) tryCatch(suppressWarnings(cor.test(x, feature, method = method, use = "c")), error = function(err) NA))
    
    # Create a data.frame with the correlation results 
    cor_df = data.frame(correlation = sapply(cor_test_results, function(x) tryCatch({x$estimate}, error = function(e) NA)), 
      p_value = sapply(cor_test_results, function(x) tryCatch({x$p.value}, error = function(e) NA)), row.names = names(cor_test_results))
    
    # Replace . with : as R automatically changes : to . in column names
    row.names(cor_df) = gsub("\\.", ":", row.names(cor_df))
    
    # Calculate FDR q-value 
    cor_df$q_value = p.adjust(cor_df$p_value, method = "fdr")
    
    return(cor_df)
}

#' Summarize CpG methylation across genomic windows 
#' 
#' @param cpg_values A data.frame with the methylation values of CpGs for different samples such as returned by extract_cpg_values_from_methrix.
#' @param window_size The size of the windows in which to bin CpGs. Default is 1000. 
#' @param summary_function Summary function to apply to CpG methylation values in windows. Default is mean. Can be any function which returns a single value from a numeric vector. 
#' @param min_cpgs_per_window The minimum number of non-missing CpG methylation values needed to assign a value to a window. Windows with less than this value are given a value of NA. Default is 1.
#' @return A data.frame with a summary of CpG methylation across genomic windows for each sample. Row names are the names of the genomic windows.
#' @export
cpg_window_summary = function(cpg_values, window_size = 1000, summary_function = mean, min_cpgs_per_window = 1) {
  
  # Make sure summary_function is a function
  if(!is.function(summary_function)){stop("summary_function should be the unquoted name of a summary statistic function e.g. mean, median, sd")}
  
  # Get chromosome of each CpG 
  chromosome = gsub(":.*", "", row.names(cpg_values))
  
  # Get start coordinate of each CpG
  cpg_position = as.numeric(gsub(".*:", "", row.names(cpg_values)))
  
  # Get genomic window of each CpG using window size and the number of different windows
  genome_window = paste0(chromosome, ":", as.integer(round(cpg_position/window_size)*window_size))
  
  # Calculate the summary value for each window in each sample.
  # Gives NA for any window where number of non-missing CpGs is not above min_cpgs_per_window.
  window_values = lapply(cpg_values, function(x) 
    sapply(split(x, genome_window), function(y) ifelse(sum(!is.na(y)) >= min_cpgs_per_window, summary_function(y, na.rm = T), NA)))
  
  # Combine results into a data.frame
  window_values = data.frame(window_values)
  return(window_values)
}

#' Summarize methylation of CpGs overalpping regions
#'
#' @param methrix_object A methrix object
#' @param genomic_regions Genomic regions to calculate mean CpG methylation for.
#' @param genomic_regions_names Names to give genomic_regions. Default is to name them region_1, region_2, etc. if no names provided.
#' @param summary_function How methylation in regions should be summarized. Can be one of "mean", "sum", "max" or "min". Default is "mean".
#' @param overlap_type Type of overlap_type to use with CpG region overlaps. Same as the options for methrix::get_region_summary(), but default is "any". 
#' @return A data.frame with the mean methylation of CpGs in each region for each sample
#' @export
summarize_region_methylation = function(methrix_object, genomic_regions, genomic_regions_names = NULL, summary_function = "mean", overlap_type = "any"){
  
  # Check that allowed values provided for summary_function and overlap_type
  match.arg(arg = summary_function, choices = c("mean", "sum", "max", "min"), several.ok = F)
  match.arg(arg = overlap_type, choices = c("within", "any"), several.ok = F)   
  
  # Found common sequences to methrix and granges object
  common_sequences = gtools::mixedsort(intersect(unique(seqnames(genomic_regions)), unique(rowData(methrix_object)$chr)))

  # Add names to genomic_regions if they are not already present
  if(is.null(genomic_regions_names)){
    names(genomic_regions) = paste0("region_", 1:length(genomic_regions))
  } else {
    if(length(genomic_regions_names) != length(genomic_regions)){
      stop("genomic_regions_names must be the same length as genomic_regions")
    } else {
      names(genomic_regions) = genomic_regions_names
    }
  }
  
  # For each sequence get methylation of the associated regions
  `%do%` = foreach::`%do%`
  region_methylation = foreach::foreach(chrom = common_sequences) %do% {

    # Print sequence which is being started
    system2("echo", paste("Summarizing regions on", chrom))

    # Filter GRanges for sequences on chrom and subset methrix object for chrom
    regions_on_chrom = genomic_regions[seqnames(genomic_regions) == chrom]
    methrix_for_chrom = suppressMessages(subset_methrix(methrix_object, contigs = chrom))
    invisible(gc())

    # Get the mean methylation of CpGs in the region on chrom
    result = suppressMessages(data.frame(get_region_summary(m = methrix_for_chrom, regions = regions_on_chrom, overlap_type = overlap_type, how = summary_function)))

    # Add the names of the regions to the result
    result$region_name = names(regions_on_chrom)[result$rid]
    gc()
    result
  }

  # Combine results for each sequence into a single data.frame
  region_methylation = dplyr::bind_rows(region_methylation)
  
  # Drop chr, start, end, n_overlap_CpGs and rid column
  region_methylation = dplyr::select(region_methylation, -c("chr", "start", "end", "n_overlap_CpGs", "rid"))
  
  # If there are duplicated region names, combine the results
  if(any(duplicated(region_methylation$region_name))){
    region_methylation = dplyr::bind_rows(lapply(split(dplyr::select(region_methylation, -region_name), region_methylation$region_name), function(x) 
      lapply(x, function(y) get(summary_function)(y, na.rm = T))), .id = "region_name")
  }
  
  # Change region_name to row.names and return 
  region_methylation = tibble::column_to_rownames(region_methylation, "region_name")
  
  gc()
  return(region_methylation)

}

#' Get names of CpGs in a methrix object which overlap a GRanges object 
#'
#' @param methrix_object A methrix object
#' @param genomic_regions A GRanges object
#' @export
find_cpgs_overlapping_granges = function(methrix_object, granges_object){
  `%do%` = foreach::`%do%`
  granges_object = reduce(granges_object, ignore.strand = T)
  cpg_df = rowData(methrix_object)
  cpg_chrom_list = split(cpg_df, cpg_df$chr)
  results = foreach::foreach(chrom = cpg_chrom_list, .combine = "c") %do% {
    as.character(subsetByOverlaps(makeGRangesFromDataFrame(chrom, end.field = "start"), granges_object, ignore.strand = T))
  }
  gc()
  return(results)
}

#' Sample random CpGs from a methrix object
#'
#' @param methrix_object A methrix object
#' @param n_cpgs Number of CpG sites to randomly select. Default is 1000.
#' @param regions Genomic regions to subset by. Could be a data.table with 3 columns (chr, start, end) or a GenomicRanges object
#' @param contigs Chromosome names to subset by.
#' @param samples Sample names to subset by.
#' @param overlap_type defines the type of the overlap of the CpG sites with the target region. 
#' Default value is 'any'. For detailed description, see the foverlaps function of the data.table package. 
#' @export
sample_cpgs = function(methrix_object, n_cpgs = 1000, regions = NULL, contigs = NULL, samples = NULL, overlap_type = "any"){
  
  # Filter methrix object for contigs, samples and regions 
  methrix_object = methrix::subset_methrix(methrix_object, contigs = contigs, samples = samples)
  if(!is.null(regions)){methrix_object = methrix::subset_methrix(methrix_object, regions = reduce(regions, ignore.strand = T), overlap_type = overlap_type)}
    
  # Get random indices from the methrix object
  random_cpgs = data.frame(rowData(methrix_object)[sample(nrow(methrix_object), n_cpgs), ])
  random_cpgs = sort(GenomicRanges::makeGRangesFromDataFrame(random_cpgs, end.field = "start"))
  methrix_object = methrix::subset_methrix(methrix_object, regions = random_cpgs, overlap_type = "any")
  
  # Get values for random CpGs
  random_cpg_values = methodical:::get_matrix_updated(methrix_object, type = "M", add_loci = T)
  
  # Add CpG names as row names and remove the loci columns
  random_cpg_values = tibble::column_to_rownames(dplyr::select(dplyr::mutate(random_cpg_values, cpg_name = paste(chr, start, sep = ":")), -chr, -start, -strand), "cpg_name") 
  
}
