#' Create a GRanges with methylation sites of interest from a BSgenome. 
#'
#' @param genome A BSgenome object or the name of one.  
#' @param pattern A pattern to match in bsgenome. Default is "CG".
#' @param plus_strand_only A logical value indicating whether to only return matches on "+" strand, 
#' avoiding returning duplicate hits for palindromic sequences e.g. CG. Default is TRUE.
#' @param meth_site_position Which position in the pattern corresponds to the methylation site of interest. 
#' Default is the first position.
#' @param standard_sequences_only A logical value indicating whether to only return sites 
#' on standard sequences (those without "-" in their names). Default is TRUE. 
#' @return A GRanges object with genomic regions matching the pattern.
#' @export
#' @examples 
#' # Get human CpG sites for hg38 genome build
#' hg38_cpgs = methodical::extract_meth_sites_from_genome("BSgenome.Hsapiens.UCSC.hg38")
#' 
#' # Find CHG sites in Arabidopsis thaliana
#' arabidopsis_cphpgs = methodical::extract_meth_sites_from_genome("BSgenome.Athaliana.TAIR.TAIR9", pattern = "CHG")
extract_meth_sites_from_genome = function(genome, pattern = "CG", plus_strand_only = TRUE, 
  meth_site_position = 1, standard_sequences_only = TRUE){
  
  # If genome is a character, try to load genome with that name
  if(is.character(genome)){genome = BSgenome::getBSgenome(genome)}
  
  # Check that meth_site_position is not greater than the length of the pattern
  if(meth_site_position > nchar(pattern)){
    stop("meth_site_position cannot be greater than the number of characters in pattern")
  }
  
  # Find sites matching pattern in genome
  meth_sites_gr = Biostrings::vmatchPattern(pattern, genome, fixed = "subject")
  
  # Filter for matches on "+" strand if specified
  if(plus_strand_only){
    meth_sites_gr = meth_sites_gr[GenomicRanges::strand(meth_sites_gr) == "+"]
  }
  
  # Add seqinfo to GRanges
  GenomeInfoDb::seqinfo(meth_sites_gr) = GenomeInfoDb::seqinfo(genome)
  
  # Subset for standard sequences if specified
  if(standard_sequences_only){
    standard_sequences = grep("_", names(genome), invert = TRUE, value = TRUE)
    GenomeInfoDb::seqlevels(meth_sites_gr, pruning.mode = "coarse") = standard_sequences
  }
  
  # Set strand as "*" if plus_strand_only is TRUE 
  if(plus_strand_only){GenomicRanges::strand(meth_sites_gr) = "*"}
  
  # Resize meth_sites_gr so that it covers just the meth site position
  meth_sites_gr = resize(meth_sites_gr, meth_site_position, fix = "start")
  start(meth_sites_gr) = end(meth_sites_gr)
  
  # Return meth_sites_gr
  return(meth_sites_gr)
  
}

#' Expand GRanges
#'
#' Expand ranges in a GRanges object upstream and downstream by specified numbers of bases, taking account of strand.
#' Unstranded ranges are treated like they on the "+" strand. 
#' If any of the resulting ranges are out-of-bounds given the seqinfo of genomic_regions, they will be trimmed using trim().
#'
#' @param genomic_regions A GRanges object
#' @param upstream Number of bases to add upstream of each region in genomic_regions. 
#' Must be numeric vector of length 1 or else equal to the length of genomic_regions. Default value is 0. 
#' Negative values result in upstream end of regions being shortened, however the width of the resulting regions cannot be less than zero. 
#' @param downstream Number of bases to add downstream of each region in genomic_regions. Negative values result in downstream end of regions being shortened. 
#' Must be numeric vector of length 1 or else equal to the length of genomic_regions. Default value is 0.
#' Negative values result in upstream end of regions being shortened, however the width of the resulting regions cannot be less than zero. 
#' @return A GRanges object
#' @export
#' @examples
#' # Create a GRanges 
#' gr = GenomicRanges::GRanges(c("chr1:1000:+", "chr1:2000:-"))
#' 
#' # Add 500 bp upstream and 200 downstream
#' gr_expanded = methodical::expand_granges(gr, 500, 200)
expand_granges = function(genomic_regions, upstream = 0, downstream = 0) {
  
  # Check that genomic_regions is a GRanges object
  if(!is(genomic_regions, "GRanges")){stop("genomic_regions must be a GRanges object")}
  
  # Check that upstream and downstream are vectors of either length 1 or with the same length as genomic_regions
  if(!length(upstream) %in% c(1, length(genomic_regions))){
    stop("upstream should be a vector of length 1 or the length of genomic_regions")}
  if(!length(downstream) %in% c(1, length(genomic_regions))){
    stop("downstream should be a vector of length 1 or the length of genomic_regions")}
  
  # Check if any regions would have negative widths after adjustment
  if(any(width(genomic_regions) + upstream + downstream < 0)){
    stop("Some regions would have a negative width after adjustment. This is not permitted.")
  }
  
  ## Save names of genomic_regions
  genomic_regions_names = names(genomic_regions)
  
  # Check for each range if it's on the negative or positive strand
  strand_is_minus = as.character(GenomicRanges::strand(genomic_regions)) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  
  # Create vectors with the start and end sites of genomic_regions
  genomic_regions_starts = start(genomic_regions)
  genomic_regions_ends = end(genomic_regions)
  
  # Adjust ranges based on whether they are on the positive or negative strand
  genomic_regions_starts[on_plus] = genomic_regions_starts[on_plus] - upstream
  genomic_regions_starts[on_minus] = genomic_regions_starts[on_minus] - downstream
  genomic_regions_ends[on_plus] = genomic_regions_ends[on_plus] + downstream
  genomic_regions_ends[on_minus] = genomic_regions_ends[on_minus] + upstream
  
  # Store metadata from genomic_regions
  genomic_regions_mcols = mcols(genomic_regions)
  
  # Recreate genomic_regions with new starts and ends
  genomic_regions = GRanges(seqnames = seqnames(genomic_regions), 
    ranges = IRanges(genomic_regions_starts, genomic_regions_ends), 
    strand = GenomicRanges::strand(genomic_regions))
  
  # Restore metadata
  mcols(genomic_regions) = genomic_regions_mcols
  
  # Remove any out-of-bounds regions, add names back and return genomic_regions
  genomic_regions = GenomicRanges::trim(genomic_regions)
  names(genomic_regions) = genomic_regions_names
  return(genomic_regions)
} 

#' Calculate distances of query GRanges upstream or downstream of subject GRanges
#' 
#' Upstream and downstream are relative to the strand of subject_gr. 
#' Unstranded regions are treated the same as regions on the "+" strand. 
#'
#' @param query_gr A GRanges object
#' @param subject_gr A GRanges object. 
#' @return A numeric vector of distances
#' @export
#' @examples 
#' # Create query and subject GRanges 
#' query_gr = GenomicRanges::GRanges(c("chr1:100-1000:+", "chr1:2000-3000:-"))
#' subject_gr = GenomicRanges::GRanges(c("chr1:1500-1600:+", "chr1:4000-4500:-"))
#' 
#' # Calculate distances between query and subject
#' methodical::stranded_distance(query_gr, subject_gr)
stranded_distance = function(query_gr, subject_gr){
  
  # Check that query_gr and subject_gr are of the correct length
  if(!length(subject_gr) %in% c(1, length(query_gr))){
    stop("subject_gr should have length 1 or the same length as query_gr")}
  
  # Initialize a vector of zeros with length equal to query_gr
  d = rep(0, length(query_gr))
  
  # Find the length of the gap between query_gr and subject_gr, leaving as 0 if they overlap
  d[end(query_gr) < start(subject_gr)] = (end(query_gr) - start(subject_gr))[end(query_gr) < start(subject_gr)]
  d[start(query_gr) > end(subject_gr)] = (start(query_gr) - end(subject_gr))[start(query_gr) > end(subject_gr)]
  
  # If ranges are on different sequences, set distance as NA
  d[as.vector(seqnames(query_gr)) != as.vector(seqnames(subject_gr))] = NA
  
  # Set the strand of subject_gr as 1 if it is on the "+" strand or unstranded ("*") and -1 if it is on the "-" strand
  subject_strand = ifelse(as.character(strand(subject_gr)) == "-", -1, 1) 
  
  # Get the signed distance of query_gr from subject_gr and return
  d = d * subject_strand
  return(d)

}

#' Find locations of genomic regions relative to transcription start sites.
#'
#' @param genomic_regions A GRanges object. 
#' @param tss_gr A GRanges object with transcription start sites. Each range should have width 1. 
#' Upstream and downstream are relative to strand of tss_gr.
#' @return A GRanges object where all regions have "relative" as the sequence names.  
#' @export
#' @examples
#' # Create query and subject GRanges 
#' genomic_regions = GenomicRanges::GRanges(c("chr1:100-1000:+", "chr1:2000-3000:-"))
#' tss_gr = GenomicRanges::GRanges(c("chr1:1500:+", "chr1:4000:-"))
#' 
#' # Calculate distances between query and subject
#' methodical::ranges_relative_to_tss(genomic_regions, tss_gr)
ranges_relative_to_tss = function(genomic_regions, tss_gr){
  
  # Check that all tss ranges have width 1 and resize them with a warning if not
  if(!all(width(tss_gr) == 1)){
    warning("All regions in tss_gr should have a width of 1. Shortening each region so that it consists of only the most upstream position")
    tss_gr = GenomicRanges::resize(tss_gr, 1, fix = "start")
  }

  # Get distances from start and end of ranges in gr from tss_gr
  relative_start = methodical::stranded_distance(query_gr = resize(genomic_regions, 1, fix = "start"), subject_gr = tss_gr)
  relative_end = methodical::stranded_distance(query_gr = resize(genomic_regions, 1, fix = "end"), subject_gr = tss_gr)
  
  # Create an IRanges with the relative distances
  relative_iranges = IRanges::IRanges(pmin(relative_start, relative_end), pmax(relative_start, relative_end))
  
  # Convert IRanges to GRanges with "relative" as seqnames
  relative_granges = GenomicRanges::GRanges(seqnames = "relative", ranges = relative_iranges)
  
  # Add metadata from gr to relative_granges
  mcols(relative_granges) = mcols(genomic_regions)
  
  # Return relative_granges
  return(relative_granges)
  
}

#' Calculate the number of unique bases covered by all regions in a GRanges object
#'
#' @param gr A GRanges object
#' @return An numeric value
.count_covered_bases = function(gr){
  
  return(sum(width(reduce(gr, ignore.strand = TRUE))))

}

#' Calculate the number of bases in the intersection of two GRanges objects
#'
#' @param gr1 A GRanges object
#' @param gr2 A GRanges object
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param overlap_measure One of "absolute", "proportion" or "jaccard" indicating whether to calculate 
#' the absolute size of the intersection in base pairs, the proportion base paris of gr1 overlapping gr2 
#' or the Jaccard index of the intersection in terms of base pairs. Default value is "absolute".
#' @return An numeric value
.calculate_regions_intersections = function(gr1, gr2, ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # Check allowed value provided for overlap_measure
  match.arg(overlap_measure, c("absolute", "proportion", "jaccard"))
  
  # Create GRanges with the intersection and union of gr1 and gr2
  intersection = GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
  union = c(gr1, gr2)
  
  # Caluclate proportion, Jaccard index or absolute overlap depending on overlap_measure
  if(overlap_measure == "proportion"){
    return(.count_covered_bases(intersection)/.count_covered_bases(gr1))
  } else if(overlap_measure == "jaccard"){
    return(.count_covered_bases(intersection)/.count_covered_bases(union))
  } else {
    return(.count_covered_bases(intersection))
  }

}

#' Create a GRanges with random regions from a genome
#' 
#' Can constrain the random regions so that they do not overlap each other and/or an optional set of masked regions.
#' Random regions which do meet these constraints will be discarded and new ones generated until the desired number 
#' of regions has been reached or a maximum allowed number of attempts has been made. 
#' After the maximum number of allowed attempts, the created random regions meeting the constraints up to that point will be returned. 
#' Any random regions that are out-of-bounds relative to their sequence length are trimmed before being returned. 
#'
#' @param genome A BSgenome object. 
#' @param n_regions Number of random regions to create. Default is 1000. 
#' @param region_widths The widths of the random regions. Widths cannot be negative. 
#' Can be just a single value if all regions are to have the same widths. Default is 1000.
#' @param sequences The names of sequences to create random regions on. Default is to use all standard sequences (those without "_" in their name)
#' @param all_sequences_equally_likely A logical value indicating if the probability of creating random regions on a sequence should be the same for each sequence.
#' Default is FALSE, indicating to make the probability proportional to a sequences length.
#' @param stranded A logical value indicating if created regions should have a strand randomly assigned. Default is FALSE, indicating to make unstranded regions. 
#' @param masked_regions An optional GRanges object which random regions will not be allowed to overlap. 
#' @param allow_overlapping_regions A logical value inidicating if created random regions should be allowed to overlap. Default is FALSE. 
#' @param ignore.strand A logical value indicating whether strand should be ignored when 
#' identifying overlaps between random regions with each other or with masked_regions. 
#' Only relevant if stranded is TRUE and either allow_overlapping_regions is FALSE or masked_regions is provided. Default is TRUE.
#' @param max_tries The maximum number of attempts to make to find non-overlapping regions which do not overlap masked_regions. Default value is 100. 
#' @return A GRanges object
#' @export 
#' @examples 
#' # Set random seed
#' set.seed(123)
#' 
#' # Create 10,000 random non-overlapping regions with width 1,000 for hg38
#' random_regions = methodical::create_random_regions(genome = "BSgenome.Hsapiens.UCSC.hg38", n_regions = 10000)
create_random_regions = function(genome, n_regions = 1000, region_widths = 1000, sequences = NULL, all_sequences_equally_likely = FALSE,
   stranded = FALSE, masked_regions = NULL, allow_overlapping_regions = FALSE, ignore.strand = TRUE, max_tries = 100){
  
  # Check that all region_widths are positive
  if(any(region_widths < 0)){stop("region_widths cannot contain negative values")}
  
  # If genome is a character, try to load genome with that name
  if(is.character(genome)){genome = BSgenome::getBSgenome(genome)}
  
  # If no sequences provided, use the standard sequences for the species from the provider
  if(is.null(sequences)){
    sequences = grep("_", seqnames(genome), invert = TRUE, value = TRUE)
  }
  
  # If all_sequences_equally_likely is false, make likelihood of sequences proportional to their lengths
  if(!all_sequences_equally_likely){
    sequence_probabilities = seqlengths(genome)[sequences]
  } else {
    sequence_probabilities = rep(1, length(sequences))
  }
  
  # Initialize an empty vector of GRanges and try_number to 1
  final_random_gr = GRanges()
  try_number = 1
  original_n_regions = n_regions
  
  # Keep attempting to find random regions meeting the criteria until the required number 
  # of regions are found or the maximum number of tries is reached
  while(n_regions > 0 & try_number <= max_tries){
    
    # Print the attempt number if using masked region or overlapping regions are not permitted
    if(!is.null(masked_regions) | !allow_overlapping_regions){
      message(paste("Attempt", try_number, "to find", original_n_regions, "random regions:"))
    }
    
    # Select random sequences
    random_sequences = sample(sequences, size = n_regions, prob = sequence_probabilities, replace = TRUE)
    
    # Create a data.frame
    random_gr_df = data.frame(seqnames = random_sequences, seqlengths = seqlengths(genome)[random_sequences], row.names = NULL)
    
    # Select a random start site on each sequence
    random_gr_df$start = sapply(random_gr_df$seqlengths, function(x) sample(1:x, 1))
    
    # Add the widths to each start site
    random_gr_df$end = random_gr_df$start + region_widths - 1
    
    # Randomly assign a strand if specified
    if(stranded){random_gr_df$strand = sample(c("+", "-"), nrow(random_gr_df), replace = TRUE)}
    
    # Create a GRanges from the data.frame and return and initialize a column indicating if they pass the contrainsts to TRUE
    temp_random_gr = makeGRangesFromDataFrame(random_gr_df, keep.extra.columns = FALSE)
    temp_random_gr$pass = TRUE
    
    # Indicate if any random regions overlap masked_regions
    if(!is.null(masked_regions)){
      temp_random_gr$pass = !overlapsAny(temp_random_gr, masked_regions, ignore.strand = ignore.strand)
    }
    
    # Indicate if any random regions overlap other random regions, including those previously found
    if(!allow_overlapping_regions){
      temp_random_gr$pass = 
        countOverlaps(temp_random_gr, c(temp_random_gr, final_random_gr), ignore.strand = ignore.strand) == 1 & temp_random_gr$pass
    }
    
    # Identify the passing_regions and add to final_random_gr
    passing_regions = temp_random_gr[temp_random_gr$pass]
    final_random_gr = c(final_random_gr, passing_regions)
    
    # Identify failing regions
    failing_regions = temp_random_gr[!temp_random_gr$pass]
    
    # Update n_regions to number of regions still remaining to be found
    n_regions = length(failing_regions)
    region_widths = width(failing_regions)
    
    # Print the number of random regions found
    if(!is.null(masked_regions) | !allow_overlapping_regions){message(paste("Found", length(final_random_gr), "regions"))}
    
    # Update try_number
    try_number = try_number + 1
    
  }
  
  # Add seqinfo to final_random_gr
  suppressWarnings({seqinfo(final_random_gr) = seqinfo(genome)[sequences]})
  
  # Trim out of bound regions, remove pass column and return final_random_gr
  final_random_gr = trim(final_random_gr)
  final_random_gr$pass = NULL
  return(final_random_gr)
  
}

