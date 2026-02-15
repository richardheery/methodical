#' Create a GRanges with methylation sites of interest from a BSgenome or DNAStringSet
#'
#' @param genome A BSgenome object or a DNAStringSet with names indicating the sequences.  
#' @param pattern A pattern to match in genome. Default is "CG".
#' @param stranded TRUE or FALSE indicating whether to return matches on 
#' both strands or else just the "+" strand. Strand will be set to "*" if FALSE. Default is TRUE.
#' @param standard_sequences_only TRUE or FALSE indicating whether to only return ranges 
#' on standard sequences (those without "_" in their names). Default is TRUE. 
#' @return A GRanges object with genomic regions matching the pattern.
#' @export
#' @examples 
#' # Get human CpG sites for chr18 from hg38 genome build
#' data(hg38_chr18, package = "methodical")
#' hg38_chr18_cpgs <- methodical::extractMethSitesFromGenome(hg38_chr18)
#' head(hg38_chr18_cpgs)
#' 
#' # Find CHG sites in Arabidopsis thaliana
#' data(arabidopsis_chr4, package = "methodical")
#' arabidopsis_chr4_CHG_sites <- methodical::extractMethSitesFromGenome(arabidopsis_chr4, pattern = "CHG")
#' head(head(arabidopsis_chr4_CHG_sites))
extractMethSitesFromGenome <- function(genome, pattern = "CG", 
  stranded = TRUE, standard_sequences_only = TRUE){
  
  # Check that inputs have the correct data type
  stopifnot(is(genome, "BSgenome") | is(genome, "DNAStringSet"), 
    is(pattern, "character"), S4Vectors::isTRUEorFALSE(stranded), 
    S4Vectors::isTRUEorFALSE(standard_sequences_only))
  
  # If genome is a DNASringSet, check that it has names
  if(is(genome, "DNASringSet") & is.null(names(genome))){
    stop("If genome is a DNASringSet, it must have names indicating the sequence")
  }
  
  # Extract seqinfo from genome
  seqinfo <- GenomeInfoDb::seqinfo(genome)
  
  # Convert genome to a DNAStringSet and subset for standard chromosomes if specified
  sequence_names = names(genome)
  if(standard_sequences_only){
    message("Searching only standard sequences (those without \"_\" in their names)")
    sequence_names <- grep("_", sequence_names, invert = T, value = TRUE)
    if(length(sequence_names) == 0){stop("There are no sequences which appear to be standard sequences")}
  }
  genome <- BSgenome::getSeq(genome, sequence_names)
  
  # Find sites matching pattern in genome, on both the + and - strands if stranded is TRUE or just the + strand otherwise
  if(stranded){
    
    # Find sites on plus and minus strands
    meth_sites_gr_plus <- GRanges(Biostrings::vmatchPattern(pattern, genome, fixed = "subject"), strand = "*")
    reverse_complement_pattern = Biostrings::reverseComplement(Biostrings::DNAString(pattern))
    meth_sites_gr_minus <- GRanges(Biostrings::vmatchPattern(reverse_complement_pattern, genome, fixed = "subject"), strand = "*")
    
    # Get the start and ends of sites and combine into a single GRanges and remove duplicated regions 
    meth_sites_gr_combined <- c(meth_sites_gr_plus, meth_sites_gr_minus)
    meth_sites_gr_start <- GRanges(resize(meth_sites_gr_combined, 1, fix = "start"), strand = "+")
    meth_sites_gr_end <- GRanges(resize(meth_sites_gr_combined, 1, fix = "end"), strand = "-")
    meth_sites_gr <- c(meth_sites_gr_start, meth_sites_gr_end)
    meth_sites_gr <- sort(disjoin(meth_sites_gr), ignore.strand = TRUE)
    
  } else {
    
    # Find sites only on the plus strand
    meth_sites_gr <- GRanges(Biostrings::vmatchPattern(pattern, genome, fixed = "subject"), strand = "*")
    
  }
    
  # Sort and resize ranges so that they cover just the first base
  meth_sites_gr <- sort(meth_sites_gr, ignore.strand = TRUE)
  meth_sites_gr <- resize(meth_sites_gr, 1)
  
  # Add seqinfo to GRanges, run garbage collection and return
  GenomeInfoDb::seqinfo(meth_sites_gr) <- seqinfo
  invisible(gc())
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
#' data(tubb6_tss, package = "methodical")
#' tubb6_tss
#' methodical::expand_granges(tubb6_tss, upstream = 5000, downstream = 5000)
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
  
  # Store strand and metadata from genomic_regions
  genomic_regions_strand = strand(genomic_regions)
  genomic_regions_mcols = mcols(genomic_regions)
  
  # Recreate genomic_regions with new starts and ends
  genomic_regions = GRanges(seqnames = seqnames(genomic_regions), 
    ranges = IRanges(genomic_regions_starts, genomic_regions_ends))
  
  # Restore strand and metadata
  strand(genomic_regions) = genomic_regions_strand
  mcols(genomic_regions) = genomic_regions_mcols
  
  # Remove any out-of-bounds regions and return genomic_regions
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
#' query_gr <- GenomicRanges::GRanges(c("chr1:100-1000:+", "chr1:2000-3000:-"))
#' subject_gr <- GenomicRanges::GRanges(c("chr1:1500-1600:+", "chr1:4000-4500:-"))
#' 
#' # Calculate distances between query and subject
#' methodical::strandedDistance(query_gr, subject_gr)
strandedDistance <- function(query_gr, subject_gr){
  
  # Check that inputs have the correct data type
  stopifnot(is(query_gr, "GRanges"), is(subject_gr, "GRanges"))
  
  # Check that query_gr and subject_gr are of the correct length
  if(!length(subject_gr) %in% c(1, length(query_gr))){
    stop("subject_gr should have length 1 or the same length as query_gr")}
  
  # Initialize a vector of zeros with length equal to query_gr
  d <- rep(0, length(query_gr))
  
  # Find the length of the gap between query_gr and subject_gr, leaving as 0 if they overlap
  d[end(query_gr) < start(subject_gr)] <- (end(query_gr) - start(subject_gr))[end(query_gr) < start(subject_gr)]
  d[start(query_gr) > end(subject_gr)] <- (start(query_gr) - end(subject_gr))[start(query_gr) > end(subject_gr)]
  
  # If ranges are on different sequences, set distance as NA
  d[as.vector(seqnames(query_gr)) != as.vector(seqnames(subject_gr))] <- NA
  
  # Set the strand of subject_gr as 1 if it is on the "+" strand or unstranded ("*") and -1 if it is on the "-" strand
  subject_strand <- ifelse(as.character(strand(subject_gr)) == "-", -1, 1) 
  
  # Get the signed distance of query_gr from subject_gr and return
  d <- d * subject_strand
  return(d)

}

#' Calculate the number of unique bases covered by all regions in a GRanges object
#'
#' @param gr A GRanges object
#' @return An numeric value
.count_covered_bases <- function(gr){
  
  # Check that inputs have the correct data type
  stopifnot(is(gr, "GRanges"))
  
  return(sum(width(reduce(gr, ignore.strand = TRUE))))

}

#' Calculate the number of bases in the intersection of two GRanges objects
#'
#' @param gr1 A GRanges object
#' @param gr2 A GRanges object
#' @param ignore.strand TRUE or FALSE indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param overlap_measure One of "absolute", "proportion" or "jaccard" indicating whether to calculate 
#' the absolute size of the intersection in base pairs, the proportion base pairs of gr1 overlapping gr2 
#' or the Jaccard index of the intersection in terms of base pairs. Default value is "absolute".
#' @return An numeric value
.calculate_regions_intersections <- function(gr1, gr2, ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # Check that inputs have the correct data type
  stopifnot(is(gr1, "GRanges"), is(gr2, "GRanges"), 
    S4Vectors::isTRUEorFALSE(ignore.strand), is(overlap_measure, "character"))
  
  # Check allowed value provided for overlap_measure
  match.arg(overlap_measure, c("absolute", "proportion", "jaccard"))
  
  # Create GRanges with the intersection and union of gr1 and gr2
  intersection <- GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
  union <- c(gr1, gr2)
  
  # Calculate proportion, Jaccard index or absolute overlap depending on overlap_measure
  if(overlap_measure == "proportion"){
    return(.count_covered_bases(intersection)/.count_covered_bases(gr1))
  } else if(overlap_measure == "jaccard"){
    return(.count_covered_bases(intersection)/.count_covered_bases(union))
  } else {
    return(.count_covered_bases(intersection))
  }

}