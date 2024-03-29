#' Extract values for methylation sites overlapping genomic regions from a methylation RSE. 
#' 
#' @param meth_rse A RangedSummarizedExperiment for methylation data.
#' @param genomic_regions A GRanges object. If set to NULL, returns all methylation sites in meth_rse
#' @param samples_subset Optional sample names used to subset meth_rse.
#' @param assay_number The assay from meth_rse to extract values from. Default is the first assay.
#' @return A data.frame with the methylation site values for all sites in meth_rse which overlap genomic_ranges. 
#' Row names are the coordinates of the sites as a character vector. 
#' @export
#' @examples 
#' # Load sample RangedSummarizedExperiment with CpG methylation data
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#' 
#' # Create a sample GRanges object to use
#' test_region <- GRanges("chr18:12305000-12310000")
#' 
#' # Get methylation values for CpG sites overlapping HDAC1 gene
#' test_region_methylation <- methodical::extractGRangesMethSiteValues(meth_rse = tubb6_meth_rse, genomic_regions = test_region)
extractGRangesMethSiteValues <- function(meth_rse, genomic_regions = NULL, samples_subset = NULL, assay_number = 1){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), 
    is(genomic_regions, "GRanges") | is.null(genomic_regions), is(samples_subset, "character") | is.null(samples_subset),
    is(assay_number, "numeric"))
  
  # If samples_subset provided, check that all samples present in meth_rse and then subset meth_rse for those samples
  if(!is.null(samples_subset)){
    if(any(!samples_subset %in% colnames(meth_rse))){
      stop("Some sample names in samples_subset not present in meth_rse")
    } else {
      meth_rse <- meth_rse[, samples_subset]
    }
  }
  
  # Set genomic_regions to all regions in meth_rse if not provided
  if(is.null(genomic_regions)){
    message("genomic_regions not provided so extracting all methylation values from meth_rse")
    genomic_regions <- rowRanges(meth_rse)
  }
  
  # Subset meth_rse for sites overlapping genomic_ranges
  meth_rse_subset <- IRanges::subsetByOverlaps(meth_rse, genomic_regions)
  
  # Extract GRanges from meth_rse_subset and convert into a character vector
  meth_sites_subset <- as.character(SummarizedExperiment::rowRanges(meth_rse_subset))
  
  # Extract methylation values from meth_rse_subset
  meth_values_subset <- as.data.frame(as.matrix(assay(meth_rse_subset, assay_number)))
  
  # Set row names as the names of the methylation sites and return
  row.names(meth_values_subset) <- meth_sites_subset
  return(meth_values_subset)
  
}

#' Randomly sample methylation sites from a methylation RSE. 
#' 
#' @param meth_rse A RangedSummarizedExperiment for methylation data.
#' @param n_sites Number of sites to randomly sample. Default is 1000.
#' @param genomic_ranges_filter An optional GRanges object used to first subset meth_rse. 
#' Sites will then be chosen randomly from those overlapping these ranges.
#' @param invert_filter TRUE or FALSE indicating whether to invert the genomic_ranges_filter so 
#' as to exclude sites overlapping these regions. Default value is FALSE.
#' @param samples_subset Optional sample names used to subset meth_rse.
#' @param assay_number The assay from meth_rse to extract values from. Default is the first assay.
#' @return A data.frame with the methylation site values for all sites in meth_rse which overlap genomic_ranges. 
#' Row names are the coordinates of the sites as a character vector. 
#' @export
#' @examples 
#' # Load sample RangedSummarizedExperiment with CpG methylation data
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#' 
#' # Create a sample GRanges object to use to mask tubb6_meth_rse
#' mask_ranges <- GRanges("chr18:12305000-12310000")
#' 
#' # Get 20 random CpG sites outside mask_ranges
#' random_cpgs <- methodical::sampleMethSites(tubb6_meth_rse, n_sites = 20, genomic_ranges_filter = mask_ranges, 
#'   invert_filter = TRUE)
#' 
#' # Check that no CpGs overlap repeats
#' intersect(rowRanges(random_cpgs), mask_ranges)
#' 
sampleMethSites <- function(meth_rse, n_sites = 1000, genomic_ranges_filter = NULL, 
  invert_filter = FALSE, samples_subset = NULL, assay_number = 1){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), is(n_sites, "numeric") & n_sites >= 1,
    is(genomic_ranges_filter, "GRanges") | is.null(genomic_ranges_filter), 
    S4Vectors::isTRUEorFALSE(invert_filter), is(samples_subset, "character") | is.null(samples_subset),
    is(assay_number, "numeric"))
  
  # If genomic_ranges_filter provided, subset meth_rse with it
  if(!is.null(genomic_ranges_filter)){
    meth_rse <- IRanges::subsetByOverlaps(meth_rse, genomic_ranges_filter, invert = invert_filter)
  }
  
  # Randomly sample specified number of sites from meth_rse
  sites <- sample(nrow(meth_rse), n_sites, replace = FALSE)
  
  # Subset meth_rse for random sites
  meth_rse_sites <- meth_rse[sites, ]
  
  # Subset for samples if specified
  if(!is.null(samples_subset)){
    meth_rse_sites <- meth_rse_sites[, samples_subset]
  }
  return(meth_rse_sites)
  
}
  
#' Liftover rowRanges of a RangedSummarizedExperiment for methylation data from one genome build to another
#' 
#' Removes methylation sites which cannot be mapped to the target genome build and those which result in 
#' many-to-one mappings. Also removes one-to-many mappings by default and can remove sites which do not
#' map to allowed regions in the target genome e.g. CpG sites. 
#' 
#' @param meth_rse A RangedSummarizedExperiment for methylation data
#' @param chain A "Chain" object to be used with rtracklayer::liftOver
#' @param remove_one_to_many_mapping TRUE or FALSE indicating whether to remove regions in the source genome 
#' which map to multiple regions in the target genome. Default is TRUE.
#' @param permitted_target_regions An optional GRanges object used to filter the rowRanges by overlaps after liftover, 
#' for example CpG sites from the target genome. Any regions which do not overlap permitted_target_regions will be removed.  
#' GRangesList to GRanges if all remaining source regions can be uniquely mapped to the target genome. 
#' @return A RangedSummarizedExperiment with rowRanges lifted over to the genome build indicated by chain. 
#' @examples
#' # Load sample RangedSummarizedExperiment with CpG methylation data
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#'   
#' # Get CpG sites for hg19
#' hg19_cpgs <- methodical::extractMethSitesFromGenome("BSgenome.Hsapiens.UCSC.hg19")
#' 
#' # Get liftover chain for mapping hg38 to hg19
#' library(AnnotationHub)
#' ah <- AnnotationHub()
#' chain <- ah[["AH14108"]]
#'   
#' # Liftover tubb6_meth_rse from hg38 to hg19, keeping only sites that were mapped to CpG sites in hg19
#' tubb6_meth_rse_hg19 <- methodical::liftoverMethRSE(tubb6_meth_rse, chain = chain, 
#'   permitted_target_regions = hg19_cpgs)
#' @export
liftoverMethRSE <- function(meth_rse, chain, remove_one_to_many_mapping = TRUE, permitted_target_regions = NULL){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_rse, "RangedSummarizedExperiment"), is(chain, "Chain"),
    S4Vectors::isTRUEorFALSE(remove_one_to_many_mapping), 
    is(permitted_target_regions, "GRanges") | is.null(permitted_target_regions))
  
  # Liftover rowRanges for meth_rse using specified liftover chain file
  liftover_ranges <- rtracklayer::liftOver(SummarizedExperiment::rowRanges(meth_rse), chain)
  
  # Put seqlevels of liftover_ranges in the same order as meth_rse
  GenomeInfoDb::seqlevels(liftover_ranges) <- GenomeInfoDb::seqlevels(meth_rse)
  
  # Initialize selected regions to all liftover_ranges
  selected_ranges <- seq_along(liftover_ranges)
  
  # Get the number of regions in the target genome each region in the source genome matches to
  mappings_count <- lengths(liftover_ranges)
  
  # Remove non-mapping regions from selected_ranges if specified
  non_mapping_regions <- which(mappings_count == 0)
  message(paste(length(non_mapping_regions), "non-mapping sites removed"))
  selected_ranges <- setdiff(selected_ranges, non_mapping_regions)
  
  # Remove one-to-many mapping regions from selected_ranges if specified
  if(remove_one_to_many_mapping){
    one_to_many_mapping_regions <- which(mappings_count > 1)
    message(paste(length(one_to_many_mapping_regions), "one-to-many mapping sites removed"))
    selected_ranges <- setdiff(selected_ranges, one_to_many_mapping_regions)
  }
  
  # Remove many-to-one mapping regions from selected_ranges if specified
  self_overlaps <- GenomicRanges::countOverlaps(liftover_ranges, liftover_ranges)
  many_to_one_mapping_regions <- which(self_overlaps > 1)
  message(paste(length(many_to_one_mapping_regions), "many-to-one mapping sites removed"))
  selected_ranges <- setdiff(selected_ranges, many_to_one_mapping_regions)
  
  # Identify regions which overlap permitted_target_regions is provided
  if(!is.null(permitted_target_regions)){
    target_regions_overlaps <- GenomicRanges::countOverlaps(liftover_ranges, permitted_target_regions)
    non_target_overlaps <- which(target_regions_overlaps < 1)
    message(paste(length(non_target_overlaps), "sites not overlapping permitted target regions removed"))
    selected_ranges <- setdiff(selected_ranges, non_target_overlaps)
  }
  
  # Subset meth_rse for selected rows and update rowRanges
  meth_rse <- meth_rse[selected_ranges, ]
  SummarizedExperiment::rowRanges(meth_rse) <- unlist(liftover_ranges[selected_ranges]) #
  
  # Return meth_rse 
  return(meth_rse)
  
}

#' Mask regions in a ranged summarized experiment
#'
#' @param rse A RangedSummarizedExperiment.
#' @param mask_ranges Either a GRanges with regions to be masked in all samples (e.g. repetitive sequences) 
#' or a GRangesList object with different regions to mask in each sample (e.g. mutations). If using a GRangesList object, names of the list
#' elements should be the names of samples in rse.
#' @param assay_number Assay to perform masking. Default is first assay
#' @return A RangedSummarizedExperiment with the regions present in mask_ranges masked
#' @export
#' @examples 
#' # Load sample RangedSummarizedExperiment with CpG methylation data
#' data(tubb6_meth_rse, package = "methodical")
#' tubb6_meth_rse <- eval(tubb6_meth_rse)
#' 
#' # Create a sample GRanges object to use to mask tubb6_meth_rse
#' mask_ranges <- GRanges("chr18:12305000-12310000")
#' 
#' # Mask regions in tubb6_meth_rse
#' tubb6_meth_rse_masked <- methodical::maskRangesInRSE(tubb6_meth_rse, mask_ranges)
#' 
#' # Count the number of NA values before and after masking
#' sum(is.na(assay(tubb6_meth_rse)))
#' sum(is.na(assay(tubb6_meth_rse_masked)))
maskRangesInRSE <- function(rse, mask_ranges, assay_number = 1){
  
  # Check that inputs have the correct data type
  stopifnot(is(rse, "RangedSummarizedExperiment"), 
    is(mask_ranges, "GRanges") | is(mask_ranges, "GRangesList"),
    is(assay_number, "numeric"))
  
  # Create a copy of rse for masking
  rse_masked <- rse
  
  # If mask_ranges is a GRangesList loop through each individual GRanges
  if(is(mask_ranges, "GRangesList")){
  
    # Check that names of mask_ranges match those of samples in rse
    if(length(setdiff(names(mask_ranges), colnames(rse))) > 0){
      stop("Names of mask_ranges should be present in colnames of rse")
    }
    
    # Loop through each GRanges object and mask the overlapping ranges in the corresponding sample in rse
    for(sample in names(mask_ranges)){
      
      # Print name of sample being masked
      message(sprintf("masking regions in sample %s", sample))
      
      # Find row indices for regions overlapping regions to be masked for sample
      mask_indices <- unique(queryHits(findOverlaps(rse_masked, mask_ranges[[sample]])))
      
      # Set values to be masked as NA
      assay(rse_masked, assay_number)[mask_indices, sample] <- NA
    }
    
  } else if(is(mask_ranges, "GRanges")){
    
    # Find row indices for regions overlapping regions to be masked in all samples
    mask_indices <- unique(queryHits(findOverlaps(rse_masked, mask_ranges)))
    
    # Set values to be masked as NA
    SummarizedExperiment::assay(rse_masked, assay_number)[mask_indices, ] <- NA
    
    
  } else {
    
    stop("mask_ranges should be either a GRanges or a GRangesList")
    
  }
  
    # Return rse_masked
    return(rse_masked)
  
}
