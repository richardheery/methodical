#' Annotate GRanges
#' 
#' @param genomic_regions A GRanges object to be annotated
#' @param annotation_ranges A GRanges object with regions for different features e.g. introns, exons, enhancers. 
#' @param annotation_column Name of the metadata column of annotation_ranges indicating the feature group that regions belong to. Default is "region_type".
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param overlap_measure One of "absolute", "proportion" or "jaccard" indicating whether to calculate 
#' the absolute size of the intersection in base pairs, the proportion base pairs of gr1 overlapping gr2 
#' or the Jaccard index of the intersection in terms of base pairs. Default value is "absolute".
#' @return A numeric vector with the overlap measure for genomic_regions with each type of region in annotation_ranges
#' @export
annotateGRanges <- function(genomic_regions, annotation_ranges, annotation_column = "region_type", ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # Check that inputs have the correct data type
  stopifnot(is(genomic_regions, "GRanges"), is(annotation_ranges, "GRanges"),
    is(annotation_column, "character"), is(ignore.strand, "logical"), is(overlap_measure, "character"))
    
  # Check that provided annotation_ranges has a metadata column matching annotation_column
  if(!annotation_column %in% names(mcols(annotation_ranges))){
      stop(paste(annotation_ranges, "is not the name of a metadata column of annotation_ranges"))
  }
  
  # Split annotation_ranges into a list using annotation_column
  annotation_ranges <- split(annotation_ranges, mcols(annotation_ranges)[[annotation_column]]); gc()
  
  # Calculate the intersection between genomic_regions and different groups of regions defined by annotation_ranges. 
  annotation_overlaps <- sapply(annotation_ranges, function(x) 
    .calculate_regions_intersections(gr1 = genomic_regions, gr2 = x, ignore.strand = ignore.strand, overlap_measure = overlap_measure))
  
  return(annotation_overlaps)
  
}
