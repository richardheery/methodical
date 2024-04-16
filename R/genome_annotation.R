#' Annotate GRanges
#' 
#' @param genomic_regions A GRanges object to be annotated
#' @param annotation_ranges A GRangesList object with GRanges for different features e.g. introns, exons, enhancers. 
#' @param ignore.strand TRUE or FALSE indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param overlap_measure One of "absolute", "proportion" or "jaccard" indicating whether to calculate 
#' the absolute size of the intersection in base pairs, the proportion of base pairs of 
#' genomic_ranges overlapping one of the component GRanges of annotation_ranges. 
#' or the Jaccard index of the intersection in terms of base pairs. Default value is "absolute".
#' @return A numeric vector with the overlap measure for genomic_regions with each type of region in annotation_ranges. 
#' @export
#' @examples 
#' 
#' # Load annotation for CpG islands and repetitive DNA
#' cpg_island_annotation <- annotatr::build_annotations(genome = "hg38", annotations = "hg38_cpgs")
#' cpg_island_annotation <- cpg_island_annotation[cpg_island_annotation$type == "hg38_cpg_islands"]
#' repeat_annotation_hg38 <- AnnotationHub::AnnotationHub()[["AH99003"]]
#' 
#' # Convert repeat_annotation_hg38 into a GRangesList
#' repeat_annotation_hg38 <- GRangesList(split(repeat_annotation_hg38, repeat_annotation_hg38$repClass))
#'  
#' # Calculate the proportion of base pairs in CpG islands overlapping different classes of repetitive elements
#' annotateGRanges(genomic_regions = cpg_island_annotation, annotation_ranges = repeat_annotation_hg38, overlap_measure = "proportion")
#' 
annotateGRanges <- function(genomic_regions, annotation_ranges, ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # Check that inputs have the correct data type
  stopifnot(is(genomic_regions, "GRanges"), is(annotation_ranges, "GRangesList"),
    S4Vectors::isTRUEorFALSE(ignore.strand), is(overlap_measure, "character"))
  
  # If annotation_ranges are missing names, print a warning and say that regions are being named 
  if(is.null(names(annotation_ranges))){
    message("annotation_ranges are missing names. They will be named granges_1, granges_2, etc. in the output")
    names(annotation_ranges) <- paste0("granges_", seq_along(annotation_ranges))
  }
  
  # Calculate the intersection between genomic_regions and different groups of regions defined by annotation_ranges. 
  annotation_overlaps <- sapply(annotation_ranges, function(x) 
    .calculate_regions_intersections(gr1 = genomic_regions, gr2 = x, ignore.strand = ignore.strand, overlap_measure = overlap_measure))
  
  return(annotation_overlaps)
  
}
