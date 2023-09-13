#' genome_annotation_hg38
#'
#' Comprehensive annotation of the standard chromosomes of the hg38 genome build including 
#' gene annotation, structural annotation, regulatory element annotation and repeat annotation. 
#' 
#'
#'@format GRanges object with 9,034,338 ranges and one metadata column region_type describing the region. 
#'@source Created by combining gene annotation from Gencode v38, regulatory feature annotation from Ensembl release 109 and
#'masked CpG island data and repeat annotation from UCSC genome browser annotation database. 
"genome_annotation_hg38"

#' Annotate GRanges
#' 
#' @param genomic_regions A GRanges object to be annotated
#' @param annotation_ranges A GRanges object with regions for different features e.g. introns, exons, enhancers. 
#' @param annotation_column Name of the metadata column of annotation_ranges indicating the feature group that regions belong to. Default is "region_type".
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating intersections. Default is TRUE.
#' @param overlap_measure One of "absolute", "proportion" or "jaccard" indicating whether to calculate 
#' the absolute size of the intersection in base pairs, the proportion base paris of gr1 overlapping gr2 
#' or the Jaccard index of the intersection in terms of base pairs. Default value is "absolute".
#' @return A numeric vector with the overlap measure for genomic_regions with each type of region in annotation_ranges
#' @export
annotate_granges = function(genomic_regions, annotation_ranges, annotation_column = "region_type", ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # If annotation_ranges is "hg38", load genome_annotation_hg38. 
  # Otherwise, check that provided annotation_ranges is a GRanges and has a metadata column matching annotation_column
  if(!is(annotation_ranges, "GRanges")){stop("annotation_ranges must be a GRanges object")}
  if(!annotation_column %in% names(mcols(annotation_ranges))){
      stop(paste(annotation_ranges, "is not the name of a metadata column of annotation_ranges"))
  }
  
  # Split annotation_ranges into a list using annotation_column
  annotation_ranges = split(annotation_ranges, mcols(annotation_ranges)[[annotation_column]]); gc()
  
  # Calculate the intersection between genomic_regions and different groups of regions defined by annotation_ranges. 
  annotation_overlaps = sapply(annotation_ranges, function(x) 
    methodical:::calculate_regions_intersections(gr1 = genomic_regions, gr2 = x, ignore.strand = ignore.strand, overlap_measure = overlap_measure))
  
  return(annotation_overlaps)
  
}
