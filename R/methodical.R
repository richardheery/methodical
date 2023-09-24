#' methodical: A one-stop shop for dealing with big DNA methylation datasets
#' 
#' DNA methylation is generally considered to be associated with transcriptional silencing. 
#' However, comprehensive, genome-wide investigation of this relationship requires the evaluation of
#' potentially millions of correlation values between the methylation of individual genomic loci and 
#' expression of associated transcripts in a relatively large numbers of samples.  
#' Methodical makes this process quick and easy while keeping a low memory footprint. 
#' It also provides a novel method for identifying regions where a number of methylation sites are consistently 
#' strongly associated with transcriptional expression. In addition, Methodical enables housing DNA methylation 
#' data from diverse sources (e.g. WGBS, RRBS and methylation arrays) with a common framework, 
#' lifting over DNA methylation data between different genome builds and creating base-resolution 
#' plots of the association between DNA methylation and transcriptional activity at transcriptional start sites. 
#' 
#' @author Richard Heery
#' @docType package
#' @name methodical-package
#' @import GenomicRanges
#' @import ggplot2
#' @import SummarizedExperiment
"_PACKAGE"