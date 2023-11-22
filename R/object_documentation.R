#' tubb6_tss
#'
#' The location of the TSS for TUBB6.
#'
#'@format GRanges object with 1 ranges for the TUBB6 TSS.
#'@source The TSS of the ENST00000591909 transcript. 
"tubb6_tss"

#' tubb6_meth_rse
#'
#' The location of the TSS for TUBB6.
#'
#'@format A call to create a RangedSummarizedExperiment with methylation data for 355 CpG sites within +/- 5,000 
#'base pairs of the TUBB6 TSS in 126 normal prostate samples. 
#'Should be evaluated after loading using `tubb6_meth_rse <- tubb6_meth_rse <- eval(tubb6_meth_rse)` to restore the RangedSummarizedExperiment.
#'@source WGBS data from 'Li, Jing, et al. "A genomic and epigenomic atlas of prostate cancer in Asian populations." 
#'Nature 580.7801 (2020): 93-99.' 
"tubb6_meth_rse"

#' tubb6_transcript_counts
#'
#' Transcript counts for TUBB6 in normal prostate samples. 
#'
#'@format A data.frame with normalized transcript counts for TUBB6 in 126 normal prostate samples. 
#'@source RNA-seq data from 'Li, Jing, et al. "A genomic and epigenomic atlas of prostate cancer in Asian populations." 
#'Nature 580.7801 (2020): 93-99.' 
"tubb6_transcript_counts"

#' tubb6_cpg_meth_transcript_cors
#'
#' A data.frame with the methylation-transcription correlation results for CpGs around the TUBB6 TSS.
#'
#'@format A ggplot object.  
"tubb6_cpg_meth_transcript_cors"

#' tubb6_tmrs
#'
#' TMRs identified for TUBB6
#'
#'@format A GRanges object with two ranges. 
"tubb6_tmrs"

#' tubb6_correlation_plot
#'
#' A plot of the correlation values between methylation-transcription correlations for CpG 
#' sites around the TUBB6 TSS.
#'
#'@format A ggplot object.  
"tubb6_correlation_plot"

#' hg38_cpgs_subset
#'
#' All the CpG sites within the first one million base pairs of chromosome 1.
#'
#'@format A GRanges object.  
"hg38_cpgs_subset"

#' TumourMethDatasets
#'
#' A table describing the datasets available from TumourMethData.
#'
#'@format A data.frame with one row for each dataset
"TumourMethDatasets"