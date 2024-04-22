
# Load required packages
library(methodical)

# Load TSS GRanges and get region +/- 5KB of the TSS
tss_gr = readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")
tubb6_tss = tss_gr[tss_gr$transcript_id == "ENST00000591909"]
names(tubb6_tss) = tubb6_tss$transcript_id
tubb6_gr = promoters(tubb6_tss, 5000, 5000 + 1)

# Get common samples
common_normal_samples = readRDS("~/wgbs/cpgea/metadata/common_normal_samples.rds")

# Load CPGEA meth RSE
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("~/wgbs/cpgea/wgbs/cpgea_wgbs_hg38/")

# Subset cpgea_rse for tubb6 CpGs and common samples
tubb6_meth_rse = subsetByOverlaps(cpgea_rse, tubb6_gr)[, common_normal_samples]
HDF5Array::saveHDF5SummarizedExperiment(x = tubb6_meth_rse, dir = "inst/extdata/tubb6_meth_rse", replace = T)

# Get methylation values for TUBB6 in CPGEA normal samples
tubb6_tss_proximal_cpg_methylation = extractGRangesMethSiteValues(meth_rse = tubb6_meth_rse, genomic_regions = tubb6_gr,
  samples_subset = common_normal_samples)

# Get transcript values for CPGEA
cpgea_transcript_values = data.frame(data.table::fread("~/wgbs/cpgea/rnaseq/kallisto_tables//kallisto_deseq2_normalized_counts_pcg.tsv.gz"), row.names = 1)

# Get TUBB6 transcript values
tubb6_transcript_counts = cpgea_transcript_values["ENST00000591909", common_normal_samples]

# Calculate correlation values between methylation values and transcript values for TUBB6
tubb6_cpg_meth_transcript_cors = methodical::calculateMethSiteTranscriptCors(meth_rse = tubb6_meth_rse, 
  transcript_expression_table = tubb6_transcript_counts, tss_gr = tubb6_tss, expand_upstream = 5000, expand_downstream = 5000)
tubb6_cpg_meth_transcript_cors = tubb6_cpg_meth_transcript_cors$ENST00000591909
  
# Plot transcript cors for tubb6_cpg_meth_transcript_cors
tubb6_correlation_plot = methodical::plotMethSiteCorCoefs(tubb6_cpg_meth_transcript_cors)

# Plot Methdodical scores for TUBB6
tubb6_methodical_plot = methodical::plotMethodicalScores(meth_site_values = tubb6_cpg_meth_transcript_cors, smooth_scores = F)

# Add smoothed scores to plot
tubb6_smoothed_methodical_plot = methodical::plotMethodicalScores(meth_site_values = tubb6_cpg_meth_transcript_cors, 
  smooth_scores = T, smoothed_curve_colour = "hotpink2", curve_alpha = 1)

# Find TUBB6 TMRs
tubb6_tmrs = methodical::findTMRs(tubb6_cpg_meth_transcript_cors)

# Add TMRs to TUBB6 plot
tubb6_correlation_plot_with_tmrs = methodical::plotTMRs(tubb6_smoothed_methodical_plot, tmrs_gr = tubb6_tmrs)

# Change tubb6_meth_rse to a call to the tubb6_meth_rse directory located in the package extdata
tubb6_meth_rse = quote(HDF5Array::loadHDF5SummarizedExperiment(system.file('extdata/tubb6_meth_rse', package = 'methodical')))

# Get hg38 CpGs
hg38_cpgs = methodical::extractMethSitesFromGenome("BSgenome.Hsapiens.UCSC.hg38")

# Subset for CpGs within first million base pairs on chromosome 1
hg38_cpgs_subset = subsetByOverlaps(hg38_cpgs, GRanges("chr1:1-1000000"))

# Add objects to package
usethis::use_data(tubb6_tss, overwrite = T, compress = "xz")
usethis::use_data(tubb6_meth_rse, overwrite = T, compress = "xz")
usethis::use_data(tubb6TranscriptCounts, overwrite = T, compress = "xz")
usethis::use_data(tubb6_cpg_meth_transcript_cors, overwrite = T, compress = "xz")
usethis::use_data(tubb6_tmrs, overwrite = T, compress = "xz")
usethis::use_data(tubb6_correlation_plot, overwrite = T, compress = "xz")
usethis::use_data(hg38_cpgs_subset, overwrite = T, compress = "xz")
