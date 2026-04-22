# Load required packages
library(methodical)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Athaliana.TAIR.TAIR9)

# Load TSS GRanges and get region +/- 5KB of the TSS
wget("")
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

# Find TUBB6 TMRs
tubb6_tmrs = methodical::findTMRs(list(ENST00000591909 = tubb6_cpg_meth_transcript_cors))

# Change tubb6_meth_rse to a call to the tubb6_meth_rse directory located in the package extdata
tubb6_meth_rse = quote(HDF5Array::loadHDF5SummarizedExperiment(system.file('extdata/tubb6_meth_rse', package = 'methodical')))

# Get CpG islands using annotatr
hg38_cpg_islands = annotatr::build_annotations(genome = "hg38", annotations = "hg38_cpgs")

# Get 1 MB sequence from human chr11 from hg38 and Arabidopsis chr4 from BSgenome.Athaliana.TAIR.TAIR9
chr1_subset_hg38 = setNames(DNAStringSet(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr1")[1:1000000]), "chr1")
chr11_hg38 = setNames(DNAStringSet(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr11")), "chr11")
chr4_subset_a_thal = setNames(DNAStringSet(Biostrings::getSeq(BSgenome.Athaliana.TAIR.TAIR9, "Chr4")[1:1000000]), "Chr4")

# Extract hg38 CpGs from chr11_subset_hg38 and filter for those in vicinity of GSTP1
chr11_hg38_cpgs = methodical::extractMethSitesFromGenome(chr11_hg38, stranded = T)
chr11_subset_hg38_cpgs = subsetByOverlaps(chr11_hg38_cpgs, GRanges("chr11:67578812-67588812"))

# Get hg19 CpGs on chr18
hg19_chr18 = setNames(DNAStringSet(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, "chr18")), "chr18")
hg19_chr18_cpgs = methodical::extractMethSitesFromGenome(hg19_chr18)

# Add objects to package
usethis::use_data(tubb6_tss, overwrite = T, compress = "xz")
usethis::use_data(tubb6_meth_rse, overwrite = T, compress = "xz")
usethis::use_data(tubb6TranscriptCounts, overwrite = T, compress = "xz")
usethis::use_data(tubb6_cpg_meth_transcript_cors, overwrite = T, compress = "xz")
usethis::use_data(tubb6_tmrs, overwrite = T, compress = "xz")
usethis::use_data(tubb6_correlation_plot, overwrite = T, compress = "xz")
usethis::use_data(chr11_subset_hg38_cpgs, overwrite = T, compress = "xz")
usethis::use_data(hg38_cpg_islands, overwrite = T, compress = "xz")
usethis::use_data(chr1_subset_hg38, overwrite = T, compress = "xz")
usethis::use_data(chr4_subset_a_thal, overwrite = T, compress = "xz")
usethis::use_data(hg19_chr18_cpgs, overwrite = T, compress = "xz")