## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 60), tidy = TRUE)
)
options(tinytex.verbose = TRUE)

## ---- message = F, warning = F------------------------------------------------
library(methodical)
library(biomaRt)
library(plotly)
library(ggpubr)

## -----------------------------------------------------------------------------
bedgraphs = list.files(system.file("extdata", "prostate_cancer_chr11_wgbs_bedgraphs",
  package = "methodical"), full.names = T)

## ---- message = F-------------------------------------------------------------
add_dummy_coverage_to_bedgraphs(bedgraphs = bedgraphs,
  output_directory = "bedgraphs_with_dummy_coverage", ncores = 4)

## ---- message = F, warning= F-------------------------------------------------
hg38_cpgs = methrix::extract_CPGs("BSgenome.Hsapiens.UCSC.hg38")

## ---- message = F-------------------------------------------------------------
# Create sample annotation
sample_annotation = data.frame(
  condition = c(rep("Normal", 10), rep("Tumour", 10)),
  patient = rep(c("1", "10", "12", "18", "19", "26", "28", "31", "32", "6"), times = 2),
  row.names = gsub("_WGBS_chr11.bg.gz", "", basename(bedgraphs))
)

## ---- message = F-------------------------------------------------------------
prostate_methrix = methrix::read_bedgraphs(files = list.files("bedgraphs_with_dummy_coverage", full.names = T),
  ref_cpgs = hg38_cpgs, chr_idx = 1, start_idx = 2, end_idx = 3, beta_idx = 4, cov_idx = 5, h5 = T,
  h5_dir = "prostate_methrix_h5", contigs = "chr11", coldata = sample_annotation)

## -----------------------------------------------------------------------------
hg38_cgis = rtracklayer::import.bed(system.file(
  "extdata", "masked_cpg_islands_hg38.bed", package = "methodical"))

## -----------------------------------------------------------------------------
prostate_chr11_cgi_methylation = extract_cpg_values_from_methrix(
  methrix = prostate_methrix,
  regions = hg38_cgis[seqnames(hg38_cgis) == "chr11"]
)

## -----------------------------------------------------------------------------
prostate_chr11_cgi_methylation_2mb_windows = cpg_window_summary(
  cpg_values = prostate_chr11_cgi_methylation,
  window_size = 2*10^6,
  min_cpgs_per_window = 10
)

## -----------------------------------------------------------------------------
sample_subset = c("N1", "N10", "N12", "T1", "T10", "T12")
plot_cpg_methylation(
  cpg_values = prostate_chr11_cgi_methylation_2mb_windows,
  geom = "geom_linepoint",
  samples = sample_subset,
  group_colours = c(RColorBrewer::brewer.pal(3, "Blues"), RColorBrewer::brewer.pal(3, "Reds")),
   title = "Chr11 CpG island Methylation in Prostate Tumour and Normal Samples")

## -----------------------------------------------------------------------------
plot_cpg_methylation(
  cpg_values = prostate_chr11_cgi_methylation_2mb_windows,
  geom = "geom_linepoint",
  samples = sample_subset,
  sample_groups = colData(prostate_methrix)[sample_subset, "condition"], group_colours = c("#00BFC4", "#F8766D"),
  title = "Chr11 CpG island Methylation in Prostate Tumour and Normal Samples"
)

## -----------------------------------------------------------------------------
plot_cpg_methylation(
  cpg_values = prostate_chr11_cgi_methylation_2mb_windows,
  geom = "geom_boxplot",
  samples = sample_subset,
  sample_groups = colData(prostate_methrix)[sample_subset, "condition"], group_colours = c("#00BFC4", "#F8766D"),
  title = "Chr11 CpG island Methylation in Prostate Tumour and Normal Samples"
)

## -----------------------------------------------------------------------------
prostate_chr11_cgi_methylation_2mb_windows$mean_tumour_normal_difference =
  rowMeans(select(prostate_chr11_cgi_methylation_2mb_windows, starts_with("T")), na.rm = T) -
  rowMeans(select(prostate_chr11_cgi_methylation_2mb_windows, starts_with("N")), na.rm = T)

plot_cpg_methylation(
  cpg_values = prostate_chr11_cgi_methylation_2mb_windows,
  geom = "geom_linepoint",
  samples = "mean_tumour_normal_difference",
  group_colours = "darkblue",
  title = "Mean Chr11 CpG island Methylation Change in Prostate Tumour Samples"
)

## -----------------------------------------------------------------------------
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gstp1_tss_biomart_result = getBM(mart = mart,
  attributes = c("chromosome_name", "transcription_start_site"),
  filters = "ensembl_transcript_id",
  values = list("ENST00000398606")
)

# Create a GRanges object for the GSTP1 TSS and set chromosome style to UCSC so that it matches the chromosome style in the bedGraphs and methrix object.
gstp1_tss_gr = GenomicRanges::GRanges(seqnames = gstp1_tss_biomart_result$chromosome_name,
  ranges = IRanges(gstp1_tss_biomart_result$transcription_start_site))
seqlevelsStyle(gstp1_tss_gr) = "UCSC"

## -----------------------------------------------------------------------------
gstp1_methylation = extract_cpg_values_from_methrix(
  methrix = prostate_methrix,
  regions = expand_granges(gstp1_tss_gr, 5000, 2000)
)

# Calculate the mean methylation  difference at GSTP1 CpG sites between tumour and normal samples and plot
gstp1_methylation_tumour_normal_difference = dplyr::transmute(gstp1_methylation,
  mean_tumour_normal_difference = rowMeans(dplyr::select(gstp1_methylation, starts_with("T")), na.rm = T) -
    rowMeans(dplyr::select(gstp1_methylation, starts_with("N")), na.rm = T))

## -----------------------------------------------------------------------------
transcript_tpm = read.table(system.file("extdata", "prostate_tumour_normal_transcript_tpm.tsv.gz", package = "methodical"), header = T, sep = "\t", row.names = 1)
gstp1_tpm = unlist(transcript_tpm["ENST00000398606", ])

## -----------------------------------------------------------------------------
gstp1_cpg_methylation_correlation = cpg_feature_cor(feature = gstp1_tpm, cpg_table = gstp1_methylation, method = "s")


# Add distance of CpGs to ENST00000398606 TSS to the table
gstp1_cpg_methylation_correlation$tss_distance = cpg_distances_to_region(query_region = gstp1_tss_gr,
  cpg_names = row.names(gstp1_cpg_methylation_correlation))

# Take a look at the top of the table
head(gstp1_cpg_methylation_correlation)

## -----------------------------------------------------------------------------
plot_cpg_methylation(cpg_values = gstp1_cpg_methylation_correlation, geom = "geom_linepoint", samples = "correlation",
  title = "Correlation of CpG Methylation with GSTP1 Transcription")

## -----------------------------------------------------------------------------
gstp1_tss_cpg_correlation_plot = plot_cpg_methylation(cpg_values = gstp1_cpg_methylation_correlation, geom = "geom_linepoint", samples = "correlation", group_colours = "#66C2A5", reference_region = gstp1_tss_gr, title = "Correlation of CpG Methylation with GSTP1 Transcription", xlabel = "Distance to TSS", ylabel = "Spearman Correlation") + geom_hline(yintercept = 0, linetype = "dashed")
gstp1_tss_cpg_correlation_plot

## -----------------------------------------------------------------------------
gstp1_tss_cpg_methylation_plot = plot_cpg_methylation(cpg_values = gstp1_methylation_tumour_normal_difference, geom = "geom_linepoint", samples = "mean_tumour_normal_difference", group_colours = "#5E4FA2", reference_region = gstp1_tss_gr, title = "CpG Methylation Close to GTSP1 TSS", xlabel = "Distance to TSS") + geom_hline(yintercept = 0, linetype = "dashed")
gstp1_tss_cpg_methylation_plot

## -----------------------------------------------------------------------------
hide_legend(ggplotly(gstp1_tss_cpg_methylation_plot , tooltip = c("position_name")))

## -----------------------------------------------------------------------------
ggarrange(gstp1_tss_cpg_methylation_plot, gstp1_tss_cpg_correlation_plot, nrow = 2, align = "hv")

## -----------------------------------------------------------------------------
methylation_feature_scatter_plot(cpg_values = gstp1_methylation, cpg_name = "chr11:67583741", feature = log(gstp1_tpm),
  title = "GSTP1 Expression and chr11:67583741 Methylation", xlab = "Methylation Value", ylab = "Log GSTP1 TPM",
  add_correlation = T, method = "p")

# We can make the same plot colouring the normal and tumour samples and showing a separate regression line for them
methylation_feature_scatter_plot(
  cpg_values = gstp1_methylation,
  cpg_name = "chr11:67583741",
  feature = log(gstp1_tpm),
  title = "GSTP1 Expression and chr11:67583741 Methylation",
  xlab = "Methylation Value", ylab = "Log GSTP1 TPM",
  regression_line = T, add_correlation = T, method = "p",
  sample_groups = colData(prostate_methrix)$condition
)

# Open as an interactive plot using plotly
ggplotly()

## -----------------------------------------------------------------------------
sessionInfo()

