# Calculate CG and AT methylation for different genomic features in host and GEVE

# Load required packages and functions
library(rtracklayer)
library(methodical)
library(HDF5Array)
library(dplyr)
source("../auxillary_scripts/ggplot_functions.R")

# Load GRanges for annotating chlamy genome
geve_gr = import.bed("../genome_files/GEVE_T2T.bed")
chlamy_genes_gr = readRDS("../genome_files/chlamy_genes_gr.rds")
chlamy_tes_gr = import.gff2("../genome_files/CC2937_T2T.TEs_v3_8.bed.gtf")
centromeres_gr = import.bed("../genome_files/CC2937_T2T.centromeres.bed")
subtelomeres_gr = import.bed("../genome_files/CC2937_T2T.subtelomeres.bed")
mat_sites = readRDS("mAT_sites_whole_genome.rds")
satellites_gr = import.bed("../genome_files/CC2937_T2T.satellites.bed")

# Filter satellites for regions which do not overlap TEs or centromeres
satellites_gr = subsetByOverlaps(satellites_gr, c(chlamy_tes_gr, subtelomeres_gr), invert = T)

# Get a list with the names of different classes of chlamy genes
gene_class_list = readRDS("../genome_files/gene_class_list.rds")

# Define promoters as -400 to +1000 bp
chlamy_tss = resize(chlamy_genes_gr, fix = "start", width = 1)
chlamy_promoters_gr = expand_granges(chlamy_tss, upstream = 400, downstream = 999)

# Separate GEVE genes and promoters into BVP and DVP
bvp_genes = chlamy_genes_gr[gene_class_list$BVP]
dvp_genes = chlamy_genes_gr[gene_class_list$DVP]
bvp_promoters = chlamy_promoters_gr[gene_class_list$BVP]
dvp_promoters = chlamy_promoters_gr[gene_class_list$DVP]

# Get host genes and promoters without GEVE genes
host_genes = chlamy_genes_gr[grep("GEVE", chlamy_genes_gr$gene_id, invert = T)]
host_promoters = chlamy_promoters_gr[grep("GEVE", chlamy_promoters_gr$gene_id, invert = T)]

# Load RSEs for 5mC and 6mA and mask methylation of sites covered by < 10 reads or >= 350 reads
chlamy_cg_meth_rse = loadHDF5SummarizedExperiment("Erazo_2025_5mC_rse")
chlamy_at_meth_rse = loadHDF5SummarizedExperiment("Erazo_2025_6ma_rse")
assay(chlamy_cg_meth_rse, 1)[assay(chlamy_cg_meth_rse, 2) < 10 | assay(chlamy_cg_meth_rse, 2) >= 350] = NA
assay(chlamy_at_meth_rse, 1)[assay(chlamy_at_meth_rse, 2) < 10 | assay(chlamy_at_meth_rse, 2) >= 350] = NA

# Liftover RSEs to T2T assembly
chain = rtracklayer::import.chain("/home/richardheery/chlamydomonas_project/chlamy_genome/t2t/original_contigs_t2t.chain")
rtracklayer::liftOver(SummarizedExperiment::rowRanges(chlamy_cg_meth_rse), chain)
chlamy_cg_meth_rse = liftoverMethRSE(chlamy_cg_meth_rse, chain)
chlamy_at_meth_rse = liftoverMethRSE(chlamy_at_meth_rse, chain)

# Subset RSEs for sites overlapping GEVE and not overlapping GEVE
chlamy_cg_meth_rse_geve = subsetByOverlaps(chlamy_cg_meth_rse, geve_gr)
chlamy_at_meth_rse_geve = subsetByOverlaps(chlamy_at_meth_rse, geve_gr)
chlamy_cg_meth_rse_host = subsetByOverlaps(chlamy_cg_meth_rse, geve_gr, invert = T)
chlamy_at_meth_rse_host = subsetByOverlaps(chlamy_at_meth_rse, geve_gr, invert = T)

# Find mean CpG methylation for different regions in host genome
host_genomic_region_cg_methylation = data.frame(
  "All Sites" = mean(assay(chlamy_cg_meth_rse_host)[, "ThirdRun"], na.rm = T),
  Promoters = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_host, host_promoters))[, "ThirdRun"], na.rm = T),
  Genes = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_host, host_genes))[, "ThirdRun"], na.rm = T),
  TEs = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_host, chlamy_tes_gr))[, "ThirdRun"], na.rm = T),
  Centromeres = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_host, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Subtelomeres = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_host, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Satellites = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_host, satellites_gr))[, "ThirdRun"], na.rm = T),
  check.names = FALSE
)

# Find mean AT methylation for different regions in host genome
host_genomic_region_at_methylation = data.frame(
  "All Sites" = mean(assay(chlamy_at_meth_rse_host)[, "ThirdRun"], na.rm = T),
  Promoters = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_host, host_promoters))[, "ThirdRun"], na.rm = T),
  Genes = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_host, host_genes))[, "ThirdRun"], na.rm = T),
  TEs = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_host, chlamy_tes_gr))[, "ThirdRun"], na.rm = T),
  Centromeres = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_host, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Subtelomeres = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_host, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Satellites = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_host, satellites_gr))[, "ThirdRun"], na.rm = T),
  check.names = FALSE
)

# Find mean CpG methylation for different regions in GEVE
geve_genomic_region_cg_methylation = data.frame(
  "All Sites" = mean(assay(chlamy_cg_meth_rse_geve)[, "ThirdRun"], na.rm = T),
  "BVP\nPromoters" = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, bvp_promoters))[, "ThirdRun"], na.rm = T),
  "BVP\nGenes" = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, bvp_genes))[, "ThirdRun"], na.rm = T),
  "DVP\nPromoters" = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, dvp_promoters))[, "ThirdRun"], na.rm = T),
  "DVP\nGenes" = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, dvp_genes))[, "ThirdRun"], na.rm = T),
  TEs = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, chlamy_tes_gr))[, "ThirdRun"], na.rm = T),
  Centromeres = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Subtelomeres = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Satellites = mean(assay(subsetByOverlaps(chlamy_cg_meth_rse_geve, satellites_gr))[, "ThirdRun"], na.rm = T),
  check.names = FALSE
)

# Find mean AT methylation for different regions in GEVE
geve_genomic_region_at_methylation = data.frame(
  "All Sites" = mean(assay(chlamy_at_meth_rse_geve)[, "ThirdRun"], na.rm = T),
  "BVP\nPromoters" = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, bvp_promoters))[, "ThirdRun"], na.rm = T),
  "BVP\nGenes" = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, bvp_genes))[, "ThirdRun"], na.rm = T),
  "DVP\nPromoters" = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, dvp_promoters))[, "ThirdRun"], na.rm = T),
  "DVP\nGenes" = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, dvp_genes))[, "ThirdRun"], na.rm = T),
  TEs = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, chlamy_tes_gr))[, "ThirdRun"], na.rm = T),
  Centromeres = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Subtelomeres = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, centromeres_gr))[, "ThirdRun"], na.rm = T),
  Satellites = mean(assay(subsetByOverlaps(chlamy_at_meth_rse_geve, satellites_gr))[, "ThirdRun"], na.rm = T),
  check.names = FALSE
)

# Combine tables with CG and AT methylation for host and GEVE
combined_genomic_region_methylation = bind_rows(
  Host_CG  = host_genomic_region_cg_methylation,
  Host_AT  = host_genomic_region_at_methylation,
  GEVE_CG  = geve_genomic_region_cg_methylation,
  GEVE_AT  = geve_genomic_region_at_methylation,
  .id = "group"
)
combined_genomic_region_methylation = tidyr::separate(combined_genomic_region_methylation, group, into = c("context", "site"), sep = "_")

# Convert combined_genomic_region_methylation to long format
combined_genomic_region_methylation_long = tidyr::pivot_longer(combined_genomic_region_methylation, -c(context, site), names_to = "region")

# Set display order of region and contexts 
regions_order = c("All Sites", "Genes", "Promoters", "BVP\nGenes", "BVP\nPromoters", "DVP\nGenes", "DVP\nPromoters",
  "Satellites", "Centromeres", "Subtelomeres", "TEs")
combined_genomic_region_methylation_long$region = factor(combined_genomic_region_methylation_long$region, levels = regions_order)
combined_genomic_region_methylation_long$context = factor(combined_genomic_region_methylation_long$context, c("Host", "GEVE"))

# Remove regions with missing values
combined_genomic_region_methylation_long = filter(combined_genomic_region_methylation_long, !is.na(value))

# Make barplot of mean CG and AT methylation in different regions
combined_genomic_region_methylation_barplot = ggplot(combined_genomic_region_methylation_long,  aes(x = region, y = value*100, fill = site)) +
  geom_col(position = "dodge", color = "black")
combined_genomic_region_methylation_barplot = customize_ggplot_theme(combined_genomic_region_methylation_barplot, 
  title = NULL, ylab = "Mean Methylation %",
  fill_colors = c("#C49E85", "#710627"), facet = "context", facet_nrow = 2, facet_scales = "free_x") +
  theme(strip.background = element_rect(fill = "grey95"), axis.text.x = element_text(size = 12)) 
combined_genomic_region_methylation_barplot
ggsave(plot = combined_genomic_region_methylation_barplot, "plots/region_methylation_barplot.pdf", width = 10, height = 10, device = cairo_pdf)