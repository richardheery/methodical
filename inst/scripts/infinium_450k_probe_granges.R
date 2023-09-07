# Make GRanges objects for the probe coordinates of the Illumina 450k array

# Load required packages
library(dplyr)
library(rtracklayer)
library(AnnotationHub)

# Download manifest file for Infinium HumanMethylation450K array from Illumina
download.file(url = "https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv?_gl=1*ocsx4f*_ga*MTk1Nzc4MDkwMy4xNjg3ODcxNjg0*_ga_VVVPY8BDYL*MTY4Nzg3MTY4My4xLjEuMTY4Nzg3MzU5Mi4xMC4wLjA.",
  destfile = "humanmethylation450_15017482_v1-2.csv", method = "wget")

# Load HumanMethylation450 probe information table
infinium_450k_probe_data_hg19 = read.delim("humanmethylation450_15017482_v1-2.csv", header = T, skip = 7, sep = ",")

# Remove any probes that are not for measuring cytosine methylation (SNP probes and controls)
infinium_450k_probe_data_hg19 = filter(infinium_450k_probe_data_hg19, grepl("c", IlmnID))

# Select the infiniumID, chromosome and genomic coordinate (MAPINFO) columns
infinium_450k_probe_data_hg19 = transmute(infinium_450k_probe_data_hg19, IlmnID, CHR = paste0("chr", CHR), MAPINFO)

# Convert into a GRanges object
infinium_450k_probe_granges_hg19 = makeGRangesFromDataFrame(infinium_450k_probe_data_hg19, 
  seqnames.field = "CHR", start.field = "MAPINFO", end.field = "MAPINFO")

# Add probe name to infinium_450k_probe_granges_hg19
infinium_450k_probe_granges_hg19$name = infinium_450k_probe_data_hg19$IlmnID

# Sort infinium_450k_probe_granges_hg19
infinium_450k_probe_granges_hg19 = sort(infinium_450k_probe_granges_hg19)

# Remove probes which do not overlap a CpG site in hg19
probe_seqs_hg19 = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, resize(infinium_450k_probe_granges_hg19, fix = "start", width = 2)))
infinium_450k_probe_granges_hg19 = infinium_450k_probe_granges_hg19[probe_seqs_hg19 == "CG"]

# Create a connection to AnnotationHub
ah = AnnotationHub(localHub = F)

# Load hg19 to hg38 liftover chain from AnnotationHub
hg19_to_hg38_liftover_chain = ah[["AH14150"]]

# Liftover infinium_450k_probe_grange to hg38
infinium_450k_probe_granges_hg38 = liftOver(infinium_450k_probe_granges_hg19, chain = hg19_to_hg38_liftover_chain)

# Remove any probes which coudn't be mapped to hg38
infinium_450k_probe_granges_hg38 = infinium_450k_probe_granges_hg38[lengths(infinium_450k_probe_granges_hg38) > 0]

# Convert infinium_450k_probe_granges_hg38 to a GRanges object
infinium_450k_probe_granges_hg38 = unlist(infinium_450k_probe_granges_hg38)

# Remove probes which do not overlap a CpG site in hg38
probe_seqs_hg38 = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, resize(infinium_450k_probe_granges_hg38, fix = "start", width = 2)))
infinium_450k_probe_granges_hg38 = infinium_450k_probe_granges_hg38[probe_seqs_hg38 == "CG"]

# Put seqlevels in correct order
seqlevels(infinium_450k_probe_granges_hg38) = seqlevels(infinium_450k_probe_granges_hg19)
