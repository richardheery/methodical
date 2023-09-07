# Make GRanges objects with the C position of all CpG sites in hg19 and hg38.
# There are 29.4 million CpGs in hg38 and 28.2 million in hg19

# Load required packages
library(BSgenome.Hsapiens.UCSC.hg19) # package version 1.4.3. UCSC version hg19, based on GRCh37.p13
library(BSgenome.Hsapiens.UCSC.hg38) # package version 1.4.4. UCSC version hg38, based on GRCh38.p14

# Create BSParams object for hg19
params_hg19 = new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19, FUN = matchPattern, exclude = "_")

# Get all CpG site in hg19 as as a list of XStringsViews for standard chromosomes
all_cpgs_hg19  = bsapply(BSParams = params_hg19, pattern = "CG")[extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="all")]

# Get seqinfo for hg19
seqinfo_hg19 = seqinfo(BSgenome.Hsapiens.UCSC.hg19)[extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="all")]

# Convert all_cpgs_hg19 into a GRanges object 
cpg_genome_ranges_hg19 = GRanges(seqnames = rep(names(all_cpgs_hg19),  lengths(all_cpgs_hg19)),
  ranges = IRanges(unlist(lapply(all_cpgs_hg19, start)), unlist(lapply(all_cpgs_hg19, start))), seqinfo = seqinfo_hg19)

# Create BSParams object for hg38
params_hg38 = new("BSParams", X = BSgenome.Hsapiens.UCSC.hg38, FUN = matchPattern, exclude = "_")

# Get all CpG sites in hg38 as a list of XStringsViews for standard chromosomes
all_cpgs_hg38  = bsapply(BSParams = params_hg38, pattern = "CG")[extractSeqlevelsByGroup(species = "Homo_sapiens", style = "UCSC", group = "all")]

# Get seqinfo for hg38
seqinfo_hg38 = seqinfo(BSgenome.Hsapiens.UCSC.hg38)[extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="all")]

# Convert all_cpgs_hg38 into a GRanges object 
cpg_genome_ranges_hg38 = GRanges(seqnames = rep(names(all_cpgs_hg38),  lengths(all_cpgs_hg38)),
  ranges = IRanges(unlist(lapply(all_cpgs_hg38, start)), unlist(lapply(all_cpgs_hg38, start))), seqinfo = seqinfo_hg38)

# Save cpg_genome_ranges_hg19 and cpg_genome_ranges_hg38
saveRDS(cpg_genome_ranges_hg19, "cpg_genome_ranges_hg19.rds")
saveRDS(cpg_genome_ranges_hg38, "cpg_genome_ranges_hg38.rds")
