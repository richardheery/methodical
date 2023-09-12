# Create an annotation file for the human genome combining information for gene bodies, transcription start sites, exons, introns,
# CpG islands, regulatory elements and repetitive elements.
# Gene annotation come from Gencode v38.
# Regulatory elements were downloaded for Ensembl release 109.
# CpG islands come from the UCSC genome browser annotation database (downloaded 26/5/2023).
# Repetitive element information comes from the UCSC genome browser annotation database (downloaded 26/5/2023) and
# was created using RepeatMasker version 4-0-3 with the RepBase library release 20130422. 

# Load required packages
library(dplyr)
library(GenomicRanges)
library(biomaRt)

# Create a directory for storing downloaded annotation files
dir.create("genome_annotation_files")

# Get the names of all standard chromosomes
standard_chroms = extractSeqlevels(species = "Homo_sapiens", style = "UCSC") 

### Get gene annotation from gencode

# Download GTF file for Gencode version 38 
download.file(url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz", 
  destfile = "genome_annotation_files/gencode.v38.annotation.gtf.gz")

# Import Gencode annotation
gencode_annotation = rtracklayer::import.gff2("genome_annotation_files/gencode.v38.annotation.gtf.gz")

# Add a column to gencode_annotation  classifying transcripts
gencode_annotation$region_type = gencode_annotation $gene_type

# rRNAs, rRNA pseudogenes and mitochondrial RNAs are all classed as rRNAs
gencode_annotation$region_type[grepl("rRNA", gencode_annotation$region_type, ignore.case = T)] = "rRNA"

# All remaining pseudogenes classed as pseudogenes
gencode_annotation$region_type[grepl("pseudogene", gencode_annotation$region_type, ignore.case = T)] = "Pseudogene"

# snoRNA, miRNAs, snRNA and misc_RNA are all classed small RNAs
gencode_annotation$region_type[grepl("snorna|mirna|snrna|misc_RNA", gencode_annotation$region_type, ignore.case = T)] = "Small RNA"

# protein_coding is changed to Protein-Coding
gencode_annotation$region_type[gencode_annotation$region_type == "protein_coding"] = "Protein-Coding"

# Filter for protein-coding genes, lncRNAs, pseudogenes, small RNAs and rRNAs
gencode_annotation  = gencode_annotation[gencode_annotation$region_type %in% c("Protein-Coding", "lncRNA", "Pseudogene", "Small RNA", "rRNA")]

# Get transcript annotation from Gencode 
transcript_ranges_hg38 = gencode_annotation[gencode_annotation$type == "transcript"]

# Get transcription start sites for transcripts
tss_sites_hg38 = resize(transcript_ranges_hg38, fix = "start", width = 1, ignore.strand = F)

# Set region_type as TSS
tss_sites_hg38$region_type = paste("TSS", transcript_ranges_hg38$region_type)

# Get exon annotation from Gencode 
exon_ranges_hg38 = gencode_annotation[gencode_annotation$type == "exon"]

# Set region_type as exon
exon_ranges_hg38$region_type = paste("Exon", exon_ranges_hg38$region_type)

# Create lists of transcripts and exons per transcript
transcript_ranges_hg38_list = split(transcript_ranges_hg38, transcript_ranges_hg38$transcript_id)
exon_ranges_hg38_list = split(exon_ranges_hg38, exon_ranges_hg38$transcript_id)

# Filter transcript_ranges_hg38_list for transcripts with exons
transcript_ranges_hg38_list = transcript_ranges_hg38_list[names(exon_ranges_hg38_list)]

# Find introns for each transcript by taking the set difference for the range of the transcript and its exons
# Takes about 5 hours to finish!
intron_ranges_hg38_list = lapply(seq_along(transcript_ranges_hg38_list), function(transcript) 
  GenomicRanges::setdiff(transcript_ranges_hg38_list[[transcript]], exon_ranges_hg38_list[[transcript]]))

# Name intron_ranges_hg38_list with transcripts introns come from
names(intron_ranges_hg38_list) = names(transcript_ranges_hg38_list)
saveRDS(intron_ranges_hg38_list, "genome_annotation_files/intron_ranges_hg38_list.rds")
intron_ranges_hg38_list = readRDS("genome_annotation_files/intron_ranges_hg38_list.rds")

# Convert intron_ranges_hg38_list into a single GRanges and set region_type as intron
intron_ranges_hg38 = unlist(GRangesList(intron_ranges_hg38_list))

# Add names of transcript as metadata column
intron_ranges_hg38$transcript_id = names(intron_ranges_hg38)

# Add region type to intron_ranges_hg38 from transcript_ranges_hg38
intron_ranges_hg38$region_type = transcript_ranges_hg38$region_type[
  match(intron_ranges_hg38$transcript_id, transcript_ranges_hg38$transcript_id)]

# Update region type to include intron
intron_ranges_hg38$region_type = paste("Intron", intron_ranges_hg38$region_type)

# save intron_ranges_hg38
saveRDS(intron_ranges_hg38, "genome_annotation_files/intron_ranges_hg38.rds")

# Remove objects that are no longer needed
rm(gencode_annotation, transcript_ranges_hg38_list, exon_ranges_hg38_list, intron_ranges_hg38_list); gc()

### Get CpG islands from UCSC

# Download a file with masked CpG islands from UCSC
download.file(url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz", 
  destfile = "genome_annotation_files/masked_cpg_islands_hg38_ucsc_hg38.txt.gz")

# Read table of CpG islands. Note that name column is not actually the name, but instead the number of CpGs in the island. 
cpg_islands_hg38 = read.table("genome_annotation_files/masked_cpg_islands_hg38_ucsc_hg38.txt.gz", sep = "\t", header = F)

# Set column names of cpg_islands_hg38 taken from URL:
# http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=cpgIslandExt&hgta_table=cpgIslandExt&hgta_doSchema=describe%20table%20schema
names(cpg_islands_hg38) = c("bin", "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")

# Set row.names of CpG islands to the name of the chromosome and the start site of the island
row.names(cpg_islands_hg38) = paste(cpg_islands_hg38$chrom, cpg_islands_hg38$chromStart, sep = "_")

# Turn CpG islands into a GRanges and sort
cpg_islands_hg38 = sort(makeGRangesFromDataFrame(cpg_islands_hg38, keep.extra.columns = F, ignore.strand = T))

# Give cpg_islands_hg38 a metadata column called region_type with "CpG Island"
cpg_islands_hg38$region_type = "CpG Island"

### Get regulatory regions from biomart  

# Get regulatory feature mart for hg38 for Ensembl version 109
regulatory_mart_hg38 = useEnsembl("ENSEMBL_MART_FUNCGEN", dataset = "hsapiens_regulatory_feature", version = 109)

# Get names of all regulatory features
regulatory_feature_names_hg38 = listFilterOptions(regulatory_mart_hg38, "regulatory_feature_type_name")
regulatory_feature_names_hg38 = setNames(regulatory_feature_names_hg38, regulatory_feature_names_hg38)

# Download GRanges for each type of regulatory feature
regulatory_features_hg38 = lapply(regulatory_feature_names_hg38, function(feature)
  GRanges(
    getBM(
      mart = regulatory_mart_hg38, 
      filters = c("chromosome_name", "regulatory_feature_type_name"), 
      values = list(c(1:22, "X", "Y"),  feature), 
      attributes = c("chromosome_start", "chromosome_end", "chromosome_name", "feature_type_description")
    )
  )
)

# Combine regulatory_features_hg38 into a single GRanges
regulatory_features_hg38 = unlist(GRangesList(regulatory_features_hg38))

# Update seqlevels style to UCSC
seqlevelsStyle(regulatory_features_hg38) = "UCSC"

# Change feature_type_description to region_class
names(mcols(regulatory_features_hg38)) = "region_type"

# Convert remaining regulatory feature names to title case
regulatory_features_hg38$region_type = stringr::str_to_title(regulatory_features_hg38$region_type)

# Change "Ctcf Binding Site" to "CTCF BS" and "Transcription Factor Binding Site" to "TF BS" 
regulatory_features_hg38$region_type[regulatory_features_hg38$region_type == "Ctcf Binding Site"] = "CTCF BS"
regulatory_features_hg38$region_type[regulatory_features_hg38$region_type == "Transcription Factor Binding Site"] = "TF BS"

### Get repetitive element annotation from UCSC

# Download a file with repetitive element annotation from UCSC
options(timeout = max(300, getOption("timeout")))
download.file(url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz", 
  destfile = "genome_annotation_files/ucsc_repetitive_sequence_annotation.tsv.gz")

# Read in table for repetitive element annotation
repetitive_element_table = data.table::fread("genome_annotation_files/ucsc_repetitive_sequence_annotation.tsv.gz", 
  header = F, fill = T, skip = 3)

# Add names to the relevant columns
names(repetitive_element_table) = c(rep(NA, 4), "chromosome", "start", "end", NA, NA, NA, "type", rep(NA, 4))

# Extract class and family from type
repetitive_element_table$class = gsub("/.*", "", repetitive_element_table$type)
repetitive_element_table$family = ifelse(grepl("/", repetitive_element_table$type), gsub(".*/", "", repetitive_element_table$type), NA)

# Discard unnecessary columns
repetitive_element_table = dplyr::select(repetitive_element_table, chromosome, start, end, class, family)

# Remove ranges where the class or family is uncertain, denoted by the presence of a ? in the class/family name. 
repetitive_element_table = filter(repetitive_element_table, !grepl("\\?", repetitive_element_table$class) & !grepl("\\?", repetitive_element_table$family))

# Convert the table into a GRanges object
repeat_ranges_hg38 = GRanges(repetitive_element_table)

# Add a region_type metadata column, initialized to NA
repeat_ranges_hg38$region_type = NA

# Add the most abundant groups of repeat element to region_type
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$class == "DNA")] = "DNA Transposon"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$class == "LTR")] = "LTR"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$class == "Satellite")] = "Satellite"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$class == "Simple_repeat")] = "Simple Repeat"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$class == "Low_complexity")] = "Low Complexity"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$family == "Alu")] = "Alu"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$family == "MIR")] = "MIR"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$family == "SVA")] = "SVA"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$family == "L1")] = "L1"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$family == "L2")] = "L2"
repeat_ranges_hg38$region_type[which(repeat_ranges_hg38$family == "CR1")] = "CR1"

# Remove repeats which do not belong to the above classes
repeat_ranges_hg38 = repeat_ranges_hg38[!is.na(repeat_ranges_hg38$region_type)]

# Remove repetitive_element_table
rm(repetitive_element_table); gc()

### Combine all annotation sources

# Create a GRanges combining all the other annotation GRanges. 
# 82% of the genome is annotated by genome_annotation_hg38. 
genome_annotation_hg38 = c(tss_sites_hg38, exon_ranges_hg38, intron_ranges_hg38, 
  cpg_islands_hg38, regulatory_features_hg38, repeat_ranges_hg38)

# Remove all ranges not on standard chromosomes and update seqlevels
seqlevels(genome_annotation_hg38, pruning.mode = "coarse") = standard_chroms

# Add seqinfo to genome_annotation_hg38 from BSgenome.Hsapiens.UCSC.hg38
seqinfo(genome_annotation_hg38) = seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[standard_chroms]

# Sort genome_annotation_hg38
genome_annotation_hg38 = sort(genome_annotation_hg38, ignore.strand = T)

# Remove ranges names and all metadata columns except region_type
mcols(genome_annotation_hg38) = mcols(genome_annotation_hg38)["region_type"]
names(genome_annotation_hg38) = NULL

# Save genome_annotation_hg38
saveRDS(genome_annotation_hg38, "genome_annotation_files/genome_annotation_hg38.rds")
