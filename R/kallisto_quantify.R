#' Create an index file for running Kallisto
#'
#' @param path_to_kallisto Path to kallisto executable
#' @param transcripts_fasta Path to a fasta file for the transcripts to be quantified.
#' @param index_name Name to give the created index file. Default is "kallisto_index.idx".
#' @return Invisibly returns TRUE. 
#' @export
#' @examples \dontrun{
#' # Download transcripts FASTA from Gencode
#' download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz")
#' 
#' # Locate the kallisto executable (provided that it is on the path)
#' kallisto_path <- system2("which", args = "kallisto", stdout = TRUE)
#' 
#' # Create transcripts index for use with Kallisto
#' methodical::kallisto_index(kallisto_path, transcripts_fasta = "gencode.v44.transcripts.fa.gz")
#' }
kallisto_index <- function(path_to_kallisto, transcripts_fasta, index_name = "kallisto_index.idx"){
  
  # Check that inputs have the correct data type
  stopifnot(is(path_to_kallisto, "character"), is(transcripts_fasta, "character"),
    is(index_name, "character"))
  
  # Get the canonical path to kallisto
  path_to_kallisto <- normalizePath(path_to_kallisto)
  
  # Check if kallisto can be executed from the given path
  if(suppressWarnings(system2(command = path_to_kallisto, args = "version", stdout = NULL, stderr = NULL)) != 0){
    stop("kallisto could can not be executed from given path")
  }
  
  # Create the index
  system2(command = path_to_kallisto,
    args = paste("index -i", index_name, transcripts_fasta))
  
  # Invisibly return TRUE
  invisible(return(TRUE))
  
}

#' Run kallisto on a set of FASTQ files and merge the results
#'
#' @param path_to_kallisto Path to kallisto executable
#' @param kallisto_index Path to a kallisto index 
#' @param forward_fastq_files A vector with the paths to forward FASTQ files. Each file should correspond to the file at the same position in reverse_fastq_files.
#' @param reverse_fastq_files A vector with the paths to reverse FASTQ files. Each file should correspond to the file at the same position in forward_fastq_files.
#' @param sample_names A vector with the names of samples for each pair of samples from forward_fastq_files and reverse_fastq_files
#' @param output_directory The name of the directory to save results in. Will be created if it doesn't exist. 
#' @param merged_output_prefix Prefix to use for names of merged output files for counts and TPM which take the form 
#' {merged_output_prefix}_counts_merged.tsv.gz and {merged_output_prefix}_tpm_merged.tsv.gz. Default prefix is "kallisto_transcript" i.e. default output merged output files are 
#' kallisto_transcript_counts_merged.tsv.gz and kallisto_transcript_tpm_merged.tsv.gz.
#' @param messages_file Name of file to save kallisto run messages. If no file name given, information is printed to stdout.
#' @param ncores The number of cores to use. Default is 1.
#' @param number_bootstraps The number of bootstrap samples. Default is 100. 
#' @return The path to the merged counts table. 
#' @export
kallisto_quantify <- function(path_to_kallisto, kallisto_index, forward_fastq_files, reverse_fastq_files, 
  sample_names, output_directory, merged_output_prefix = "kallisto_transcript", messages_file = "", ncores = 1, number_bootstraps  = 100){
  
  # Check that inputs have the correct data type
  stopifnot(is(path_to_kallisto, "character"), is(kallisto_index, "character"), 
    is(forward_fastq_files, "character"), is(reverse_fastq_files, "character"),
    is(sample_names, "character"), is(output_directory, "character"),
    is(merged_output_prefix, "character"), is(messages_file, "character"),
    is(ncores, "numeric") & ncores >= 1, is(number_bootstraps, "numeric") & number_bootstraps >= 1)
  
  # Get the canonical path to kallisto
  path_to_kallisto <- normalizePath(path_to_kallisto)
  
  # Check if kallisto can be executed from the given path
  if(suppressWarnings(system2(command = path_to_kallisto, args = "version", stdout = NULL, stderr = NULL)) != 0){
    stop("kallisto could can not be executed from given path")
  }
  
  # Check that path to kallisto_index exits
  if(!file.exists(kallisto_index)){stop("kallisto_index could not be found")}
  
  # Check that at least one pair of FASTQs provided
  if(length(forward_fastq_files) == 0 | length(forward_fastq_files) == 0){stop(
    "At least one pair of forward and reverse FASTQ files must be provided")
  }
  
  # Check that forward_fastq_files, reverse_fastq_files and sample names have the same lengths
  if(length(forward_fastq_files) != length(reverse_fastq_files)){
    stop("forward_fastq_files and reverse_fastq_files have different lengths")
  }
  if(length(forward_fastq_files) != length(sample_names)){
    stop("sample_names does not have same length as forward_fastq_files and reverse_fastq_files")
  }
  
  # Check if messages_file already exists
  if(messages_file != ""){if(file.exists(messages_file)){stop(paste(messages_file, "already exists"))}}
  
  # Create output_directory if it doesn't exist
  if(!dir.exists(output_directory)){dir.create(output_directory)}
  
  # Create output subdirectories for each pair of FASTQ samples by pasting sample_names to output_directory
  sample_directories <- paste(output_directory, sample_names, sep = "/")
  
  # Check that sample directories do not already exist
  for(directory in sample_directories) {
    if(file.exists(directory)){stop(paste(directory, "already exists"))}
  }
  
  # Create names for merged output files and check that they do not already exist
  merged_counts_file <- paste(merged_output_prefix, "counts_merged.tsv.gz", sep = "_")
  merged_tpm_file <- paste(merged_output_prefix, "tpm_merged.tsv.gz", sep = "_")
  
  # Check if merged_counts_file already exists
  if(file.exists(paste(output_directory, merged_counts_file, sep = "/"))){
    stop(paste(paste(output_directory, merged_counts_file, sep = "/"), "already exists"))}
  
  # Check if merged_tpm_file already exists
  if(file.exists(paste(output_directory, merged_tpm_file, sep = "/"))){
    stop(paste(paste(output_directory, merged_tpm_file, sep = "/"), "already exists"))}
  
  # Loop through each pair of FASTQs and run kallisto
  for(pair in seq_along(forward_fastq_files)){
    message(paste("Starting to process FASTQ pair", pair))
    system2(command = path_to_kallisto,
      args = paste(
        "quant -i", kallisto_index, "-t", ncores, "-b", number_bootstraps, "-o", 
        sample_directories[pair], forward_fastq_files[pair], reverse_fastq_files[pair]
        ),
      stderr = messages_file)
  }
  
  # Get paths to all abundance files
  abundance_files <- paste0(sample_directories, "/abundance.tsv")
  
  # Create a data.frame with the counts calculated using kallisto for each sample
  kallisto_counts <- data.frame(setNames(lapply(abundance_files, function(x) 
    data.table::fread(x, sep = "\t", header = TRUE)$est_counts), basename(sample_directories)))
  
  # Create a data.frame with the TPM values calculated using kallisto for each sample
  kallisto_tpm <- data.frame(setNames(lapply(abundance_files, function(x) 
    data.table::fread(x, sep = "\t", header = TRUE)$tpm), basename(sample_directories)))
  
  # Get names of transcripts and add to output tables
  transcript_names <- data.table::fread(abundance_files[1], sep = "\t", header = TRUE)$target_id
  row.names(kallisto_counts) <- transcript_names
  kallisto_counts <- tibble::rownames_to_column(kallisto_counts, "transcript_id")
  row.names(kallisto_tpm) <- transcript_names
  kallisto_tpm <- tibble::rownames_to_column(kallisto_tpm, "transcript_id")
      
  # Write tables to output directory
  data.table::fwrite(kallisto_counts, paste(output_directory, merged_counts_file, sep = "/"), sep = "\t", row.names = FALSE, quote = FALSE)
  data.table::fwrite(kallisto_tpm, paste(output_directory, merged_tpm_file, sep = "/"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Return path to counts table
  return(paste(output_directory, merged_counts_file, sep = "/"))
  
}

#' Combine the expression values of transcripts to get overall expression of their associated genes
#' 
#' @param transcript_expression_table A table where rows are transcripts and columns are samples. Row names should be the names of transcripts. 
#' @param gene_to_transcript_list A list of vectors where the name of each list entry is a gene name and its elements are the names of transcripts.
#' Can alternatively be a GRangeList where the name of each list element is a gene and the names of the individual ranges are transcripts.
#' @return A data.frame with the sum of transcript expression values for genes where rows are genes and columns are samples
#' @export
sum_transcript_values_for_genes <- function(transcript_expression_table, gene_to_transcript_list){
  
  # Check that inputs have the correct data type
  stopifnot(is(transcript_expression_table, "data.frame") | is(transcript_expression_table, "matrix"), 
    is(gene_to_transcript_list, "list"), all(sapply(gene_to_transcript_list, class) = "character"))
  
  # If gene_to_transcript_list is a GRangesList, extract a list of vectors matching genes to transcripts
  if(is(gene_to_transcript_list, "GRangesList")){
    gene_to_transcript_list <- lapply(gene_to_transcript_list, names)
  }
  
  # Get the sum of the expression values for all transcripts associated with each gene in each sample
  `%do%` <- foreach::`%do%`
  results_list <- foreach::foreach(gene_transcripts = gene_to_transcript_list) %do% {
    
    # Sum the transcript values for each transcript associated with a gene
    gene_results <- colSums(transcript_expression_table[gene_transcripts, ], na.rm = TRUE)
    
    # Samples where all values for gene_transcripts are NA are given a value of NA for the gene. 
    # This is done as colSums returns a value of 0 if all values in a column are NA and na.rm <- TRUE.  
    gene_results[apply(transcript_expression_table[gene_transcripts, ], 2, function(x) all(is.na(x)))] <- NA
    gene_results
  }
  
  # Set names for results_list
  names(results_list) <- names(gene_to_transcript_list)
  
  # Turn results_list into a data.frame
  results_table <- data.frame(dplyr::bind_rows(results_list, .id = "gene_name"))
  
  # Set gene names as row.names and return
  return(tibble::column_to_rownames(results_table, "gene_name"))
}
