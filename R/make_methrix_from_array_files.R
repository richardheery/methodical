#' Create a Methrix object from Illumina methylation array files.
#'
#' Needs to use 1.4 times the space used by array_files as temporary space 
#'
#' @param array_files A list of paths to Illumina methylation array files
#' @param probe_name_column The number of the column which corresponds to the name of the probes. Default is 1st column. 
#' @param beta_value_column The number of the column which corresponds to the beta values . Default is 2nd column.  
#' @param methylation_array Name of the Illumina methylation array. At the moment can just be "HumanMethylation450K".
#' @param genome Name of a BSgenome. Allowed values are "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38" and "BSgenome.Hsapiens.NCBI.GRCh38". 
#' @param coldata A data.frame with sample metadata. row.names should be the names of the samples. 
#' @param h5_dir Directory to store the HDF5 methrix object. Cannot already exist.  
#' @param tmpdir Directory in which to create a subdirectory to store temporary files. Defaults to tempdir().
#' @return The methrix object created from array_files
#' @export
make_methrix_from_array_files = function(array_files, probe_name_column = 1, beta_value_column = 2, methylation_array, genome, coldata = NULL, h5_dir, tmpdir = tempdir()){
  
  # Check that genome is one of the allowed BSgenomes
  match.arg(methylation_array, choices = c("HumanMethylation450K"))
  match.arg(genome, choices = c("BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.NCBI.GRCh38"), several.ok = F)
  
  # Extract genome build from provided BSgenome name
  if(genome %in% c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.NCBI.GRCh38")){
    genome_build = "hg38"
  } else {probe_coordinates = genome_build = "hg19"
  }
  
  # Check that h5_dir does not already exist
  if(file.exists(h5_dir)){stop(paste(h5_dir, "already exists"))}
  
  # Get appropriate probe coordinates for methylation array and genome build
  if(methylation_array == "HumanMethylation450K"){
    if(genome_build == "hg38"){
      probe_coordinates = readRDS(system.file("extdata", "illumina_450k_probe_to_cpg_vector_hg38.rds", package = "methodical"))
    } else {probe_coordinates = readRDS(system.file("extdata", "illumina_450k_probe_to_cpg_vector_hg19.rds", package = "methodical"))}
  }
  
  # Set samples names equal to the basenames of the input files if no coldata provided.
  if(is.null(coldata)){
    sample_names = tools::file_path_sans_ext(basename(array_files))
  } else {
    sample_names = row.names(coldata)
  }
  
  # Check that sample_names are unique and that they are equal in length to array_files
  if(nrow(coldata) != length(array_files)){stop("coldata must have one row for each file in array_files")}
  if(any(duplicated(sample_names))){stop("There cannot be duplicate row.names for coldata")}
  
  # Create a directory for intermediate files
  temporary_directory = tempfile(pattern = "make_methrix_from_array_files_", tmpdir = tmpdir)
  dir.create(temporary_directory, recursive = T)
  
  # Create bedGraphs for all files. 
  for(f in seq_along(array_files)){
    
    # Read in file with probe beta values. 
    beta_file = setNames(data.frame(data.table::fread(array_files[f]))[c(probe_name_column, beta_value_column)], c("probe_name", "beta_value"))
    
    # Remove any probes which have missing beta values or whose names were not in the probe_coordinates
    beta_file = dplyr::filter(beta_file, probe_name %in% names(probe_coordinates), !is.na(beta_value))
    
    # Get CpG name from probe using probe coordinates
    beta_file$cpg_name = probe_coordinates[beta_file$probe_name]
    
    # Create a GRanges object for probes with beta value
    beta_file_granges = GenomicRanges::GRanges(beta_file$cpg_name, score = beta_file$beta_value)
    
    # Put the seqlevels in order and sort the GRanges
    GenomeInfoDb::seqlevels(beta_file_granges) = gtools::mixedsort(seqlevels(beta_file_granges))
    beta_file_granges = sort(beta_file_granges)
    
    # Export GRanges as a bedGraph file
    rtracklayer::export.bedGraph(object = beta_file_granges, con = paste0(temporary_directory, "/", sample_names[f], ".bg"))
  }
  
  # Add dummy coverage for the bedGraphs
  methodical::add_dummy_coverage_to_bedgraphs(bedgraphs = list.files(temporary_directory, full.names = T), output_directory = temporary_directory, tmpdir = tmpdir)
  
  # Create ref_cpgs from probe_coordinates
  cpg_df = data.table(data.frame(resize(GRanges(probe_coordinates), width = 2)))
  cpg_df = dplyr::mutate(cpg_df, "chr" = as.character(seqnames), .keep = "unused")
  contig_df = data.table(contig = gtools::mixedsort(unique(cpg_df$chr)), length = lengths(split(cpg_df$chr, cpg_df$chr)[gtools::mixedsort(unique(cpg_df$chr))]))
  ref_cpgs = list(cpgs = cpg_df, contig_lens = contig_df, release_name = "hg38")
  
  # Combine bedGraphs into methrix. TOOK 39 MINUTES
  methrix_object = methrix::read_bedgraphs(list.files(temporary_directory, full.names = T, pattern = "*.bg$"), ref_cpgs = ref_cpgs, chr_idx = 1, start_idx = 2, end_idx = 3, beta_idx = 4, cov_idx = 5, 
    h5 = T, h5_dir = paste(temporary_directory, h5_dir, sep = "/"), coldata = coldata, h5temp = temporary_directory)
  
  # Remove uncovered CpGs
  methrix_object = methrix::remove_uncovered(methrix_object)
  
  # Save updated methrix object. 
  methrix::save_HDF5_methrix(m = methrix_object, dir = h5_dir, replace = F)
  
  # Remove files in temporary directory
  unlink(temporary_directory, recursive = T)
  
  # Return methrix object
  return(methrix_object)
}

