#' Create a HDF5-backed RangedSummarizedExperiment for methylation values in meth_files
#'
#' @param meth_files A vector of paths to input methylation files. All sites in each file are assumed to be for the 
#' same sequence context e.g. CG or CHG. Automatically detects if meth_files contain a header if every field in the first line is a character. 
#' @param seqnames_col The column number in meth_files which corresponds to the sequence names. 
#' @param start_col The column number in meth_files which corresponds to the genomic start coordinate. 
#' @param total_reads_col The column number in meth_files which corresponds to the total number of reads for the position. 
#' @param meth_reads_col The column number in meth_files which corresponds to the number of methylated reads for the position.
#' @param unmeth_reads_col The column number in meth_files which corresponds to the number of unmethylated reads for the position.
#' @param meth_fraction_col The column number in meth_files which corresponds to the fraction of reads that support methylation at the position.
#' Will be converted to a proportion if it appears to be a fraction.
#' @param zero_based TRUE or FALSE indicating if files are zero-based. Default value is TRUE. 
#' @param meth_sites A GRanges object with non-overlapping locations of methylation sites of interest e.g. CpG sites. 
#' Any methylation sites in meth_files that are not in meth_sites are ignored.
#' @param sequence_context A single character string or DNAString with the 
#' sequence context of the methylation sites e.g. CG or CHG. If a character, must be coercible to a DNAString. Default is "CG".
#' @param collapse_strands TRUE or FALSE indicating whether or not to collapse data on + and - strands. 
#' Only makes sense for symmetrically methylated contexts e.g. CG or CHG and meth_sites should include ranges on both the + and - strands if TRUE. 
#' @param decimal_places Optional integer indicating the number of decimal places to round beta values to. Default is not to round.
#' @param sample_metadata Sample metadata to be used as colData for the RangedSummarizedExperiment.
#' @param hdf5_dir Directory to save HDF5 file. Is created if it doesn't exist. HDF5 file is called assays.h5. 
#' @param overwrite TRUE or FALSE indicating whether to allow overwriting if hdf5_dir already exists. Default is FALSE.
#' @param chunkdim The dimensions of the chunks for the HDF5 file. Should be a vector of length 2 giving the number of rows and then the number of columns in each chunk.
#' Uses HDF5Array::getHDF5DumpChunkDim(length(meth_sites), length(meth_files))) by default. 
#' @param temporary_dir Name to give temporary directory created to store intermediate files. A directory with this name cannot already exist. 
#' Default is to create a subdirectory named "temporary_meth_chunks_" inside the directory given by `tempdir()`. Will be deleted after completion. 
#' @param BPPARAM A BiocParallelParam object for parallel processing. Defaults to `BiocParallel::SerialParam()`. 
#' @param ... Additional arguments to be passed to HDF5Array::HDF5RealizationSink() for controlling the physical properties of the created HDF5 file, 
#' such as compression level. Uses the defaults for any properties that are not specified. 
#' @return A RangedSummarizedExperiment with two assays for all methylation sites in meth_sites, 
#' beta with the proportion of methylated reads and Cov with the total  number of reads for each site. methylation sites will be in the same order as sort(meth_sites). 
#' @export
#' @examples
#' 
#' # Load CpGs within first million base pairs of chromosome 1 as a GRanges object
#' data("hg38_cpgs_subset", package = "methodical")
#' 
#' # Get paths to meth_files
#' meth_files <- list.files(path = system.file('extdata', package = 'methodical'), 
#'   pattern = ".bg.gz", full.names = TRUE)
#' 
#' # Create sample metadata
#' sample_metadata <- data.frame(
#'   tcga_project = gsub("_.*", "", gsub("TCGA_", "", basename(meth_files))),
#'   sample_type = ifelse(grepl("N", basename(meth_files)), "Normal", "Tumour"),
#'   row.names = tools::file_path_sans_ext(basename(meth_files))
#' )
#' 
#' # Create a HDF5-backed RangedSummarizedExperiment from meth_files
#' meth_rse <- makeMethRSEFromInputFiles(meth_files = meth_files, 
#'   seqnames_col = 1, start_col = 2, meth_reads_col = 4, unmeth_reads_col = 5, 
#'   zero_based = TRUE, meth_sites = hg38_cpgs_subset, sample_metadata = sample_metadata, 
#'   hdf5_dir = paste0(tempdir(), "/bedgraph_hdf5_1"))
#'   
makeMethRSEFromInputFiles <- function(meth_files, seqnames_col, start_col, 
  total_reads_col = NULL, meth_reads_col = NULL, unmeth_reads_col = NULL, meth_fraction_col = NULL, 
  zero_based, meth_sites, sequence_context = "CG", collapse_strands = TRUE, decimal_places = NA, sample_metadata = NULL, 
  hdf5_dir, overwrite = FALSE, chunkdim = NULL, temporary_dir = NULL, BPPARAM = BiocParallel::SerialParam(), ...){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_files, "character"), 
    is.numeric(seqnames_col) && length(seqnames_col) == 1 && seqnames_col > 0 && seqnames_col %% 1 == 0,
    is.numeric(start_col) && length(start_col) == 1 && start_col > 0 && start_col %% 1 == 0,
    is.null(total_reads_col) | is.numeric(total_reads_col) && length(total_reads_col) == 1 && total_reads_col > 0 && total_reads_col %% 1 == 0,
    is.null(meth_reads_col) | is.numeric(meth_reads_col) && length(meth_reads_col) == 1 && meth_reads_col > 0 && meth_reads_col %% 1 == 0,
    is.null(unmeth_reads_col) | is.numeric(unmeth_reads_col) && length(unmeth_reads_col) == 1 && unmeth_reads_col > 0 && unmeth_reads_col %% 1 == 0,
    is.null(meth_fraction_col) | is.numeric(meth_fraction_col) && length(meth_fraction_col) == 1 && meth_fraction_col > 0 && meth_fraction_col %% 1 == 0,
    S4Vectors::isTRUEorFALSE(zero_based), is(meth_sites, "GRanges"), 
    (is(sequence_context, "DNAString") || is(sequence_context, "character")) && nchar(sequence_context) > 0,
    S4Vectors::isTRUEorFALSE(collapse_strands), is(decimal_places, "numeric") | is.na(decimal_places), 
    is(sample_metadata, "data.frame") | is.null(sample_metadata), is(hdf5_dir, "character"),
    is(temporary_dir, "character") | is.null(temporary_dir), is(BPPARAM, "BiocParallelParam"))
  
  # Check that if sequence_context is a character that it is coercible to a DNAString
  tryCatch({invisible(Biostrings::DNAString(sequence_context))}, 
    error = function(e) stop(paste(sequence_context, "is not coercible to a DNAString")))
    
  # Check that total_reads_col and at least one of meth_fraction_col, meth_reads_col or unmeth_reads_col or
  # else both meth_reads_col and unmeth_reads_col are provided 
  if((is.null(total_reads_col) | all(sapply(list(meth_reads_col, unmeth_reads_col, meth_fraction_col), is.null))) &&
    (is.null(meth_reads_col) | is.null(unmeth_reads_col))){
      stop("Either both meth_reads_col and unmeth_reads_col should be provided or else total_reads_col and one other column")
  }
  
  # Check that different columns given for meth_files columns
  if(anyDuplicated(c(seqnames_col, start_col, total_reads_col, meth_reads_col, unmeth_reads_col, meth_fraction_col))){
    stop("Duplicate column indices given for seqnames, start, total reads, meth reads, unmeth reads or meth fraction")
  }
  
  # Check that all ranges in meth_sites have a width of 1, that they are disjoint 
  # and that they are stranded if collapse_strands is TRUE
  if(any(GenomicRanges::width(meth_sites) != 1)){
    stop("All meth_sites should have a width of 1")
  }
  if(!GenomicRanges::isDisjoint(meth_sites)){
    stop("There cannot be overlapping regions in meth_sites")
  }
  if(collapse_strands && "*" %in% strand(meth_sites)){
    stop("If collapse_strands is TRUE, all ranges in meth_sites must be stranded (on + or - strand")
  }
    
  # Check if meth_sites is sorted and print a message if it is not. 
  if(!all(meth_sites == sort(meth_sites, ignore.strand = T))){
    message("meth_sites is not sorted (ignoring strand). It will be sorted and this sorted order used for methylation sites in the HDF5 file")
  }
  
  # Indicate the final meth_sites to use based on whether collapse_strands is TRUE
  if(collapse_strands){
    meth_sites_final <- meth_sites[strand(meth_sites) != "-"] 
    strand(meth_sites_final) <- "*"
  } else {
    meth_sites_final <- meth_sites
  }
  
  # Convert meth_sites into a data.table
  meth_sites_df <- data.table::data.table(data.frame(meth_sites)[c("seqnames", "start", "strand")])
  
  # Set seqnames and start as keys for meth_sites_df and convert back to a sorted GRanges 
  data.table::setkey(meth_sites_df, seqnames, start)
  meth_sites <- GenomicRanges::makeGRangesFromDataFrame(meth_sites_df, end.field = "start")
  
  # If temporary_dir not provided, set it to a directory in tempdir()
  if(is.null(temporary_dir)){
    temporary_dir <- tempfile("temporary_meth_chunks_")
  }
  
  # Check temporary directory doesn't already exist and create it if it doesn't
  if(dir.exists(temporary_dir)){
    stop(paste("Directory", temporary_dir, "already exists. Please provide a temporary directory name that isn't already in use."))
  } else {
    dir.create(temporary_dir)
  }
  
  # Perform setup
  setup <- .make_meth_rse_setup(meth_files = meth_files, meth_sites = meth_sites_final, sample_metadata = sample_metadata, 
    hdf5_dir = hdf5_dir, overwrite = overwrite, chunkdim = chunkdim, 
    temporary_dir = temporary_dir, ...)
  
  # Create a list which specifies the columns in meth_files
  meth_files_columns <- list(seqnames_col = seqnames_col, start_col = start_col, total_reads_col = total_reads_col,
    meth_reads_col = meth_reads_col, unmeth_reads_col = unmeth_reads_col, meth_fraction_col = meth_fraction_col)
  
  # Read in meth_files and write data from chunks to appropriate temporary directory
  .split_meth_files_into_chunks(meth_files = meth_files, meth_files_columns,
    file_grid_columns = setup$file_grid_columns, meth_sites_df = meth_sites_df, collapse_strands = collapse_strands,
    sequence_context = sequence_context, meth_site_groups = setup$meth_site_groups, temp_chunk_dirs = setup$temp_chunk_dirs, 
    zero_based = zero_based, decimal_places = decimal_places, BPPARAM = BPPARAM)
  
  # Write the chunks to the HDF5 file
  .write_chunks_to_hdf5(temp_chunk_dirs = setup$temp_chunk_dirs, files_in_chunks = setup$files_in_chunks, 
    beta_sink = setup$beta_sink, Cov_sink = setup$Cov_sink, hdf5_grid = setup$hdf5_grid)
  
  # Create a RangedSummarizedExperiment
  rse <- .create_meth_rse_from_hdf5(hdf5_filepath = setup$hdf5_filepath, hdf5_dir = hdf5_dir,
    meth_sites = meth_sites_final, sample_metadata = sample_metadata)
  
  return(rse)
  
}

#' Create a HDF5-backed RangedSummarizedExperiment for methylation values in array files
#'
#' @param array_files A vector of paths to input files. Automatically detects if array_files contain a header if every field in the first line is a character. 
#' @param probe_name_column The number of the column which corresponds to the name of the probes. Default is 1st column. 
#' @param beta_value_column The number of the column which corresponds to the beta values . Default is 2nd column.  
#' @param normalization_factor An optional numerical value to divide methylation values by to convert them to fractions e.g. 100 if they are percentages. 
#' Default is not to leave values as they are in the input files. 
#' @param decimal_places Integer indicating the number of decimal places to round beta values to. Default is 2. 
#' @param probe_ranges A GRanges object giving the genomic locations of probes where each region corresponds to a separate probe. 
#' There should be a metadata column called name with the name of the probe associated with each region. 
#' Any probes in array_files that are not in probe_ranges are ignored. 
#' @param sample_metadata Sample metadata to be used as colData for the RangedSummarizedExperiment
#' @param hdf5_dir Directory to save HDF5 file. Is created if it doesn't exist. HDF5 file is called assays.h5. 
#' @param dataset_name Name to give data set in HDF5 file. Default is "beta".
#' @param overwrite TRUE or FALSE indicating whether to allow overwriting if dataset_name already exists in assays.h5. Default is FALSE.
#' @param chunkdim The dimensions of the chunks for the HDF5 file. Should be a vector of length 2 giving the number of rows and then the number of columns in each chunk.
#' @param temporary_dir Name to give a temporary directory to store intermediate files. A directory with this name cannot already exist. 
#' Default is to create a name using tempfile("temporary_meth_chunks_"). 
#' @param BPPARAM A BiocParallelParam object for parallel processing. Defaults to `BiocParallel::SerialParam()`. 
#' @param ... Additional arguments to be passed to HDF5Array::HDF5RealizationSink() for controlling the physical properties of the created HDF5 file, 
#' such as compression level. Uses the defaults for any properties that are not specified. 
#' @return A RangedSummarizedExperiment with methylation values for all methylation sites in meth_sites. Methylation sites will be in the same order as sort(meth_sites). 
#' @export
#' @examples
#' # Get human CpG sites for hg38 genome build
#' data("infinium_450k_probe_granges_hg19", package = "methodical")
#' 
#' # Get paths to array files
#' array_files <- list.files(path = system.file('extdata', package = 'methodical'), 
#'   pattern = ".txt.gz", full.names = TRUE)
#' 
#' # Create sample metadata
#' sample_metadata <- data.frame(
#'   tcga_project = "LUAD",
#'   sample_type = "Tumour", submitter = gsub("_01.tsv.gz", "", basename(array_files)),
#'   row.names = gsub(".tsv.gz", "", basename(array_files))
#' )
#' 
#' # Create a HDF5-backed RangedSummarizedExperiment from array files using default chumk dimensions
#' meth_rse <- makeMethRSEFromArrayFiles(array_files = array_files, 
#'  probe_ranges = infinium_450k_probe_granges_hg19, 
#'  sample_metadata = sample_metadata, hdf5_dir =  paste0(tempdir(), "/array_file_hdf5_1"))
#'
makeMethRSEFromArrayFiles <- function(array_files, probe_name_column = 1, beta_value_column = 2, 
  normalization_factor = NULL, decimal_places = NA, probe_ranges, sample_metadata = NULL, hdf5_dir, dataset_name = "beta", 
  overwrite = FALSE, chunkdim = NULL, temporary_dir = NULL, BPPARAM = BiocParallel::SerialParam(), ...){
  
  # Check that inputs have the correct data type
  stopifnot(is(array_files, "character"), is(probe_name_column, "numeric") & probe_name_column >= 1,
    is(beta_value_column, "numeric") & beta_value_column >= 1, 
    is(normalization_factor, "numeric") | is.null(normalization_factor),
    is(decimal_places, "numeric") | is.na(decimal_places), is(probe_ranges, "GRanges"),
    is(sample_metadata, "data.frame") | is.null(sample_metadata), is(hdf5_dir, "character"),
    is(dataset_name, "character"), S4Vectors::isTRUEorFALSE(overwrite), is(chunkdim, "numeric") | is.null(chunkdim),
    is(temporary_dir, "character") | is.null(temporary_dir), is(BPPARAM, "BiocParallelParam"))
  
  # Check that normalization_factor is a whole integer if provided
  if(!is.null(normalization_factor)){
    if(length(normalization_factor) != 1 | normalization_factor %% 1 != 0 | normalization_factor < 0){
      stop("normalization_factor should be single whole number")
    }
  }
  
  # Check that probe_ranges has a metadata column called name and that there are no duplicate names
  if(!"name" %in% names(mcols(probe_ranges))){
    stop("probe_ranges must have a metadata column called name")
  } else if(anyDuplicated(probe_ranges$name)){
    stop("probe_ranges$name cannot contain any duplicates")
  }
  
  # If temporary_dir not provided, set it to a directory in tempdir()
  if(is.null(temporary_dir)){
    temporary_dir <- tempfile("temporary_meth_chunks_")
  }
  
  # Check temporary directory doesn't already exist and create it if it doesn't
  if(dir.exists(temporary_dir)){
    stop(paste("Directory", temporary_dir, "already exists. Please provide a temporary directory name that isn't already in use."))
  } else {
    dir.create(temporary_dir)
  }
  
  # Check if probe_ranges is sorted and print a message if it is not. 
  if(!S4Vectors::isSorted(probe_ranges)){
    message("probe_ranges is not sorted. It will be sorted and this sorted order will be used for methylation sites in the HDF5 file")
  }
  
  # Perform setup
  setup <- .make_meth_rse_setup(meth_files = array_files, meth_sites = probe_ranges, sample_metadata = sample_metadata, 
    hdf5_dir = hdf5_dir, dataset_name = dataset_name, overwrite = overwrite, chunkdim = chunkdim, 
    temporary_dir = temporary_dir, ...)
  
  # Read in array files and write data from chunks to appropriate temporary directory
  probe_sites_df <- .split_meth_array_files_into_chunks(array_files = array_files, probe_name_column = probe_name_column, 
    beta_value_column = beta_value_column, file_grid_columns = setup$file_grid_columns, probe_ranges = probe_ranges,
    probe_groups = setup$meth_site_groups, temp_chunk_dirs = setup$temp_chunk_dirs, 
    normalization_factor = normalization_factor, decimal_places = decimal_places, BPPARAM = BPPARAM)
  
  # Write the chunks to the HDF5 file
  .write_chunks_to_hdf5(hdf5_sink = setup$hdf5_sink, hdf5_grid = setup$hdf5_grid, 
    temp_chunk_dirs = setup$temp_chunk_dirs, files_in_chunks = setup$files_in_chunks)
  
  # Create a RangedSummarizedExperiment
  rse <- .create_meth_rse_from_hdf5(hdf5_filepath = setup$hdf5_filepath, hdf5_dir = hdf5_dir,
    meth_sites_df = probe_sites_df, sample_metadata = sample_metadata)
  
  # Delete temporary_dir if it is empty
  if(length(list.files(temporary_dir)) == 0){
    unlink(temporary_dir, recursive = TRUE)
  }
  
  return(rse)
  
}

#' infinium_450k_probe_granges_hg19
#'
#' The hg19 genomic coordinates for methylation sites analysed by the Infinium HumanMethylation450K array.
#'
#'@format GRanges object with 482,421 ranges and one metadata column name giving the name of the associated probe. 
#'@source Derived from the manifest file downloaded from https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv?_gl<-1*ocsx4f*_ga*MTk1Nzc4MDkwMy4xNjg3ODcxNjg0*_ga_VVVPY8BDYL*MTY4Nzg3MTY4My4xLjEuMTY4Nzg3MzU5Mi4xMC4wLjA.
"infinium_450k_probe_granges_hg19"

#' Convert a Methrix object into a RangedSummarizedExperiment
#'
#' @param methrix A methrix object
#' @param assays A vector indicating the names of assays in methrix used to create a RangedSummarizedExperiment. Can be one or both of "beta" and "cov". 
#' Default is both "beta" and "cov" assays. 
#' @return A RangedSummarizedExperiment 
#' @export
#' @examples
#' # Load a sample methrix object
#' data("methrix_data", package = "methrix")
#'   
#' # Convert methrix to a RangedSummarizedExperiment with one assay for the methylation beta values
#' meth_rse <- methodical::methrixToRSE(methrix_data, assays = "beta")
#' print(meth_rse)
methrixToRSE <- function(methrix, assays = c("beta", "cov")){
  
  # Check that inputs have the correct data type
  stopifnot(is(methrix, "methrix"), is(assays, "character"))
  
  # Check that allowed values are provided for assays
  if(any(!assays %in% c("beta", "cov"))){stop("assays should only be one or both of \"beta\" and \"cov\"")}
  
  # Extract rowdata from methrix
  rowdata <- SummarizedExperiment::rowData(methrix)
  
  # Extract GRanges of methylation sites from methrix
  methrix_ranges <- GenomicRanges::makeGRangesFromDataFrame(rowdata, start.field = "start", end.field = "start")
  
  # Create a RangedSummarizedExperiment from methrix
  rse <- SummarizedExperiment::SummarizedExperiment(
    assays = SummarizedExperiment::assays(methrix)[assays], 
    colData = SummarizedExperiment::colData(methrix), 
    rowRanges = methrix_ranges)
  
  # Sort rse and return
  rse <- sort(rse)
  return(rse)
}