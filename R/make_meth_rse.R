#' Create a HDF5-backed RangedSummarizedExperiment for methylation values in bedGraphs
#'
#' @param bedgraphs A vector of paths to bedGraph files. Automatically detects if bedGraphs contain a header if every field in the first line is a character. 
#' @param seqnames_column The column number in bedgraphs which corresponds to the sequence names. Default is 1st column.  
#' @param start_column The column number in bedgraphs which corresponds to the start positions. Default is 2nd column. 
#' @param end_column The column number in bedgraphs which corresponds to the end positions. Default is 3rd column. 
#' @param value_column The column number in bedgraphs which corresponds to the methylation values. Default is 4th column. 
#' @param zero_based TRUE or FALSE indicating if files are zero-based. Default value is TRUE. 
#' @param normalization_factor An optional numerical value to divide methylation values by to convert them to fractions e.g. 100 if they are percentages. 
#' Default is not to leave values as they are in the input files. 
#' @param decimal_places Optional integer indicating the number of decimal places to round beta values to. Default is not to round. 
#' @param meth_sites A GRanges object with the locations of the methylation sites of interest. Any methylation sites in bedGraphs that are not in meth_sites are ignored.
#' @param sample_metadata Sample metadata to be used as colData for the RangedSummarizedExperiment.
#' @param hdf5_dir Directory to save HDF5 file. Is created if it doesn't exist. HDF5 file is called assays.h5. 
#' @param dataset_name Name to give data set in HDF5 file. Default is "beta".
#' @param overwrite TRUE or FALSE indicating whether to allow overwriting if dataset_name already exists in assays.h5. Default is FALSE.
#' @param chunkdim The dimensions of the chunks for the HDF5 file. Should be a vector of length 2 giving the number of rows and then the number of columns in each chunk.
#' Uses HDF5Array::getHDF5DumpChunkDim(length(meth_sites), length(bedgraphs))) by default. 
#' @param temporary_dir Name to give temporary directory created to store intermediate files. A directory with this name cannot already exist. 
#' Default is to create a name using tempfile("temporary_meth_chunks_"). 
#' Will be deleted after completion. 
#' @param BPPARAM A BiocParallelParam object for parallel processing. Defaults to `BiocParallel::bpparam()`. 
#' @param ... Additional arguments to be passed to HDF5Array::HDF5RealizationSink() for controlling the physical properties of the created HDF5 file, 
#' such as compression level. Uses the defaults for any properties that are not specified. 
#' @return A RangedSummarizedExperiment with methylation values for all methylation sites in meth_sites. methylation sites will be in the same order as sort(meth_sites). 
#' @export
#' @examples
#' 
#' # Load CpGs within first million base pairs of chromosome 1 as a GRanges object
#' data("hg38_cpgs_subset", package = "methodical")
#' 
#' # Get paths to bedGraphs
#' bedgraphs <- list.files(path = system.file('extdata', package = 'methodical'), 
#'   pattern = ".bg.gz", full.names = TRUE)
#' 
#' # Create sample metadata
#' sample_metadata <- data.frame(
#'   tcga_project = gsub("_.*", "", gsub("TCGA_", "", basename(bedgraphs))),
#'   sample_type = ifelse(grepl("N", basename(bedgraphs)), "Normal", "Tumour"),
#'   row.names = tools::file_path_sans_ext(basename(bedgraphs))
#' )
#' 
#' # Create a HDF5-backed RangedSummarizedExperiment from bedGraphs
#' meth_rse <- makeMethRSEFromBedgraphs(bedgraphs = bedgraphs, 
#'   meth_sites = hg38_cpgs_subset, sample_metadata = sample_metadata, 
#'   hdf5_dir = paste0(tempdir(), "/bedgraph_hdf5_1"))
#'   
makeMethRSEFromBedgraphs <- function(bedgraphs, 
  seqnames_column = 1, start_column = 2, end_column = 3, value_column = 4,
  zero_based = TRUE, normalization_factor = NULL, decimal_places = NA, 
  meth_sites, sample_metadata = NULL, hdf5_dir, dataset_name = "beta", overwrite = FALSE, chunkdim = NULL, 
  temporary_dir = NULL, BPPARAM = BiocParallel::bpparam(), ...){
  
  # Check that inputs have the correct data type
  stopifnot(is(bedgraphs, "character"), is(seqnames_column, "numeric") & seqnames_column >= 1,
    is(start_column, "numeric") & start_column >= 1, is(end_column, "numeric") & end_column >= 1,
    is(value_column, "numeric") & value_column >= 1, S4Vectors::isTRUEorFALSE(zero_based),
    is(normalization_factor, "numeric") | is.null(normalization_factor),
    is(decimal_places, "numeric") | is.na(decimal_places), is(meth_sites, "GRanges"),
    is(sample_metadata, "data.frame") | is.null(sample_metadata), is(hdf5_dir, "character"),
    is(dataset_name, "character"), S4Vectors::isTRUEorFALSE(overwrite), is(chunkdim, "numeric") | is.null(chunkdim),
    is(temporary_dir, "character") | is.null(temporary_dir), is(BPPARAM, "BiocParallelParam"))
    
  # Check that normalization_factor is a whole integer if provided
  if(!is.null(normalization_factor)){
    if(length(normalization_factor) != 1 | normalization_factor %% 1 != 0 | normalization_factor < 0){
      stop("normalization_factor should be single whole number")
    }
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
  
  # Check if meth_sites is sorted and print a message if it is not. 
  if(!S4Vectors::isSorted(meth_sites)){
    message("meth_sites is not sorted. It will be sorted and this sorted order will be used for methylation sites in the HDF5 file")
  }
  
  # Perform setup
  setup <- .make_meth_rse_setup(meth_files = bedgraphs, meth_sites = meth_sites, sample_metadata = sample_metadata, 
    hdf5_dir = hdf5_dir, dataset_name = dataset_name, overwrite = overwrite, chunkdim = chunkdim, 
    temporary_dir = temporary_dir, ...)
  
  # Read in bedGraphs and write data from chunks to appropriate temporary directory
  meth_sites_df <- .split_bedgraphs_into_chunks(bedgraphs = bedgraphs, 
    seqnames_column = seqnames_column, start_column = start_column, end_column = end_column, value_column = value_column,
    file_grid_columns = setup$file_grid_columns, meth_sites = meth_sites, meth_site_groups = setup$meth_site_groups, temp_chunk_dirs = setup$temp_chunk_dirs, 
    zero_based = zero_based, normalization_factor = normalization_factor, decimal_places = decimal_places, BPPARAM = BPPARAM)
  
  # Write the chunks to the HDF5 file
  .write_chunks_to_hdf5(temp_chunk_dirs = setup$temp_chunk_dirs, files_in_chunks = setup$files_in_chunks, 
    hdf5_sink = setup$hdf5_sink, hdf5_grid = setup$hdf5_grid)
  
  # Create a RangedSummarizedExperiment
  rse <- .create_meth_rse_from_hdf5(hdf5_filepath = setup$hdf5_filepath, hdf5_dir = hdf5_dir,
    meth_sites_df = meth_sites_df, sample_metadata = sample_metadata)
  
  # Delete temporary_dir if it is empty
  if(length(list.files(temporary_dir)) == 0){
    unlink(temporary_dir, recursive = TRUE)
  }
  
  return(rse)
  
}

#' Create a HDF5-backed RangedSummarizedExperiment for methylation values in array files
#'
#' @param array_files A vector of paths to bedGraph files. Automatically detects if array_files contain a header if every field in the first line is a character. 
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
#' @param BPPARAM A BiocParallelParam object for parallel processing. Defaults to `BiocParallel::bpparam()`. 
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
  overwrite = FALSE, chunkdim = NULL, temporary_dir = NULL, BPPARAM = BiocParallel::bpparam(), ...){
  
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
#' 
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