#' Perform setup for makeMethRSEFromBedgraphs or makeMethRSEFromArrayFiles
#'
#' @param meth_files A vector of paths to files with methylation values. 
#' Automatically detects if meth_files contain a header if every field in the first line is a character. 
#' @param meth_sites A GRanges object with the locations of the methylation sites of interest. Any regions in meth_files that are not in meth_sites are ignored. 
#' @param sample_metadata Sample metadata to be used as colData for the RangedSummarizedExperiment.
#' @param hdf5_dir Directory to save HDF5 file. Is created if it doesn't exist. HDF5 file is called assays.h5. 
#' @param dataset_name Name to give data set in HDF5 file. 
#' @param overwrite TRUE or FALSE indicating whether to allow overwriting if dataset_name already exists in assays.h5. 
#' @param chunkdim The dimensions of the chunks for the HDF5 file.
#' @param temporary_dir Name to give a temporary directory to store intermediate files. A directory with this name cannot already exist. 
#' @param ... Additional arguments to be passed to HDF5Array::HDF5RealizationSink. 
#' @return A list describing the setup to be used for makeMethRSEFromBedgraphs or makeMethRSEFromArrayFiles.
.make_meth_rse_setup <- function(meth_files, meth_sites, sample_metadata, hdf5_dir, 
  dataset_name, overwrite, chunkdim, temporary_dir, ...){
  
  # If chunkdim not provided, use default values
  if(is.null(chunkdim)){
    chunkdim <- HDF5Array::getHDF5DumpChunkDim(c(length(meth_sites), length(meth_files)))
  } else {
    if(length(chunkdim) != 2){
      stop("chunkdim must be a numeric vector of length 2 if provided")
    }
  }
  
  # Set chunk_rows and chunk_cols from chunkdim
  chunk_rows <- chunkdim[1]
  chunk_cols <- chunkdim[2]

  # Define %do% from foreach
  `%do%` <- foreach::`%do%`

  # If hdf5_dir doesn't exist, it is created
  if(!dir.exists(hdf5_dir)){
    dir.create(hdf5_dir)
  }
  
  # Set hdf5_filepath as assays.h5 in hdf5_dir
  hdf5_filepath <- paste0(hdf5_dir, "/assays.h5")
  
  # Check if dataset_name is already present in HDF5 file and allow overwriting only if it is specified
  if(file.exists(hdf5_filepath)){
    if(dataset_name %in% rhdf5::h5ls(hdf5_filepath)$name & !overwrite){
      stop(paste("dataset named", dataset_name, "already present in HDF5 file and overwrite is set to FALSE"))
    } else if(dataset_name %in% rhdf5::h5ls(hdf5_filepath)$name & overwrite){
      rhdf5::h5delete(file = hdf5_filepath, dataset_name)
    }
  }
  
  # Create sample_metadata if it doesn't exist and if it does check that the number of rows equals the length of meth_files
  if(is.null(sample_metadata)){
    sample_metadata <- data.frame(
      row.names = tools::file_path_sans_ext(gsub("\\.gz$", "", basename(meth_files)))
    )
  } else {
    if(nrow(sample_metadata) != length(meth_files)){
      stop("Number of rows of sample_metadata must equal the number of methylation files")
    }
  }
  
  # Create a HDF5 realization sink
  hdf5_sink <- HDF5Array::HDF5RealizationSink(dim = as.integer(c(length(meth_sites), length(meth_files))), 
    filepath = hdf5_filepath, name = dataset_name, chunkdim = chunkdim, ...)
  
  # Make HDF5 grid
  hdf5_grid <- DelayedArray::RegularArrayGrid(
    refdim = as.integer(c(length(meth_sites), length(meth_files))), 
    spacings = chunkdim)
  
  # Create subdirectories in the temporary directory for each chunk
  temp_chunk_dirs <- sapply(seq_along(hdf5_grid), function(x) paste0(temporary_dir, "/chunk", x))
  lapply(temp_chunk_dirs, invisible(dir.create))
  
  # Get the grid row number associated with each meth site and then split meth_sites into groups
  # Each chunk will consist of the values for the methylation sites from one methylation site group for the meth_files in one file group
  meth_site_grid_rows <- ceiling(seq_along(meth_sites)/chunk_rows)
  meth_site_groups <- split(seq_along(meth_sites), meth_site_grid_rows)
  
  # Get the grid column number associated with each methylation file and then split meth_files into groups
  file_grid_columns <- ceiling(seq_along(meth_files)/chunk_cols)
  file_groups <- split(meth_files, file_grid_columns)
  
  # Create a list with the meth_files present in each chunk
  files_in_chunks <- rep(file_groups, each = length(meth_site_groups))

  # Create paths in each temporary directory to save data from the appropriate methylation file
  files_in_chunks <- lapply(seq_along(files_in_chunks), function(x)
    paste(temp_chunk_dirs[x], basename(files_in_chunks[[x]]), sep = "/"))
  
  # Create a list with all setup parameters and return
  setup_list <- list(hdf5_filepath = hdf5_filepath, hdf5_sink = hdf5_sink, hdf5_grid = hdf5_grid, temp_chunk_dirs = temp_chunk_dirs,
    meth_site_groups = meth_site_groups, file_grid_columns = file_grid_columns, files_in_chunks = files_in_chunks)
  
  return(setup_list)
  
}

#' Split data from a single methylation array files into chunks
#'
#' @param bg_file Path to a bedgraph file.
#' @param column The current grid column being processed. 
#' @param file_count The number of the current file being processed.
#' @param parameters A list of parameters for processing the bedgraph.
#' @return Invisibly returns NULL. 
.split_bedgraph <- function(bg_file, column, file_count, parameters){
  
  # Attach the parameters
  attach(parameters)
  
  # Set the current chunk to the first chunk of the current grid column
  current_chunk <- 1 + (column - 1) * length(meth_site_groups)
  
  # Print count of bedGraph being processed
  message(paste0("Processing bedGraph ", file_count, " out of ", total_files, ": ", bg_file, "\n"))
  
  # Initialize a data.frame for all methylation sites
  meth_site_values <- meth_sites_df
  
  # Read in bedGraph file
  bg <- setNames(data.table::fread(bg_file, 
    select = c(seqnames_column, start_column, end_column, value_column), nThread = dt_threads), 
    c("seqnames", "start", "end", "value"))
  
  # Add 1 to start of regions if zero_based is TRUE
  if(zero_based){
    bg$start <- bg$start + 1
  }
  
  # Convert values from percentages to proportions if specified
  if(!is.null(normalization_factor)){
    if(max(bg$value, na.rm = TRUE) > 1){
      bg$value <- bg$value/normalization_factor
    }
  }
  
  # Round values if specified
  if(!is.na(decimal_places)){
    bg$value <- round(bg$value, decimal_places)
  }
  
  # Ensure seqlevels of bg are in the same order as meth_sites_df
  bg$seqnames <- factor(bg$seqnames, levels = levels(meth_sites_df$seqnames))
  
  # Set seqnames and start as keys for bg
  data.table::setkey(bg, seqnames, start)
  
  # Add values from bg to meth_site_values
  meth_site_values <- merge(meth_site_values, bg, by = c("seqnames", "start"), all.x = TRUE, sort = FALSE)
  
  # Remove bg and run the garbage collection
  rm(bg); invisible(gc())
    
  # Loop through each chunk of methylation sites
  `%do%` <- foreach::`%do%`
  foreach::foreach(mg = meth_site_groups) %do% {
    
    # Subset meth_site_values for methylation sites in chunk
    meth_site_group_values <- data.table::as.data.table(meth_site_values[mg, "value"])
    
    # Write values to appropriate file
    data.table::fwrite(x = meth_site_group_values, 
      file = paste0(temp_chunk_dirs[current_chunk], "/", basename(bg_file)),
      row.names = FALSE, quote = FALSE, na = "NA", compress = "none", nThread = dt_threads)
    
    # Increase current chunk number
    current_chunk <- current_chunk + 1
    
    # Return NULL
    return(invisible(NULL))
    
  }
}

#' Split data from bedGraph files into chunks
#'
#' @param bedgraphs Paths to bedgraph files.
#' @param seqnames_column The column number in bedgraphs which corresponds to the sequence names. 
#' @param start_column The column number in bedgraphs which corresponds to the start positions. 
#' @param end_column The column number in bedgraphs which corresponds to the end positions. 
#' @param value_column The column number in bedgraphs which corresponds to the methylation values. 
#' @param file_grid_columns The grid column number for each file. 
#' @param meth_sites A GRanges object with the locations of the methylation sites of interest.
#' @param meth_site_groups A list with the indices of the methylation sites in each group. 
#' @param temp_chunk_dirs A vector giving the temporary directory associated with each chunk.
#' @param zero_based TRUE or FALSE indicating if files are zero-based. 
#' @param normalization_factor An optional numerical value to divide methylation values by to convert them to fractions e.g. 100 if they are percentages. 
#' Default is not to leave values as they are in the input files.  
#' @param decimal_places Integer indicating the number of decimal places to round beta values to. 
#' @param BPPARAM A BiocParallelParam object. 
#' @return A data.table with the methylation sites sorted by seqnames and start.
.split_bedgraphs_into_chunks <- function(bedgraphs, seqnames_column, start_column, end_column, value_column,
  file_grid_columns, meth_sites, meth_site_groups, temp_chunk_dirs, zero_based, normalization_factor, decimal_places, BPPARAM){
  
  # Set dt_threads to 1 if more than one core being used. 
  if(BiocParallel::bpnworkers(BPPARAM) > 1){
    dt_threads <- 1
  } else {
    dt_threads <- data.table::getDTthreads()
  }
  
  # Convert meth_sites into a data.table
  meth_sites_df <- data.table::data.table(data.frame(meth_sites)[seq_len(3)])
  
  # Set seqnames and start as keys for meth_sites_df
  data.table::setkey(meth_sites_df, seqnames, start)
  
  # Create a list with parameters to pass to .split_bedgraph
  parameters_list <- list(total_files = length(bedgraphs), meth_site_groups = meth_site_groups,
    meth_sites_df = meth_sites_df, seqnames_column = seqnames_column, start_column = start_column, 
    end_column = end_column, value_column = value_column, dt_threads = dt_threads, 
    zero_based = zero_based, normalization_factor = normalization_factor, 
    decimal_places = decimal_places, temp_chunk_dirs = temp_chunk_dirs)

  # Loop through each chunk of bedgraphs
  BiocParallel::bpmapply(.split_bedgraph, bg_file = bedgraphs, column = file_grid_columns, 
    file_count = seq_along(bedgraphs), MoreArgs = list(parameters = parameters_list), BPPARAM = BPPARAM)
  
  # Run the garbage collection
  invisible(gc())
  
  # Return meth_sites_df
  return(meth_sites_df)
  
}

#' Write chunks of data to a HDF5 sink
#'
#' @param temp_chunk_dirs A vector giving the temporary directory associated with each chunk.
#' @param files_in_chunks A list of files associated with each chunk in the order they should be placed.
#' @param hdf5_sink A HDF5RealizationSink.
#' @param hdf5_grid A RegularArrayGrid.
#' @return Invisibly returns TRUE. 
.write_chunks_to_hdf5 <- function(temp_chunk_dirs, files_in_chunks, hdf5_sink, hdf5_grid){
  
  # Define %do% from foreach
  `%do%` <- foreach::`%do%`
  
  # Loop through each chunk and write it the the HDF5 sink
  foreach::foreach(chunk = seq_along(temp_chunk_dirs)) %do% {
    
    # Print chunk being written
    message(paste0("Writing chunk ", chunk, " out of ", length(temp_chunk_dirs), "\n"))
    
    # Get chunk temporary directory
    chunk_dir <- temp_chunk_dirs[chunk]
    
    # Get the files associated with each chunk
    files <- files_in_chunks[[chunk]]
    
    # Read in all files in temporary directory as a data.frame of chunk data
    chunk_data <- as.matrix(data.frame(lapply(files, data.table::fread)))
    invisible(gc())
    
    # Write values to HDF5 file
    invisible(HDF5Array::write_block(block = chunk_data, sink = hdf5_sink, viewport = hdf5_grid[[as.integer(chunk)]]))
    
    # Delete chunk temporary directory
    unlink(chunk_dir, recursive = TRUE)
      
  }
  
  # Remove chunk_data and run gc
  rm(chunk_data); invisible(gc())
  
  # Invisibly return TRUE
  invisible(return(TRUE))
  
}

#' Create a RangedSummarizedExperiment for methylation values already deposited in HDF5
#'
#' @param hdf5_filepath Path to HDF5 file
#' @param meth_sites_df A data.frame with the positions of methylation sites
#' @param sample_metadata A data.frame with sample metadata
#' @param hdf5_dir The path to the HDF5 directory. 
#' @return A RangedSummarizedExperiment with methylation values
.create_meth_rse_from_hdf5 <- function(hdf5_filepath, hdf5_dir, meth_sites_df, sample_metadata){
  
  # Get the names of the assays in hdf5_filepath
  assay_names <- rhdf5::h5ls(hdf5_filepath)$name
  
  # Create a list of data sets present in hdf5_filepath
  assay_list <- S4Vectors::SimpleList(setNames(lapply(assay_names, function(x) 
    HDF5Array::HDF5Array(filepath = hdf5_filepath, name = x)), assay_names))
  
  # Update meth_sites to make sure they are in the same order as meth_sites_df
  meth_sites <- GenomicRanges::makeGRangesFromDataFrame(meth_sites_df, 
    keep.extra.columns = TRUE, seqinfo = levels(meth_sites_df$seqnames))
  
  # Create a RangedSummarizedExperiment using the data sets in hdf5_dir, sample_metadata and meth_sites
  rse <- SummarizedExperiment::SummarizedExperiment(assays = assay_list, colData = sample_metadata, rowRanges = meth_sites)
  
  # Save rse in hdf5_dir
  HDF5Array:::.serialize_HDF5SummarizedExperiment(x = rse, rds_path = paste0(hdf5_dir, "/se.rds"), verbose = TRUE)
  
  # Return rse
  return(rse)
  
}

#' Split data from a single methylation array files into chunks
#'
#' @param file Path to a methylation array file.
#' @param column The current grid column being processed. 
#' @param file_count The number of the file being processed
#' @param parameters A list of parameters for processing the bedgraph.
#' @return Invisibly returns NULL.
.split_meth_array_file <- function(file, column, file_count, parameters){
  
  # Attach the parameters
  attach(parameters)
    
  # Set the current chunk to the first chunk of the current grid column
  current_chunk <- 1 + (column - 1) * length(probe_groups)
  
  # Print count of bedGraph being processed
  message(paste0("Processing file ", file_count, " out of ", total_files, ": ", file, "\n"))
  
  # Initialize a data.frame for all probe sites
  probe_values <- probe_sites_df
  
  # Read in array file
  array_file <- setNames(data.table::fread(file, 
    select = c(probe_name_column, beta_value_column)), c("name", "value"))
  
  # Convert values from percentages to proportions if specified
  if(!is.null(normalization_factor)){
    if(max(array_file$value, na.rm = TRUE) > 1){
      array_file$value <- array_file$value/normalization_factor
    }
  }
  
  # Round values if specified
  if(!is.na(decimal_places)){
    array_file$value <- round(array_file$value, decimal_places)
  }
  
  # Ensure name of probes of array_file have levels with the same order as array_file
  array_file$name <- factor(array_file$name, levels = levels(probe_sites_df$name))
    
  # Set name as key for array_file
  data.table::setkey(array_file, name)
  
  # Add values from array_file to probe_values
  probe_values <- merge(probe_values, array_file, by = "name", all.x = TRUE, sort = FALSE)
  
  # Remove array_file and run the garbage collection
  rm(array_file); invisible(gc())
  
  # Loop through each chunk of methylation sites
  `%do%` <- foreach::`%do%`
  foreach::foreach(pg = probe_groups) %do% {
    
    # Subset probe_values for methylation sites in chunk
    probe_group_values <- data.table::as.data.table(probe_values[pg, "value"])
    
    # Write values to appropriate file
    data.table::fwrite(x = probe_group_values, 
      file = paste0(temp_chunk_dirs[current_chunk], "/", basename(file)),
      row.names = FALSE, quote = FALSE, na = "NA", compress = "none", nThread = dt_threads)
    
    # Increase current chunk number
    current_chunk <- current_chunk + 1
    
  }
  
  # Return TRUE
  return(TRUE)
  
}

#' Split data from methylation array files into chunks
#'
#' @param array_files Paths to methylation array files.
#' @param probe_name_column The column number in array_files which corresponds to the name of the probes. Default is 1st column. 
#' @param beta_value_column The column number in array_files which corresponds to the beta values. Default is 2nd column. 
#' @param file_grid_columns The grid column number for each file. 
#' @param probe_ranges A GRanges object giving the genomic locations of probes where each region corresponds to a separate probe.
#' @param probe_groups A list with the indices of the probes in each group.  
#' @param temp_chunk_dirs A vector giving the temporary directory associated with each chunk.
#' @param normalization_factor An optional numerical value to divide methylation values by to convert them to fractions e.g. 100 if they are percentages. 
#' Default is not to leave values as they are in the input files. 
#' @param decimal_places Integer indicating the number of decimal places to round beta values to. 
#' @param BPPARAM A BiocParallelParam object. 
#' @return A data.table with the probe sites sorted by seqnames, start and probe name.
.split_meth_array_files_into_chunks <- function(array_files, probe_name_column, beta_value_column, 
  file_grid_columns, probe_ranges, probe_groups, temp_chunk_dirs, normalization_factor, decimal_places, BPPARAM){
  
  # Set dt_threads to 1 if more than one core being used. 
  if(BiocParallel::bpnworkers(BPPARAM) > 1){
    dt_threads <- 1
  } else {
    dt_threads <- data.table::getDTthreads()
  }
  
  # Create a data.frame from probe_ranges
  probe_sites_df <- data.table::data.table(data.frame(probe_ranges)[c("seqnames", "start", "end", "name")])
  probe_sites_df$name <- factor(probe_sites_df$name, levels = probe_sites_df$name)
  
  # Set name as key for meth_sites_df
  data.table::setkey(probe_sites_df, seqnames, start, name)
  
  # Create a list with parameters to pass to .split_meth_array_file
  parameters_list <- list(total_files = length(array_files), probe_groups = probe_groups, 
    probe_sites_df = probe_sites_df, probe_name_column = probe_name_column, 
    beta_value_column = beta_value_column, dt_threads = dt_threads,
    normalization_factor = normalization_factor, decimal_places = decimal_places,
    temp_chunk_dirs = temp_chunk_dirs)
  
  # Loop through each chunk of array_files
  BiocParallel::bpmapply(.split_meth_array_file, file = array_files, column = file_grid_columns, 
    file_count = seq_along(array_files), MoreArgs = list(parameter = parameters_list), BPPARAM = BPPARAM)
  
  # Run the garbage collection
  invisible(gc())
  
  # Return meth_sites_df
  return(probe_sites_df)
  
}