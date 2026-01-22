#' Check input files have the correct number of columns specified and that columns seem to be of correct type
#'
#' @param input_files A vector of input filepaths.
#' @param meth_files_columns A list specifying the columns in meth_files.
.check_input_files = function(input_files, meth_files_columns){
  
  # Loop through each input_file and check that its format is correct
  for(file in input_files){
  
    # Read in head of file
    file_head <- data.table::fread(file, nrows = 10, data.table = F)
    
    # Check that file has at least the required number of columns
    max_col_number = max(unlist(meth_files_columns))
    if(ncol(file_head) < max_col_number){
      stop(paste(file, "has only", ncol(file_head), "columns"))
    }
    
    # Check that columns are of specified type
    with(meth_files_columns, {
      if(!is(file_head[[start]], "numeric")) stop(paste("start in", file, "is not numeric"))
      if(!is.null(total_reads) && 
          !is(file_head[[total_reads]], "numeric")) stop(paste("total_reads in", file, "is not numeric"))
      if(!is.null(meth_reads) && 
          !is(file_head[[meth_reads]], "numeric")) stop(paste("meth_reads in", file, "is not numeric"))
      if(!is.null(unmeth_reads) && 
          !is(file_head[[unmeth_reads]], "numeric")) stop(paste("unmeth_reads in", file, "is not numeric"))
      if(!is.null(meth_fraction) && !is(file_head[[meth_fraction]], "numeric")) stop(paste("meth_fraction in", file, "is not numeric"))
    
    
      # Check that count columns all all integers
      count_columns = c(total_reads, meth_reads, unmeth_reads)
      if(sum(file_head[, count_columns] %% 1, na.rm = TRUE) > 0){
        stop(paste("total_reads, meth_reads and unmeth_reads columns for", file, "are not all integers"))
      }
    })
    
  }
  
  # Print message saying all checks have been passed
  message("All files passed initial checks")
  
}

#' Perform setup for makeMethRSEFromInputFiles or makeMethRSEFromArrayFiles
#'
#' @param meth_files A vector of paths to files with methylation values. 
#' Automatically detects if meth_files contain a header if every field in the first line is a character. 
#' @param meth_sites A GRanges object with the locations of the methylation sites of interest. Should contain separate ranges 
#' for each stand if meth_files are stranded (i.e. separate ranges for the C and G positions of CpG sites). 
#' Any positions in meth_files that are not in meth_sites are ignored. 
#' @param sample_metadata A data.frame with sample metadata to be used as colData for the RangedSummarizedExperiment.
#' @param hdf5_dir Directory to save HDF5 file. Is created if it doesn't exist. HDF5 file is called assays.h5. 
#' @param overwrite TRUE or FALSE indicating whether to allow overwriting if hdf5_dir already exists. 
#' @param chunkdim The dimensions of the chunks for the HDF5 file.
#' @param temporary_dir Name to give a temporary directory to store intermediate files. A directory with this name cannot already exist. 
#' @param ... Additional arguments to be passed to HDF5Array::HDF5RealizationSink. 
#' @return A list describing the setup to be used for makeMethRSEFromInputFiles or makeMethRSEFromArrayFiles.
.make_meth_rse_setup <- function(meth_files, meth_sites, sample_metadata, 
  hdf5_dir, overwrite, chunkdim, temporary_dir, ...){
  
  # If chunkdim not provided, use default values. Otherwise check that chunkdim is a numeric vector of length 2.
  if(is.null(chunkdim)){
    chunkdim <- HDF5Array::getHDF5DumpChunkDim(c(length(meth_sites), length(meth_files)))
  } else {
    if(!is.numeric(chunkdim) && length(chunkdim) != 2){
      stop("chunkdim must be a numeric vector of length 2 if provided")
    }
  }
  
  # Set chunk_rows and chunk_cols from chunkdim
  chunk_rows <- chunkdim[1]
  chunk_cols <- chunkdim[2]

  # If hdf5_dir doesn't exist, it is created
  if(!dir.exists(hdf5_dir)){
    dir.create(hdf5_dir)
  } 
  
  # Set hdf5_filepath as assays.h5 in hdf5_dir and if the file already exists remove it if overwrite is TRUE
  hdf5_filepath <- paste0(hdf5_dir, "/assays.h5")
  if(file.exists(hdf5_filepath) && overwrite){
    file.remove(hdf5_filepath)
  }
  
  # Check that the number of rows in sample_metadata equals the length of meth_files
  if(nrow(sample_metadata) != length(meth_files)){
    stop("Number of rows of sample_metadata must equal the number of methylation files")
  }
  
  # Create HDF5 realization sinks for methylation proportion and coverage named beta and Cov 
  beta_sink <- HDF5Array::HDF5RealizationSink(dim = as.integer(c(length(meth_sites), length(meth_files))), 
    filepath = hdf5_filepath, name = "beta", chunkdim = chunkdim, ...)
  Cov_sink <- HDF5Array::HDF5RealizationSink(dim = as.integer(c(length(meth_sites), length(meth_files))), 
    filepath = hdf5_filepath, name = "Cov", chunkdim = chunkdim, ...)
  
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
  setup <- list(hdf5_filepath = hdf5_filepath, beta_sink = beta_sink, Cov_sink = Cov_sink, 
    hdf5_grid = hdf5_grid, temp_chunk_dirs = temp_chunk_dirs, meth_site_groups = meth_site_groups, 
    file_grid_columns = file_grid_columns, files_in_chunks = files_in_chunks)
  
  return(setup)
  
}

#' Process a data.frame with methylation data so that it contains the correct columns
#'
#' @param meth_df A data.frame with methylation data.
#' @param zero_based TRUE or FALSE indicating if files are zero-based. 
.set_meth_df_columns = function(meth_df, zero_based){
  
  # Ensure seqnames is a character vector
  meth_df[["seqnames"]] <- as.character(meth_df[["seqnames"]]) 
  
  # Add 1 to start of regions if zero_based is TRUE
  if(zero_based){
    meth_df[["start"]] <- meth_df[["start"]] + 1
  }
  
  # Convert meth_fraction to a proportion if its appears to be a percentage
  if(!is.null(meth_df[["meth_fraction"]])){
    if(max(meth_df[["meth_fraction"]], na.rm = TRUE) > 1){
      message("meth_fraction appears to be percentages and so converting to proportions")
      meth_df[["meth_fraction"]] <- meth_df[["meth_fraction"]]/100
    }
  }
  
  # Add columns with meth_fraction and total_reads to meth_df, depending on which columns are present in meth_df and return meth_df
  if(!is.null(meth_df[["total_reads"]]) && !is.null(meth_df[["meth_fraction"]])){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads, meth_fraction)
  } else if(!is.null(meth_df[["total_reads"]]) && !is.null(meth_df[["meth_reads"]])){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads, meth_fraction = meth_reads/total_reads)
  } else if(!is.null(meth_df[["total_reads"]]) && !is.null(meth_df[["unmeth_reads"]])){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads, meth_fraction = 1 - unmeth_reads/total_reads)
  } else if(!is.null(meth_df[["meth_reads"]]) && !is.null(meth_df[["unmeth_reads"]])){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads = meth_reads + unmeth_reads, meth_fraction = meth_reads/total_reads)
  }
  
  return(meth_df)
  
}

#' Combine values for stranded data
#'
#' @param meth_df A data.frame with methylation data.
#' @param sequence_context A single character string or DNAString with the sequence context of the methylation sites e.g. CG or CHG.
.collapse_strands = function(meth_df, sequence_context){
  
  # Adjust start of sites on - strand so that they corresponds to start of sites on + strand
  meth_df[meth_df$strand == "-", ]$start <- meth_df[meth_df$strand == "-", ]$start - (nchar(sequence_context) - 1)
  
  # Combine counts from + and - strand, remove strand column and return
  meth_df_collapsed <- dplyr::summarise(dplyr::group_by(meth_df, seqnames, start),
    meth_fraction = sum(round(total_reads * meth_fraction))/sum(total_reads),
    total_reads = sum(total_reads), 
    .groups = "drop"
  )
  meth_df_collapsed$strand <- NULL
  return(data.table::as.data.table(meth_df_collapsed))
  
}

#' Split data from a single input methylation file into chunks
#'
#' @param meth_file Path to an input methylation file.
#' @param meth_files_columns A list specifying the columns in meth_files.
#' @param grid_column The current grid column being processed. 
#' @param file_count The number of the current file being processed.
#' @param parameters A list of parameters for processing meth_file.
#' @return Invisibly returns NULL. 
.split_meth_file <- function(meth_file, meth_files_columns, grid_column, file_count, parameters){
  
  # Attach the parameters locally using with
  with(parameters, {
  
    # Set the current chunk to the first chunk of the current grid column
    current_chunk <- 1 + ((grid_column - 1) * length(meth_site_groups))
    
    # Print count of meth_file being processed
    message(paste0("Processing file ", file_count, " out of ", total_files, ": ", meth_file, "\n"))
    
    # Initialize a data.frame for all methylation sites
    meth_site_values <- meth_sites_df
    
    # Remove NULL elements from meth_files_columns
    meth_files_columns = meth_files_columns[!sapply(meth_files_columns, is.null)]
    
    # Read in input methylation file with just columns in meth_files_columns
    meth_df <- data.table::fread(meth_file, nThread = dt_threads, 
      select = unname(unlist(meth_files_columns)), col.names = names(unlist(meth_files_columns)))
    
    # Check that no sites in meth_df overlap
    if(!GenomicRanges::isDisjoint(GenomicRanges::makeGRangesFromDataFrame(meth_df, 
      seqnames.field = "seqnames", start.field = "start", end.field = "start", 
      keep.extra.columns = F, starts.in.df.are.0based = FALSE))){
      stop(paste("There are overlapping sites in", meth_file))
    }
    
    # Adjust meth_df so that it has total_reads and meth_fraction column
    meth_df = .set_meth_df_columns(meth_df = meth_df, zero_based = zero_based)
    
    # If collapse_strands is TRUE add a column with strand to meth_df and combine values for strands
    if(collapse_strands){
      meth_df <- merge(meth_df, meth_sites_df, by = c("seqnames", "start"), all.x = FALSE, sort = FALSE)
      meth_df <- .collapse_strands(meth_df, sequence_context = sequence_context)
      meth_site_values <- dplyr::filter(meth_site_values, strand != "-")
      meth_site_values[["strand"]] <- "*"
    }
    
    # Round values if specified
    if(!is.na(decimal_places)){
      meth_df$meth_fraction <- round(meth_df$meth_fraction, decimal_places)
    }
    
    # Ensure seqlevels of meth_df are in the same order as meth_sites_df
    meth_df$seqnames <- factor(meth_df$seqnames, levels = levels(meth_sites_df$seqnames))
    
    # Set seqnames and start as keys for meth_df
    data.table::setkey(meth_df, seqnames, start)
    
    # Add values from meth_df to meth_site_values
    meth_site_values <- merge(meth_site_values, meth_df, by = c("seqnames", "start"), all.x = TRUE, sort = FALSE)
    
    # Remove meth_df and run the garbage collection
    rm(meth_df); invisible(gc())
      
    # Loop through each group of methylation sites
    `%do%` <- foreach::`%do%`
    foreach::foreach(mg = meth_site_groups) %do% {
      
      # Subset meth_site_values for methylation sites in chunk
      meth_site_group_values <- meth_site_values[mg, c("meth_fraction", "total_reads")]
      
      # Write values to appropriate file
      data.table::fwrite(x = meth_site_group_values, 
        file = paste0(temp_chunk_dirs[current_chunk], "/", basename(meth_file)),
        row.names = FALSE, quote = FALSE, na = "NA", compress = "none", nThread = dt_threads)
      
      # Increase current chunk number
      current_chunk <- current_chunk + 1
      
      # Return NULL
      return(invisible(NULL))
      
    }
  
  })
  
}

#' Split data from input methylation files into chunks
#'
#' @param meth_files Paths to input methylation files.
#' @param meth_files_columns A list specifying the columns in meth_files.
#' @param file_grid_columns The grid column number for each file. 
#' @param meth_sites_df A data.table with the positions of methylation sites.
#' @param collapse_strands TRUE or FALSE indicating whether or not to collapse data on + and - strands.  
#' @param sequence_context A single character string or DNAString with the sequence context of the methylation sites e.g. CG or CHG.
#' @param meth_site_groups A list with the indices of the methylation sites in each group. 
#' @param temp_chunk_dirs A vector giving the temporary directory associated with each chunk.
#' @param zero_based TRUE or FALSE indicating if files are zero-based. 
#' @param decimal_places Integer indicating the number of decimal places to round beta values to. 
#' @param BPPARAM A BiocParallelParam object. 
#' @return Invisibly returns NULL.
.split_meth_files_into_chunks <- function(meth_files, meth_files_columns, file_grid_columns, 
  meth_sites_df, collapse_strands, sequence_context, meth_site_groups, temp_chunk_dirs, zero_based, decimal_places, BPPARAM){
  
  # Set dt_threads to 1 if more than one core being used. 
  if(BiocParallel::bpnworkers(BPPARAM) > 1){
    dt_threads <- 1
  } else {
    dt_threads <- data.table::getDTthreads()
  }
  
  # Create a list with parameters to pass to .split_meth_file
  parameters_list <- list(total_files = length(meth_files), meth_site_groups = meth_site_groups,
    meth_sites_df = meth_sites_df, sequence_context = sequence_context, 
    collapse_strands = collapse_strands, meth_files_columns = meth_files_columns, dt_threads = dt_threads, 
    zero_based = zero_based, decimal_places = decimal_places, temp_chunk_dirs = temp_chunk_dirs)

  # Loop through each chunk of meth_files
  BiocParallel::bpmapply(.split_meth_file, meth_file = meth_files, grid_column = file_grid_columns, 
    file_count = seq_along(meth_files), MoreArgs = list(parameters = parameters_list), BPPARAM = BPPARAM)
  
  # Run the garbage collection
  invisible(gc())
  
  # Return NULL
  return(invisible(NULL))
  
}

#' Write chunks of data to a HDF5 sink
#'
#' @param temp_chunk_dirs A vector giving the temporary directory associated with each chunk.
#' @param files_in_chunks A list of files associated with each chunk in the order they should be placed.
#' @param beta_sink A HDF5RealizationSink for methylation proportions.
#' @param Cov_sink A HDF5RealizationSink for coverage.
#' @param hdf5_grid A RegularArrayGrid.
#' @return Invisibly returns TRUE. 
.write_chunks_to_hdf5 <- function(temp_chunk_dirs, files_in_chunks, beta_sink, Cov_sink, hdf5_grid){
  
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
    beta_chunk_data <- as.matrix(data.frame(lapply(files, data.table::fread, select = 1, colClasses = "numeric")))
    Cov_chunk_data <- as.matrix(data.frame(lapply(files, data.table::fread, select = 2, colClasses = "numeric")))
    invisible(gc())
    
    # Replace all NA values with 0
    Cov_chunk_data[is.na(Cov_chunk_data)] <- 0
    
    # Write M and Cov values to HDF5 file
    invisible(HDF5Array::write_block(block = beta_chunk_data, sink = beta_sink, viewport = hdf5_grid[[as.integer(chunk)]]))
    invisible(HDF5Array::write_block(block = Cov_chunk_data, sink = Cov_sink, viewport = hdf5_grid[[as.integer(chunk)]]))
    
    # Delete chunk temporary directory
    unlink(chunk_dir, recursive = TRUE)
      
  }
  
  # Remove chunk_data and run gc
  rm(beta_chunk_data, Cov_chunk_data); invisible(gc())
  
  # Invisibly return TRUE
  invisible(return(TRUE))
  
}

#' Create a RangedSummarizedExperiment for methylation values already deposited in HDF5
#'
#' @param hdf5_filepath Path to HDF5 file
#' @param meth_sites A sorted GRanges object with the locations of the methylation sites of interest.
#' @param sample_metadata A data.frame with sample metadata
#' @param hdf5_dir The path to the HDF5 directory. 
#' @return A RangedSummarizedExperiment with methylation values
.create_meth_rse_from_hdf5 <- function(hdf5_filepath, hdf5_dir, meth_sites, sample_metadata){
  
  # Create a list of data sets present in hdf5_filepath
  assay_list <- S4Vectors::SimpleList(setNames(lapply(c("beta", "Cov"), function(x) 
    HDF5Array::HDF5Array(filepath = hdf5_filepath, name = x)), c("beta", "Cov")))
  
  # Create a RangedSummarizedExperiment using the data sets in hdf5_dir, sample_metadata and meth_sites
  rse <- SummarizedExperiment::SummarizedExperiment(assays = assay_list, colData = sample_metadata, rowRanges = meth_sites)
  
  # Save rse in hdf5_dir
  HDF5Array:::.serialize_HDF5SummarizedExperiment(x = rse, rds_path = paste0(hdf5_dir, "/se.rds"), verbose = TRUE)
  
  # Return rse
  return(rse)
  
}

###

#' Split data from a single methylation array files into chunks
#'
#' @param file Path to a methylation array file.
#' @param column The current grid column being processed. 
#' @param file_count The number of the file being processed
#' @param parameters A list of parameters for processing the meth_file.
#' @return Invisibly returns NULL.
.split_meth_array_file <- function(file, column, file_count, parameters){
  
  # Attach the parameters
  attach(parameters)
    
  # Set the current chunk to the first chunk of the current grid column
  current_chunk <- 1 + (column - 1) * length(probe_groups)
  
  # Print count of meth_file being processed
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
#' @param decimal_places Integer indicating the number of decimal places to round beta values to. 
#' @param BPPARAM A BiocParallelParam object. 
#' @return A data.table with the probe sites sorted by seqnames, start and probe name.
.split_meth_array_files_into_chunks <- function(array_files, probe_name_column, beta_value_column, 
  file_grid_columns, probe_ranges, probe_groups, temp_chunk_dirs, decimal_places, BPPARAM){
  
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
    decimal_places = decimal_places, temp_chunk_dirs = temp_chunk_dirs)
  
  # Loop through each chunk of array_files
  BiocParallel::bpmapply(.split_meth_array_file, file = array_files, column = file_grid_columns, 
    file_count = seq_along(array_files), MoreArgs = list(parameter = parameters_list), BPPARAM = BPPARAM)
  
  # Run the garbage collection
  invisible(gc())
  
  # Return meth_sites_df
  return(probe_sites_df)
  
}
