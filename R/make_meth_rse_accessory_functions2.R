#' Split data from a single methylation array files into chunks
#'
#' @param file Path to a methylation array files.
#' @param column The current grid column being processed. 
#' @param file_count The number of the file being processed
#' @param total_files The total number of files to be processed. 
#' @return A data.table with the probe sites sorted by seqnames, start and probe name.
.split_meth_array_file = function(file, column, file_count, total_files, probe_sites_df, normalization_factor){
    
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
  file_grid_columns, probe_ranges, probe_groups, temp_chunk_dirs, normalization_factor, decimal_places, ncores){
  
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
  
  # Loop through each chunk of array_files
  BiocParallel::bpmapply(.split_meth_array_file, file = array_files, column = file_grid_columns, 
    file_count = seq_along(array_files), MoreArgs = list(total_files = length(array_files)), BPPARAM = BPPARAM)
  
  # Run the garbage collection
  invisible(gc())
  
  # Return meth_sites_df
  return(probe_sites_df)
  
}

#' Split data from a single methylation array files into chunks
#'
#' @param bg_file Path to a bedgraph file.
#' @param column The current grid column being processed. 
#' @param file_count The number of the file being processed
#' @param total_files The total number of files to be processed. 
#' @return A data.table with the probe sites sorted by seqnames, start and probe name.
.split_bedgraph = function(file, column, file_count, total_files){
  
  # Set the current chunk to the first chunk of the current grid column
  current_chunk <- 1 + (column - 1) * length(meth_site_groups)
  
  # Print count of bedGraph being processed
  message(paste0("Processing bedGraph ", file_count, " out of ", length(bedgraphs), ": ", bg_file, "\n"))
  
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
  foreach::foreach(mg = meth_site_groups) %do% {
    
    # Subset meth_site_values for methylation sites in chunk
    meth_site_group_values <- data.table::as.data.table(meth_site_values[mg, "value"])
    
    # Write values to appropriate file
    data.table::fwrite(x = meth_site_group_values, 
      file = paste0(temp_chunk_dirs[current_chunk], "/", basename(bg_file)),
      row.names = FALSE, quote = FALSE, na = "NA", compress = "none", nThread = dt_threads)
    
    # Increase current chunk number
    current_chunk <- current_chunk + 1
    
  }
}