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
      if(!is(file_head[[start_column]], "numeric")) stop(paste("start_column in", file, "is not numeric"))
      if(!is.null(total_reads_col) && 
          !is(file_head[[total_reads_col]], "numeric")) stop(paste("total_reads_col in", file, "is not numeric"))
      if(!is.null(meth_reads_col) && 
          !is(file_head[[meth_reads_col]], "numeric")) stop(paste("meth_reads_col in", file, "is not numeric"))
      if(!is.null(unmeth_reads_col) && 
          !is(file_head[[unmeth_reads_col]], "numeric")) stop(paste("unmeth_reads_col in", file, "is not numeric"))
      if(!is.null(meth_fraction_col) && !is(file_head[[meth_fraction_col]], "numeric")) stop(paste("meth_fraction_col in", file, "is not numeric"))
    })
    
    # Check that count columns all all integers
    count_columns = c("total_reads_col", "meth_reads_col", "unmeth_reads_col")
    if(sum(file_head[, unlist(meth_files_columns[count_columns])] %% 1, na.rm = TRUE) > 0){
      stop(paste("Columns specified by total_reads_col, meth_reads_col and unmeth_reads_col for", file, "are not all integers"))
    }
  }
  
  # Print message saying all checks have been passed
  message("All files passed initial checks")
  
}

#' Process a data.frame with methylation data so that it contains the total number of reads and the fraction of methylated reads as columns
#'
#' @param meth_df A data.frame with methylation data.
#' @param meth_files_columns A list specifying the columns in meth_files.
.calculate_meth_fraction_and_total_reads = function(meth_df, meth_files_columns){
  
  # Ensure seqnames_column and start_column are named seqnames and start
  names(meth_df)[c(seqnames_column, start_column)] <- c("seqnames", "start")
  names(meth_df)[c(total_reads_col, meth_reads_col, unmeth_reads_col, meth_fraction_col)] = 
    c("total_reads", "meth_reads", "unmeth_reads", "meth_fraction")[!sapply(list(total_reads_col, meth_reads_col, unmeth_reads_col, meth_fraction_col), is.null)]
  
  # Convert meth_fraction to a proportion if its appears to be a percentage
  if(!is.null(meth_fraction_col)){
    if(max(meth_df[[meth_fraction_col]], na.rm = TRUE) > 1){
      message("meth_fraction appears to be percentages and so converting to proportions")
      meth_df[[meth_fraction_col]] <- meth_df[[meth_fraction_col]]/100
    }
  }
  
  # Add columns with meth_fraction and total_reads to meth_df, depending on which columns are present in meth_df and return meth_df
  if(!is.null(total_reads_col) && !is.null(meth_fraction_col)){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads, meth_fraction)
  } else if(!is.null(total_reads_col) && !is.null(meth_reads_col)){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads, meth_fraction = meth_reads/total_reads)
  } else if(!is.null(total_reads_col) && !is.null(unmeth_reads_col)){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads, meth_fraction = 1 - unmeth_reads/total_reads)
  } else if(!is.null(meth_reads_col) && !is.null(unmeth_reads_col)){
    meth_df = dplyr::transmute(meth_df, seqnames, start, total_reads = meth_reads + unmeth_reads, meth_fraction = meth_reads/total_reads)
  }
  
  return(meth_df)
  
}

#' Combine values for stranded data
#'
#' @param meth_df A data.frame with methylation data.
#' @param meth_site_width An integer giving the width of the methylation sites being studied e.g. 2 for CG sites. 
#' for each stand if meth_files are stranded (i.e. separate ranges for the C and G positions of CpG sites).
.collapse_strands = function(meth_df, meth_site_width){
  
  # Adjust start of sites on - strand so that they corresponds to start of sites on + strand
  meth_df[meth_df$strand == "-", ]$start <- meth_df[meth_df$strand == "-", ]$start - meth_site_width
  
  # Combine counts from + and - strand and set strand as * and return
  meth_df_collapsed <- dplyr::summarise(dplyr::group_by(meth_df, seqnames, start),
    meth_fraction = sum(round(total_reads * meth_fraction))/sum(total_reads),
    total_reads = sum(total_reads)
  )
  meth_df_collapsed$strand = "*"
  return(meth_df_collapsed)
  
}