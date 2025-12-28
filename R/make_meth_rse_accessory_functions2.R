.check_input_files = function(input_files, seqnames_column, start_column, 
  total_reads_col = NULL, meth_reads_col = NULL, unmeth_reads_col = NULL, meth_fraction_col = NULL){
  
  # Check that total_reads_col and at least one of meth_fraction_col, meth_reads_col or unmeth_reads_col or
  # else both meth_reads_col and unmeth_reads_col are provided 
  if((is.null(total_reads_col) | all(sapply(list(meth_fraction_col, meth_reads_col, unmeth_reads_col), is.null))) &&
    (is.null(meth_reads_col) | is.null(unmeth_reads_col))){
      stop("Either both meth_reads_col and unmeth_reads_col should be provided or else total_reads_col and one other column")
  }
  
  # Loop through each input_file and check that its format is correct
  for(file in input_files){
  
    # Read in head of file
    file_head <- data.table::fread(file, nrows = 10, data.table = F)
    
    # Check that file has at least the required number of columns
    max_col_number = max(c(seqnames_column, start_column, total_reads_col, meth_reads_col, unmeth_reads_col, meth_fraction_col))
    if(ncol(file_head) < max_col_number){
      stop(paste(file, "has only", ncol(file_head), "columns"))
    }
    
    # Check that columns are of specified type
    if(!is(file_head[[start_column]], "numeric")) stop(paste("start_column in", file, "is not numeric"))
    if(!is.null(total_reads_col) && !is(file_head[[total_reads_col]], "numeric")) stop(paste("total_reads_col in", file, "is not numeric"))
    if(!is.null(meth_reads_col) && !is(file_head[[meth_reads_col]], "numeric")) stop(paste("meth_reads_col in", file, "is not numeric"))
    if(!is.null(unmeth_reads_col) && !is(file_head[[unmeth_reads_col]], "numeric")) stop(paste("unmeth_reads_col in", file, "is not numeric"))
    if(!is.null(meth_fraction_col) && !is(file_head[[meth_fraction_col]], "numeric")) stop(paste("meth_fraction_col in", file, "is not numeric"))
    
    # Check that count columns all all integers
    if(sum(file_head[, c(total_reads_col, meth_reads_col, unmeth_reads_col)] %% 1, na.rm = TRUE) > 0){
      stop(paste("Columns specified by total_reads_col, meth_reads_col and unmeth_reads_col for", file, "are not all integers"))
    }
  }
  
  # Print message saying all checks have been passed
  message("All files passed initial checks")
  
}

.calculate_meth_fraction_and_total_reads = function(df, seqnames_column, start_column, 
  total_reads_col = NULL, meth_reads_col = NULL, unmeth_reads_col = NULL, meth_fraction_col = NULL){
  
  # Ensure seqnames_column and start_column are named seqnames and start
  names(df)[c(seqnames_column, start_column)] <- c("seqnames", "start")
  names(df)[c(total_reads_col, meth_reads_col, unmeth_reads_col, meth_fraction_col)] = 
    c("total_reads", "meth_reads", "unmeth_reads", "meth_fraction")[!sapply(list(total_reads_col, meth_reads_col, unmeth_reads_col, meth_fraction_col), is.null)]
  
  # Convert meth_fraction to a proportion if its appears to be a percentage
  if(!is.null(meth_fraction_col)){
    if(max(df[[meth_fraction_col]], na.rm = TRUE) > 1){
      message("meth_fraction appears to be percentages and so converting to fraction")
      df[[meth_fraction_col]] <- df[[meth_fraction_col]]/100
    }
  }
  
  # Add columns with meth_fraction and total_reads to df, depending on which columns are present in df and return df
  if(!is.null(total_reads_col) && !is.null(meth_fraction_col)){
    df = dplyr::transmute(df, seqnames, start, total_reads, meth_fraction)
  } else if(!is.null(total_reads_col) && !is.null(meth_reads_col)){
    df = dplyr::transmute(df, seqnames, start, total_reads, meth_fraction = total_reads/meth_reads)
  } else if(!is.null(total_reads_col) && !is.null(unmeth_reads_col)){
    df = dplyr::transmute(df, seqnames, start, total_reads, meth_fraction = 1 - total_reads/unmeth_reads)
  } else if(!is.null(meth_reads_col) && !is.null(unmeth_reads_col)){
    df = dplyr::transmute(df, seqnames, start, total_reads = meth_reads + unmeth_reads, meth_fraction = meth_reads/total_reads)
  }
  
  return(df)
  
}

.collapse_strands = function(df, meth_sites, starts.in.df.are.0based){
  
  # Set meth_site_width using the first site in meth_sites
  meth_site_width <- width(meth_sites)[1]
  
  # Separate meth_sites into sites on the + and - strand
  meth_sites_plus <- meth_sites[strand(meth_sites) == "+"]
  meth_sites_minus <- meth_sites[strand(meth_sites) == "-"]
  
  # Make a GRanges from df
  df_gr <- GenomicRanges::makeGRangesFromDataFrame(df, end.field = "start", 
    starts.in.df.are.0based = starts.in.df.are.0based)
  
  # Add strand to df and remove rows where strand is missing
  df$strand <- NA
  df$strand[df_gr %over% meth_sites_plus] <- "+"
  df$strand[df_gr %over% meth_sites_minus] <- "-"
  df <- dplyr::filter(df, !is.na(strand))
  
  # Adjust start of sites on - strand so that they corresponds to start of sites on + strand
  df[df$strand == "-", ]$start <- df[df$strand == "-", ]$start - meth_site_width
  
  # Combine counts from + and - strand and set strand as * and return
  df_collapsed <- dplyr::summarise(dplyr::group_by(df, seqnames, start),
    meth_fraction = sum(round(total_reads * meth_fraction))/sum(total_reads),
    total_reads = sum(total_reads)
  )
  df_collapsed$strand = "*"
  return(df_collapsed)
  
}