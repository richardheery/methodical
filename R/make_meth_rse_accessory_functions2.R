.calculate_meth_fraction_and_total_reads = function(df, meth_fraction_col = NULL, total_reads_col = NULL, meth_reads_col = NULL, unmeth_reads_col = NULL){

  # Check that at least two of meth_reads, unmeth_reads, total_reads and meth_fraction are given
  if(sum(sapply(list(meth_reads_col, unmeth_reads_col, total_reads_col, meth_fraction_col), is.null)) > 2){
    stop("At least two of meth_fraction_col, total_reads_col, meth_reads_col and unmeth_reads_col must be provided")
  }
  
  # Convert meth_fraction to a proportion if its appears to be a percentage
  if(!is.null(meth_fraction_col)){
    if(max(df[[meth_fraction_col]], na.rm = TRUE) > 1)
      df[[meth_fraction_col]] = df[[meth_fraction_col]]/100
  }
  
  # Add columns with meth_fraction and total_reads to df, depending on which columns are present in df and return df
  if(is.null(meth_fraction_col) & !is.null(total_reads_col)){
    df$meth_fraction = df[[meth_fraction_col]]
    df$total_reads = df[[total_reads_col]]
  } else if(!is.null(meth_fraction_col) & !is.null(meth_reads_col)){
    df$meth_fraction = df[[meth_fraction_col]]
    df$total_reads = df[[meth_reads_col]]/df[[meth_fraction_col]]
  } else if(!is.null(meth_fraction_col) & !is.null(unmeth_reads_col)){
    df$meth_fraction = df[[meth_fraction_col]]
    df$total_reads = df[[unmeth_reads_col]]/(1 - df[[meth_fraction_col]])
  } else if(!is.null(total_reads_col) & !is.null(meth_reads_col)){
    df$total_reads = df[[total_reads_col]]
    df$meth_fraction = df[[meth_reads_col]]/df[[total_reads_col]]
  } else if(!is.null(total_reads_col) & !is.null(unmeth_reads_col)){
    df$total_reads = df[[total_reads_col]]
    df$meth_fraction = 1-df[[unmeth_reads_col]]/df[[total_reads_col]]
  } else if(!is.null(meth_reads_col) & !is.null(unmeth_reads_col)){
    df$total_reads = df[[meth_reads_col]] + df[[unmeth_reads_col]]
    df$meth_proportion = df[[meth_reads_col]] + df[[unmeth_reads_col]]
  }
  
  return(df)
  
}


