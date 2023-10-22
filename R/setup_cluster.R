#' Create an iterator function for use with bpiterate
#'
#' @param tss_for_chunk A list of GRanges with the TSS for the current chunk
#' @return A function which returns a list with the 
.tss_iterator = function(meth_values_chunk, tss_region_indices_list, transcript_values, tss_for_chunk){
  
  n <- length(tss_for_chunk)
  i <- 0L
  function() {
    if(i >= n)
      NULL
    else {
      i <<- i + 1L
      list(
        meth_table = meth_values_chunk[tss_region_indices_list[[i]], , drop = FALSE], 
        transcript_table = transcript_values[[i]],
        transcript_tss = tss_for_chunk[[i]], 
        transcript_name = names(tss_for_chunk)[i]
        )
    }
  }
}

#' Calculate meth site-transcript correlations for given TSS
#'
#' @param correlation_objects A list with a table of methylation values, 
#' expression values for transcripts, a GRanges for the TSS and the name of the transcript. 
#' @return A data.frame with the correlation values
.tss_correlations = function(correlation_objects){
      
  meth_table <- correlation_objects[["meth_table"]]
  transcript_table <- correlation_objects[["transcript_table"]]
  transcript_tss <- correlation_objects[["transcript_tss"]]
  transcript_name <- correlation_objects[["transcript_name"]]
      
    # Transpose meth_table
    meth_table <- t(meth_table)
    
    # Transpose transcript_table and name it with transcript name
    transcript_table <- setNames(data.frame(t(transcript_table)), transcript_name)
        
    tryCatch({
      transcript_meth_site_cors <- methodical::rapidCorTest(
        table1 = meth_table, table2 = transcript_table, 
        table1_name = "meth_site", table2_name = "transcript_name", 
        cor_method = cor_method, p_adjust_method = "none")
  
        # Add meth site distance to region if specified
        if(add_distance_to_region){transcript_meth_site_cors$distance_to_tss <- 
          methodical::strandedDistance(query_gr = GenomicRanges::GRanges(transcript_meth_site_cors$meth_site), 
            subject_gr = transcript_tss)
        }
        
      # Add transcript TSS as an attribute to correlation data.frame
      attributes(transcript_meth_site_cors)$tss_range <- transcript_tss
      
      # Save results to a RDS file in results_dir if it is provided
      if(!is.null(results_dir)){
          
        # Set name of RDS file as the transcript_name with .rds appended
        rds_file <- paste0(results_dir, "/", transcript_name, ".rds")
        
        # Write the results to the specified RDS file and return the RDS filepath
        saveRDS(transcript_meth_site_cors, rds_file)
        transcript_meth_site_cors <- rds_file
      }
        
      transcript_meth_site_cors
      
    }, error = function(x) data.frame())
    
}
