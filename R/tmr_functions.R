#' Calculate methodical score and smooth it using a exponential weighted moving average
#' 
#' @param correlation_df A data.frame with correlation values between methylation sites and a transcript as returned by calculate_meth_site_transcript_cors. 
#' @param offset_length Number of methylation sites added upstream and downstream of a central methylation site to form a window, resulting in a window size of 2*offset_length + 1.
#' Default value is 10.
#' @param smoothing_factor Smoothing factor for exponential moving average. Should be a value between 0 and 1 with higher 
#' values resulting in a greater degree of smoothing. Default is 0.75. 
#' @return A GRanges object
#' @export
calculate_smoothed_methodical_scores = function(correlation_df, offset_length = 10, smoothing_factor = 0.75){
  
  # Check that smoothing_factor is between 0 and 1
  if(smoothing_factor < 0 | smoothing_factor > 1){stop("smoothing_factor should be between 0 and 1")}
  
  # Return an empty vector if there are less than 2*offset_length + 1 missing correlation values
  if(sum(!is.na(correlation_df$cor)) < (2*offset_length + 1)){
    return(numeric())
  }
  
  # Replace p-values of 0 with the next smallest p-value
  correlation_df$p_val = ifelse(correlation_df$p_val == 0, 
    min(correlation_df$p_val[correlation_df$p_val > 0]), correlation_df$p_val)
  
  # Calculate Methodical score which is the absolute log of the p_values multiplied by the inverse of the correlation sign
  correlation_df$score = sign(correlation_df$cor) * -log10(correlation_df$p_val)
  
  # Calculate exponentially decreasing weights
  wts = sapply(c(offset_length:0, 1:offset_length), function(x) smoothing_factor^x)
  
  # Calculate smoothed Methodical score with a weighted rolling mean
  smoothed_score = RcppRoll::roll_mean(x = correlation_df$score, n = (offset_length*2)+1, 
    weights = wts, normalize = T, na.rm = T, fill = NA, align = "center")
  
  # Return smoothed_score 
  return(smoothed_score)
  
}

#' Find TMRs where smoothed methodical scores exceed thresholds
#' 
#' @param meth_sites_gr A GRanges object with the location of methylation sites. 
#' @param smoothed_methodical_scores A numeric vector with the smoothed methodical scores associated with each methylation site. 
#' @param p_value_threshold The p_value cutoff to use. Default value is 0.005.
#' @param tss_gr An optional GRanges object giving the location of the TSS meth_sites_gr is associated with. 
#' @param transcript_id Name of the transcript associated with the TSS. 
#' @return A GRanges object with the location of TMRs. 
#' @export
test_tmrs = function(meth_sites_gr, smoothed_methodical_scores, p_value_threshold = 0.005, 
  tss_gr = NULL, transcript_id = NULL){
  
  # Create upper and lower bounds using the log of p_value_threshold
  upper_bound = -log10(p_value_threshold)
  lower_bound = log10(p_value_threshold)

  # Identify row indices of smoothed_methodical_score_df where smoothed Methodical score is outside the thresholds
  negative_tmr_meth_sites = which(smoothed_methodical_scores < lower_bound)
  positive_tmr_meth_sites = which(smoothed_methodical_scores > upper_bound)
  
  # Create IRanges from the row indices of negative and positive TMRs
  negative_tmr_meth_sites_ir = IRanges::IRanges(negative_tmr_meth_sites)
  positive_tmr_meth_sites_ir = IRanges::IRanges(positive_tmr_meth_sites)
  
  # Create data.frames from the IRanges
  negative_tmr_meth_sites_df = data.frame(negative_tmr_meth_sites_ir)
  positive_tmr_meth_sites_df = data.frame(positive_tmr_meth_sites_ir)
  
  # Combine negative_tmr_meth_sites_df and positive_tmr_meth_sites_df into a single data.frame and 
  # annotate if TMRs are negative or positive
  combined_tmr_df = dplyr::bind_rows(list(Negative = negative_tmr_meth_sites_df, 
    Positive = positive_tmr_meth_sites_df), .id = "direction")
  
  # Create a data.frame with the coordinates of methylation sites at the start and end of tmrs
  tmr_gr = GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(meth_sites_gr[combined_tmr_df$start]),
    ranges = IRanges::IRanges(
      start = GenomicRanges::start(meth_sites_gr[combined_tmr_df$start]), 
      end = GenomicRanges::start(meth_sites_gr[combined_tmr_df$end])
    )
  )
  
  # Return an empty GRanges if there are no TMRs
  if(length(tmr_gr) == 0){return(tmr_gr)}
  
  # Add direction to TMRs
  tmr_gr$direction = factor(combined_tmr_df$direction, levels = c("Negative", "Positive"))
  
  # Split tmr_gr into a list of negative and positive tmrs and return
  tmr_gr_list = split(tmr_gr, tmr_gr$direction)
  return(tmr_gr_list)
  
}

#' Calculate methodical score and smooth it using a exponential weighted moving average
#' 
#' @param correlation_df A data.frame with correlation values between methylation sites and a transcript 
#' or a path to an RDS file containing such a data.frame as returned by calculate_meth_site_transcript_cors. 
#' @param offset_length Number of methylation sites added upstream and downstream of a central methylation site to form a window, resulting in a window size of 2*offset_length + 1.
#' Default value is 10.
#' @param smoothing_factor Smoothing factor for exponential moving average. Should be a value between 0 and 1 with higher 
#' values resulting in a greater degree of smoothing. Default is 0.75. 
#' @param p_value_threshold The p_value cutoff to use. Default value is 0.005
#' @param min_gapwidth Merge TMRs with the same direction separated by less than this number of base pairs. default value is 150. 
#' @param min_meth_sites Minimum number of methylation sites that TMRs can contain. Default value is 5. 
#' @return A GRanges object with the location of TMRs.
#' @export
find_tmrs = function(correlation_df, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005, min_gapwidth = 150, min_meth_sites = 5){
  
  # If correlation_df is a character vector, try to read it as an RDS file and 
  # check that it is a data.frame with the correct columns
  if(is.character(correlation_df)){
    tryCatch({correlation_df = readRDS(correlation_df)}, error = function() 
      stop("Character string provided for correlation_df but could not find an RDS file at indicated file path"))
    if(!is.data.frame(correlation_df) | 
        !all(names(correlation_df) == c("meth_site", "transcript_name", "cor", "p_val", "distance_to_tss"))){
      stop("Object loaded from provided filepath is not a data.frame with the expected column names")
    } 
  }
  
  # If all p-values are above p_value_threshold or are NA can immediately return an empty GRanges
  if(sum(correlation_df$p_val > p_value_threshold, na.rm = T) == 0){return(GenomicRanges::GRanges())}
  
  # Calculate smoothed methodical scores 
  smoothed_methodical_scores = methodical:::calculate_smoothed_methodical_scores(
    correlation_df = correlation_df, offset_length = offset_length, smoothing_factor = smoothing_factor)
  
  # Create a GRanges with methylation sites from correlation_df
  meth_sites_gr = GenomicRanges::GRanges(correlation_df$meth_site)
  
  # Extract TSS range from correlation_df
  tss_gr = attributes(correlation_df)$tss_range
  transcript_id = tss_gr$transcript_id
  
  # Find TMRs where smoothed methodical scores exceed thresholds
  tmr_gr_list = methodical:::test_tmrs(meth_sites_gr = meth_sites_gr, smoothed_methodical_scores = smoothed_methodical_scores, 
    p_value_threshold = 0.005, tss_gr = tss_gr, transcript_id = transcript_id)
  
  # If there are no TMRs, return an empty GRanges
  if(sum(lengths(tmr_gr_list)) == 0){return(GenomicRanges::GRanges())}
  
  # Extract negative and positive TMRs from tmr_gr_list
  negative_tmrs = IRanges::append(tmr_gr_list$Negative, GenomicRanges::GRanges())
  positive_tmrs = IRanges::append(tmr_gr_list$Positive, GenomicRanges::GRanges())
  
  # Merge TMRs of the same direction in close proximity
  negative_tmrs_merged = GenomicRanges::reduce(negative_tmrs, ignore.strand = T, min.gapwidth = min_gapwidth)
  positive_tmrs_merged = GenomicRanges::reduce(positive_tmrs, ignore.strand = T, min.gapwidth = min_gapwidth)
  
  # Find negative and positive extensions
  negative_extensions = GenomicRanges::setdiff(negative_tmrs_merged, negative_tmrs)
  positive_extensions = GenomicRanges::setdiff(positive_tmrs_merged, positive_tmrs)
  
  # Remove extensions that overlap opposing TMRs
  negative_extensions = IRanges::subsetByOverlaps(negative_extensions, positive_tmrs, invert = T)
  positive_extensions = IRanges::subsetByOverlaps(positive_extensions, negative_tmrs, invert = T)
  
  # Combine the extensions with the original TMRs
  negative_tmrs = GenomicRanges::reduce(c(negative_tmrs, negative_extensions))
  positive_tmrs = GenomicRanges::reduce(c(positive_tmrs, positive_extensions))
  
  # Add direction
  if(length(negative_tmrs) > 0){negative_tmrs$direction = "Negative"}
  if(length(positive_tmrs) > 0){positive_tmrs$direction = "Positive"}
  
  # Turn into a single GRanges
  tmr_gr = c(negative_tmrs, positive_tmrs)
  
  # Sort TMRs depending on the strand of the associated TSS
  tmr_gr = sort(tmr_gr, ignore.strand = T, decreasing = as.character(strand(tss_gr)) == "-")
  
  # Number TMRs starting from most upstream to most downstream
  tmr_gr$tmr_name = trimws(paste(transcript_id, "tmr", 1:length(tmr_gr), sep = "_"), whitespace = "_")
  
  # Change direction to a factor
  tmr_gr$direction = factor(tmr_gr$direction, levels = c("Negative", "Positive"))
  
  # Add number of methylation sites that they overlap
  tmr_gr$meth_site_count = GenomicRanges::countOverlaps(tmr_gr, meth_sites_gr)
  
  # Add the distance to the TSS
  tmr_gr$distance_to_tss = methodical::stranded_distance(query_gr = tmr_gr, subject_gr = tss_gr)
  
  # Add the TSS location as a character vector
  tmr_gr$tss_location = as.character(tss_gr)
  
  # Add other metadata from the TSS and return TMRs
  mcols(tmr_gr) = cbind(mcols(tmr_gr), mcols(tss_gr))
  
  # Remove any TMRs with less than min_meth_sites
  tmr_gr = tmr_gr[tmr_gr$meth_site_count >= min_meth_sites]
  
  # If tmr_gr is now empty, return an empty GRanges without metadata or else return tmr_gr
  if(length(tmr_gr) == 0){return(GenomicRanges::GRanges())}
  return(tmr_gr)
  
}