#' Create plot of methylation values
#'
#' @param meth_site_values A data.frame with values associated with methylation sites. 
#' First column should be called meth_site and give the names of methylation sites. 
#' All methylation sites must be located on the same sequence. 
#' @param column_name Name of column in meth_site_values to plot. 
#' @param reference_tss TRUE or FALSE indicating whether to show distances on the X-axis
#' relative to the TSS stored as an attribute `tss_range` of meth_site_values. 
#' Alternatively, can provide a GRanges object with a single range for such a TSS site. 
#' In either case, will show the distance of methylation sites to the start of this region with methylation sites upstream 
#' relative to the reference_tss shown first. 
#' If FALSE (the default), the x-axis will instead show the start site coordinate of the methylation site. 
#' @param title Title of the plot. Default is no title. 
#' @param xlabel Label for the X axis in the plot. Defaults to "Distance to TSS" if reference_tss is used or
#' "seqname position" where seqname is the name of the relevant sequence.
#' @param ylabel Label for the Y axis in the plot. Default is "Value".
#' @param value_colours A vector with two colours to use, one for low values and the other for high values. 
#' Alternatively, can use one of two predefined colour sets by providing either "set1" or "set2":
#' set1 uses "#53868B" (blue) for low values and "#CD2626" (red) for high values 
#' while set2 uses "#7B5C90" (purple) for low values and ""#bfab25" (gold) for high values. Default is "set1". 
#' @return A ggplot object 
#' @examples 
#' # Load methylation-transcript correlation results for TUBB6 gene
#' data("tubb6_cpg_meth_transcript_cors", package = "methodical")
#' 
#' # Plot methylation-transcript correlation values around TUBB6 TSS
#' methodical::plotMethSiteValues(tubb6_cpg_meth_transcript_cors, column_name = "cor", ylabel = "Spearman Correlation")
#' 
#' # Create same plot but showing the distance to the TUBB6 TSS on the x-axis
#' methodical::plotMethSiteValues(tubb6_cpg_meth_transcript_cors, column_name = "cor", 
#'   ylabel = "Spearman Correlation", reference_tss = attributes(tubb6_cpg_meth_transcript_cors)$tss_range)
#' 
#' @export
plotMethSiteValues <- function(meth_site_values, column_name, reference_tss = FALSE, 
  title = NULL, xlabel = NULL, ylabel = "Value", value_colours = "set1"){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_site_values, "data.frame"), is(column_name, "character"),
    S4Vectors::isTRUEorFALSE(reference_tss) | is(reference_tss, "GRanges"), 
    is(title, "character") | is.null(title), is(xlabel, "character") | is.null(xlabel),
    is(ylabel, "character") | is.null(ylabel), is(value_colours, "character"))
  
  # Check that suitable input provided for value_colours and set low and high value colours if so
  if(length(value_colours) == 1){
    if(value_colours == "set1"){
      low_colour <- "#53868B"; high_colour <- "#CD2626"
    } else if(value_colours == "set2"){
      low_colour <- "#7B5C90"; high_colour <- "#bfab25"
    } else {
      stop("value_colours should be one of either set1 or set2 if only a single value provided")
    }
  } else {
    if(length(value_colours) == 2){
      low_colour <- value_colours[1]; high_colour <- value_colours[2]
    } else {
      stop("value_colours should be either a vector with two colours or else 
        one of either 'set1' or 'set2' if a single value is provided")
    }
  }
  
  # Change meth_site column to row names
  meth_site_values <- tibble::column_to_rownames(meth_site_values, "meth_site")
  
  # If reference_tss is TRUE, try to extract tss_range from meth_site_values
  if(is(reference_tss, "logical")){
    if(reference_tss){
      reference_tss <- attributes(meth_site_values)$tss_range 
      if(is.null(reference_tss)){
        stop("reference_tss was set to TRUE, but meth_site_values does not have an attribute called tss_range")
      }
    } else {
      reference_tss <- NULL
    }
  }
  
  # Check that reference_tss has a length of 1 if provided   
  if(!is.null(reference_tss) & (length(reference_tss) > 1 | !is(reference_tss, "GRanges"))){
    stop("GRanges indicated by reference_tss should have length of 1")
  }
  
  # Check that all methylation sites are on the same sequence
  if(length(unique(seqnames(GenomicRanges::GRanges(row.names(meth_site_values))))) > 1){
    stop("All methylation sites must be located on the same sequence")
  }
  
  # Check that column_name has length 1 
  if(length(column_name) > 1){stop("column_name should just be a character of length 1")}
  
  # Check that column_name is in the names of meth_site_values
  if(!column_name %in% names(meth_site_values)){stop(paste(column_name, "not the name of a column in meth_site_values"))}
  
  # Create a data.frame with the selected column
  plot_df <- dplyr::select(meth_site_values, values = !!column_name)
  
  # Add meth_site_start position to plot_df
  plot_df$meth_site_start <- GenomicRanges::start(GenomicRanges::GRanges(row.names(plot_df)))
  
  # Decide x-axis values for methylation sites depending on whether reference_tss provided
  if(!is.null(reference_tss)){
    plot_df$meth_site_plot_position <- methodical::strandedDistance(query_gr = GRanges(row.names(plot_df)), subject_gr = reference_tss)
  } else {
    plot_df$meth_site_plot_position <- plot_df$meth_site_start 
  }
  
  # Subset plot_df for complete rows
  plot_df <- plot_df[complete.cases(plot_df), ]
  
  # Create xlabel for plot if not provided
  if(is.null(xlabel)){
    if(!is.null(reference_tss)){
      xlabel <- "Distance to TSS"
    } else {
      xlabel <- paste(seqnames(GenomicRanges::GRanges(row.names(meth_site_values)))[1], "Position")
    }
  }
  
  # Create a scatter plot of Value and return
  meth_site_plot <- ggplot(data = plot_df, mapping = aes(x = meth_site_plot_position, y = values)) +
    geom_line(color = "black", alpha = 0.75) +
    geom_point(shape = 21, colour = "black", size = 4, alpha = 1, aes(fill = values)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 12),
      axis.title = element_text(size = 20), axis.text = element_text(size = 18), legend.position = "None") +
    scale_x_continuous(expand = c(0.005, 0.005), labels = scales::comma) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
    scale_fill_gradient2(low = low_colour, high = high_colour, mid = "white", midpoint = 0) +
    labs(x = xlabel, y = ylabel, title = title, color = NULL)
  
  # Add reference_tss as an attribute to plot if it was provided
  if(!is.null(reference_tss)){attributes(meth_site_plot)$tss_range <- reference_tss}
  
  return(meth_site_plot)
}

#' Create a plot with genomic annotation for a plot returned by plotMethSiteValues()
#'
#' Can combine the meth site values plot and genomic annotation together into a single plot or return the annotation plot separately. 
#'
#' @param meth_site_plot A plot of Value (generally methylation level or correlation of methylation with transcription) around a TSS
#' @param annotation_grl A GRangesList object (or list coercible to a GRangesList) where each component GRanges gives 
#' the locations of different classes of regions to display. Each class of region will 
#' be given a seprate colour in the plot, with regions ordered by the order of `names(annotation_grl)`. 
#' @param region_class_colours An optional vector of colours to use with different region classes. Must have same length as annotation_grl. 
#' @param reference_tss TRUE or FALSE indicating whether to show distances on the X-axis
#' relative to the TSS stored as an attribute `tss_range` of meth_site_plot. 
#' Alternatively, can provide a GRanges object with a single range for such a TSS site. 
#' In either case, will show the distance of methylation sites to the start of this region with methylation sites upstream 
#' relative to the reference_tss shown first. 
#' If FALSE (the default), the x-axis will instead show the start site coordinate of the methylation site. 
#' relative to the reference_tss shown first. If not, the x-axis will show the start site coordinate of the methylation site.
#' @param annotation_line_size Linewidth for annotation plot. Default is 5. 
#' @param annotation_plot_height A value giving the proportion of the height of the plot devoted to the annotation. Default is 0.5. 
#' @param keep_meth_site_plot_legend TRUE or FALSE indicating whether to retain the legend of meth_site_plot, if it has one. Default value is FALSE. 
#' @param annotation_plot_only TRUE or FALSE indicating whether to return only the annotation plot. Default is to combine meth_site_plot with the annotation. 
#' @return A ggplot object
#' @export
#' @examples 
#' # Get CpG islands from UCSC
#' cpg_island_annotation <- annotatr::build_annotations(genome = "hg38", annotations = "hg38_cpgs")
#' cpg_island_annotation <- GRangesList(split(cpg_island_annotation, cpg_island_annotation$type))
#'
#' # Load plot with CpG methylation correlation values for TUBB6
#' data("tubb6_correlation_plot", package = "methodical")
#' 
#' # Add positions of CpG islands to tubb6_correlation_plot
#' methodical::annotateMethSitePlot(tubb6_correlation_plot, annotation_grl = cpg_island_annotation, annotation_plot_height = 0.3)
#' 
annotateMethSitePlot <- function(meth_site_plot, annotation_grl, reference_tss = F, region_class_colours = NULL, 
  annotation_line_size = 5, annotation_plot_height = 0.5, keep_meth_site_plot_legend = FALSE, annotation_plot_only = FALSE){
  
  # If annotation_grl is a list, attempt to coerce it to a GRangesList
  if(is(annotation_grl, "list")){
    tryCatch(annotation_grl <- GRangesList(annotation_grl), 
      error = function(e) stop("annotation_grl is a list but cannot be coerced to a GRangesList"))
  }
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_site_plot, "ggplot"), is(annotation_grl, "GRangesList"),
    S4Vectors::isTRUEorFALSE(reference_tss) | is(reference_tss, "GRanges"), 
    is(region_class_colours, "character") | is.null(region_class_colours),
    is(annotation_line_size, "numeric"), is(annotation_line_size, "numeric"),
    is(annotation_plot_height, "numeric"), S4Vectors::isTRUEorFALSE(keep_meth_site_plot_legend),
    S4Vectors::isTRUEorFALSE(annotation_plot_only))
  
  # Create colours for region classes if region_class_colours not provided and 
  # check that it has the same length as annotation_grl if it is
  if(is.null(region_class_colours)){
    palette <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
      "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
    region_class_colours <- colorRampPalette(palette)(length(annotation_grl))
  } else if(length(annotation_grl) != length(region_class_colours)){
    stop("region_class_colours must have the same length as annotation_grl")
  }
    
  # If reference_tss is TRUE, try to extract tss_range from meth_site_plot
  if(is(reference_tss, "logical")){
    if(reference_tss){
      reference_tss <- attributes(meth_site_plot)$tss_range 
      if(is.null(reference_tss)){
        stop("reference_tss was set to TRUE, but meth_site_plot does not have an attribute called tss_range")
      }
    } else {
      reference_tss <- NULL
    }
  }
  
  # Check that reference_tss has a length of 1 if provided   
  if(!is.null(reference_tss) & (length(reference_tss) > 1 | !is(reference_tss, "GRanges"))){
    stop("GRanges indicated by reference_tss should have length of 1")
  }
  
  # Check that annotation_grl contains a metadata column called region_type
  if(is.null(names(annotation_grl))){
    message("names(annotation_grl) is NULL. Setting to genomic_regions_1, genomic_regions_2, etc.")
    names(annotation_grl) = paste0("genomic_regions_", seq_along(annotation_grl))
  }
  
  # Flatten annotation_grl
  annotation_grl = unlist(annotation_grl)
  annotation_grl$region_type = factor(names(annotation_grl), unique(names(annotation_grl)))
  
  # Check that annotation_plot_height is between 0 and 1
  if(annotation_plot_height < 0 | annotation_plot_height > 1){
    stop("annotation_plot_height should be between 0 and 1")
  }
  
  # Get most extreme methylation sites in plot
  meth_site_min <- row.names(meth_site_plot$data)[which.min(start(GRanges(row.names(meth_site_plot$data))))]
  meth_site_max <- row.names(meth_site_plot$data)[which.max(start(GRanges(row.names(meth_site_plot$data))))]
  
  # Create a GRanges object which covers the plot
  plot_region <- reduce(GRanges(c(meth_site_min, meth_site_max)), min.gapwidth = .Machine$integer.max)
  
  # Filter annotation regions for those which overlap plot_region
  annotation_grl <- subsetByOverlaps(annotation_grl, plot_region)
  
  # Update start and end of annotation_grl so that they lie within plot_region
  start(annotation_grl) <- pmax(start(annotation_grl), start(plot_region))
  end(annotation_grl) <- pmin(end(annotation_grl), end(plot_region))
  
  # Decide x-axis values for methylation sites depending on whether reference_tss provided
  if(!is.null(reference_tss)){
    meth_site_plot$data$meth_site_plot_position <- methodical::strandedDistance(query_gr = GRanges(row.names(meth_site_plot$data)), subject_gr = reference_tss)
    annotation_grl <- methodical::rangesRelativeToTSS(genomic_regions = annotation_grl, tss_gr = reference_tss)
  } else {
    meth_site_plot$data$meth_site_plot_position <- meth_site_plot$data$meth_site_start 
  }
  
  # Convert annotation_grl to a data.frame
  annotation_df <- data.frame(annotation_grl)
  
  # Ensure region_type is a factor and reverse their order so that they are displayed from top to bottom in the plot
  annotation_df$region_type <- factor(annotation_df$region_type)
  annotation_df$region_type <- factor(annotation_df$region_type, levels = rev(levels(annotation_df$region_type)))
  
  # Extract axis text size, axis title size and x-axis title from meth_site_plot
  axis_text_size <- theme(meth_site_plot)[[1]]$theme$axis.text$size
  axis_title_size <- theme(meth_site_plot)[[1]]$theme$axis.title$size
  x_axis_title <- meth_site_plot$labels$x
  
  # Create a linerange plot showing the positions of different genomic elements
  annotation_plot <- ggplot(annotation_df, aes(xmin = start, xmax = end, x = NULL, y = region_type,  group = region_type, color = region_type)) + 
    geom_linerange(linewidth = annotation_line_size, linetype = "dashed", position = position_dodge(0.06)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 24),
      axis.title = element_text(size = axis_title_size), axis.text = element_text(size = axis_text_size ), legend.position = "None")  +
    labs(x = x_axis_title, y = "Genome Annotation") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), labels = scales::comma, limits = ggplot_build(meth_site_plot)$layout$panel_params[[1]]$x.range) + 
    scale_color_manual(values = region_class_colours, guide = guide_legend(override.aes = list(color = "white"))) +
    # The following code makes the legend invisible
    theme(
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white")
    )
  
  # Return annotation_plot is return_only_annotation_plot is TRUE
  if(annotation_plot_only){return(annotation_plot)}
  
  # Get legend from meth_site_plot
  meth_site_plot_legend <- cowplot::get_legend(meth_site_plot)
  legends <- cowplot::plot_grid(meth_site_plot_legend, NULL, nrow = 2, rel_heights = c(1 - annotation_plot_height, annotation_plot_height))
  
  # Combine meth_site_plot and annotation_plot
  annotated_meth_site_plot <- cowplot::plot_grid(meth_site_plot + theme(legend.position = "none", 
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()), annotation_plot, 
    nrow = 2, align = "v", rel_heights = c(1 - annotation_plot_height, annotation_plot_height))
  
  # Add legend if specified
  if(keep_meth_site_plot_legend){
    annotated_meth_site_plot <- cowplot::plot_grid(annotated_meth_site_plot, legends, rel_widths = c(1, 0.2))
  }
  
  # Return annotated_meth_site_plot
  return(annotated_meth_site_plot)
}

#' Add TMRs to a methylation site value plot
#'
#' @param meth_site_plot A plot of Value around a TSS.
#' @param tmrs_gr A GRanges object giving the position of TMRs.
#' @param reference_tss An optional GRanges object with a single range. If provided, the x-axis will 
#' show the distance of methylation sites to the start of this region with methylation sites upstream
#' relative to the reference_tss shown first. If not, the x-axis will show the start site coordinate of the methylation site.
#' @param transcript_id An optional transcript ID. If provided, will attempt to filter tmrs_gr and reference_tss using a metadata column called transcript_id with 
#' a value identical to the provided one. 
#' @param tmr_colours A vector with colours to use for negative and positive TMRs. Defaults to "#7B5C90" for negative and "#BFAB25" for positive TMRs. 
#' @param linewidth A numeric value to be provided as linewidth for geom_segment(). 
#' @return A ggplot object
#' @export
#' @examples 
#' # Load methylation-transcript correlation results for TUBB6 gene
#' data("tubb6_cpg_meth_transcript_cors", package = "methodical")
#' 
#' # Plot methylation-transcript correlation values around TUBB6 TSS
#' tubb6_correlation_plot <- methodical::plotMethSiteValues(tubb6_cpg_meth_transcript_cors, column_name = "cor", ylabel = "Spearman Correlation")
#'   
#' # Find TMRs for TUBB6
#' tubb6_tmrs <- find_TMRs(correlation_df = tubb6_cpg_meth_transcript_cors)
#' 
#' # Plot TMRs on top of tubb6_correlation_plot
#' methodical::plotTMRs(tubb6_correlation_plot, tmrs_gr = tubb6_tmrs)
plotTMRs <- function(meth_site_plot, tmrs_gr, reference_tss = NULL, transcript_id = NULL, tmr_colours = c("#A28CB1", "#D2C465"), linewidth = 5){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_site_plot, "ggplot"), is(tmrs_gr, "GRanges"),
    is(reference_tss, "GRanges") | is.null(reference_tss), 
    is(transcript_id, "character") | is.null(transcript_id),
    is(tmr_colours, "character"), is(linewidth, "numeric"))
    
  # Filter tmrs_gr and reference_tss for transcript_id if provided
  if(!is.null(transcript_id)){
    tmrs_gr <- tmrs_gr[tmrs_gr$transcript_id == transcript_id]
    reference_tss <- reference_tss[reference_tss$transcript_id == transcript_id]
  }
  
  # Check that if reference_tss is provided, it has a length of 1
  if(!is.null(reference_tss) & length(reference_tss) > 1){
    stop("reference_tss should have length of 1 if provided")
  }
  
  # Decide positions for tmrs depending on whether reference_tss provided
  if(!is.null(reference_tss)){
      tmrs_df <- data.frame(methodical::rangesRelativeToTSS(
        genomic_regions = tmrs_gr, tss_gr = reference_tss))
  } else {
      tmrs_df <- data.frame(tmrs_gr)
  }
  
  # Make direction a factor with both Negative and Positive
  tmrs_df$direction = factor(tmrs_df$direction, levels = c("Negative", "Positive"))
  
  # Add TMRs to meth_site_plot
  meth_site_plot_with_tmrs <- meth_site_plot +
  geom_segment(data = tmrs_df, aes(x = start, xend = end, y = 0, yend = 0, color = direction), 
    linewidth = linewidth) + 
  scale_color_manual(values = setNames(tmr_colours, levels(tmrs_df$direction)), drop = FALSE) + 
  labs(color = "TMR Direction")
  
  # Return meth_site_plot_with_tmrs
  return(meth_site_plot_with_tmrs)
  
}

#' Create plot of Methodical score values for methylation sites around a TSS
#'
#' @param meth_site_values A data.frame with correlation values for methylation sites. There should be one column called "cor".
#' and another called "p_val" which are used to calculate the Methodical score. row.names should be the names of methylation sites and all methylation sites must be located on the same sequence. 
#' @param reference_tss An optional GRanges object with a single range. If provided, the x-axis will show the distance of methylation sites to the start of this region with methylation sites upstream.
#' relative to the reference_tss shown first. If not, the x-axis will show the start site coordinate of the methylation site. 
#' @param p_value_threshold The p-value threshold used to identify TMRs. Default value is 0.005. Set to NULL to turn off significance thresholds.
#' @param smooth_scores TRUE or FALSE indicating whether to display a curve of smoothed Methodical scores on top of the plot. Default is TRUE.
#' @param offset_length Offset length to be supplied to calculateSmoothedMethodicalScores.
#' @param smoothing_factor Smoothing factor to be provided to calculateSmoothedMethodicalScores.
#' @param smoothed_curve_colour Colour of the smoothed curve. Default is "black".
#' @param linewidth Line width of the smoothed curve. Default value is 1.
#' @param curve_alpha Alpha value for the curve. Default value is 0.75. 
#' @param title Title of the plot. Default is no title. 
#' @param xlabel Label for the X axis in the plot. Default is "Genomic Position".
#' @param low_colour Colour to use for low values. Default value is "#7B5C90".
#' @param high_colour Colour to use for high values. Default value is "#BFAB25".
#' @return A ggplot object 
#' @export
#' @examples 
#' # Load methylation-transcript correlation results for TUBB6 gene
#' data("tubb6_cpg_meth_transcript_cors", package = "methodical")
#'   
#' # Calculate and plot Methodical scores from correlation values
#' methodical::plotMethodicalScores(tubb6_cpg_meth_transcript_cors, reference_tss = attributes(tubb6_cpg_meth_transcript_cors)$tss_range)
plotMethodicalScores <- function(meth_site_values, reference_tss = NULL, p_value_threshold = 0.005,
  smooth_scores = TRUE, offset_length = 10, smoothing_factor = 0.75, 
  smoothed_curve_colour = "black", linewidth = 1, curve_alpha = 0.75, 
  title = NULL, xlabel = "Genomic Position", low_colour = "#7B5C90", high_colour = "#BFAB25"){
  
  # Check that inputs have the correct data type
  stopifnot(is(meth_site_values, "data.frame"), 
    is(reference_tss, "GRanges") | is.null(reference_tss), 
    is(p_value_threshold, "numeric") | is.null(p_value_threshold), S4Vectors::isTRUEorFALSE(smooth_scores),
    is(offset_length, "numeric"), is(smoothing_factor, "numeric"),
    is(smoothed_curve_colour, "character"), is(linewidth, "numeric"),
    is(curve_alpha, "numeric"), is(title, "character") | is.null(title),
    is(xlabel, "character") | is.null(xlabel), is(low_colour, "character"), 
    is(high_colour, "character"))
  
  # Change meth_site column to row names
  meth_site_values_plot_df <- tibble::column_to_rownames(meth_site_values, "meth_site")
  
  # Check that if reference_tss is provided, it has a length of 1
  if(!is.null(reference_tss) & length(reference_tss) > 1){stop("reference_tss should have length of 1 if provided")}
  
  # Check that all methylation sites are on the same sequence
  if(length(seqlevels(GenomicRanges::GRanges(row.names(meth_site_values_plot_df)))) > 1){
    stop("All methylation sites must be located on the same sequence")
  }
  
  # Add meth_site_start position to meth_site_values_plot_df
  meth_site_values_plot_df$meth_site_start <- GenomicRanges::start(GenomicRanges::GRanges(row.names(meth_site_values_plot_df)))
  
  # Decide x-axis values for methylation sites depending on whether reference_tss provided
  if(!is.null(reference_tss)){
    meth_site_values_plot_df$meth_site_plot_position <- methodical::strandedDistance(query_gr = GRanges(row.names(meth_site_values_plot_df)), subject_gr = reference_tss)
  } else {
    meth_site_values_plot_df$meth_site_plot_position <- meth_site_values_plot_df$meth_site_start 
  }
  
  # Convert p-values into methodical score
  meth_site_values_plot_df$methodical_score <- log10(meth_site_values_plot_df$p_val) * -sign(meth_site_values_plot_df$cor)
  
  # Subset meth_site_values_plot_df for necessary columns
  meth_site_values_plot_df <- dplyr::select(meth_site_values_plot_df, meth_site_start, meth_site_plot_position, methodical_score, cor)
  
  # Create a scatter plot of Value and return
  meth_site_plot <- ggplot(data = meth_site_values_plot_df, mapping = aes(x = meth_site_plot_position, y = methodical_score)) +
    geom_line(color = "black", alpha = 0.75) +
    geom_point(shape = 21, colour = "black", size = 4, alpha = 1, aes(fill = cor)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 12),
      axis.title = element_text(size = 20), axis.text = element_text(size = 18), legend.position = "None") +
    scale_x_continuous(expand = c(0.005, 0.005), labels = scales::comma) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
    scale_fill_gradient2(low = low_colour, high = high_colour, mid = "white", midpoint = 0) +
    labs(x = xlabel, y = "Methodical Score", title = title, color = NULL) 
  
  # Add TMR thresholds if specified
  if(!is.null(p_value_threshold)){
    meth_site_plot <- meth_site_plot + 
      geom_hline(yintercept = log10(p_value_threshold), linetype = "dashed", colour = low_colour) +
      geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", colour = high_colour) 
  }
  
  # Add smoothed Methodical scores if specified
  if(smooth_scores){
    smoothed_methodical_scores <- calculateSmoothedMethodicalScores(correlation_df = meth_site_values)
    meth_site_plot <- meth_site_plot +
    geom_line(mapping = aes(y = smoothed_methodical_scores), 
      color = smoothed_curve_colour, alpha = curve_alpha, linewidth = linewidth)
  }
  
  return(meth_site_plot)
}
