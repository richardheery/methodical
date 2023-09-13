#' Create plot of methylation values
#'
#' @param meth_site_values A data.frame with values associated with methylation sites. row.names should be the names of methylation sites and all methylation sites must be located on the same sequence. 
#' @param column_name Name of column in meth_site_values to plot. If not provided, will attempt to use the first column. 
#' @param reference_region An optional GRanges object with a single range. If provided, the x-axis will 
#' show the distance of methylation sites to the start of this region with methylation sites upstream 
#' relative to the reference_region shown first. If not, the x-axis will show the start site coordinate of the methylation site. 
#' @param title Title of the plot. Default is no title. 
#' @param xlabel Label for the X axis in the plot. Default is "Genomic Position".
#' @param ylabel Label for the Y axis in the plot. Default is "Methylation Site Values".
#' @param low_colour Colour to use for low values. Default value is "cadetblue4".
#' @param high_colour Colour to use for high values. Default value is "firebrick3".
#' @return A ggplot object 
#' @export
plot_meth_site_values = function(meth_site_values, column_name = NULL, reference_region = NULL, 
  title = NULL, xlabel = "Genomic Position", ylabel = "Methylation Site Values", low_colour = "cadetblue4", high_colour = "firebrick3") {
  
  # Check that if reference_region is provided, it has a length of 1
  if(!is.null(reference_region) & length(reference_region) > 1){stop("reference_region should have length of 1 if provided")}
  
  # Check that all methylation sites are on the same sequence
  if(length(seqlevels(GenomicRanges::GRanges(row.names(meth_site_values)))) > 1){
    stop("All methylation sites must be located on the same sequence")
  }
  
  # If column_name not provided, set to name of first column in meth_site_values
  if(is.null(column_name)){
    column_name = names(meth_site_values)[1]
  }
  
  # Check that column_name has length 1 
  if(length(column_name) > 1){stop("column_name should just be a character of length 1")}
  
  # Check that column_name is in the names of meth_site_values
  if(!column_name %in% names(meth_site_values)){stop(paste(column_name, "not the name of a column in meth_site_values"))}
  
  # Create a data.frame with the selected column
  plot_df = dplyr::select(meth_site_values, values = !!column_name)
  
  # Add meth_site_start position to plot_df
  plot_df$meth_site_start = GenomicRanges::start(GenomicRanges::GRanges(row.names(plot_df)))
  
  # Decide x-axis values for methylation sites depending on whether reference_region provided
  if(!is.null(reference_region)){
    plot_df$meth_site_plot_position = methodical::stranded_distance(query_gr = GRanges(row.names(plot_df)), subject_gr = reference_region)
  } else {
    plot_df$meth_site_plot_position = plot_df$meth_site_start 
  }
  
  # Subset plot_df for complete rows
  plot_df = plot_df[complete.cases(plot_df), ]
  
  # Create a scatter plot of methylation site values and return
  meth_site_plot = ggplot(data = plot_df, mapping = aes(x = meth_site_plot_position, y = values)) +
    geom_line(color = "black", alpha = 0.75) +
    geom_point(shape = 21, colour = "black", size = 4, alpha = 1, aes(fill = values)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 12),
      axis.title = element_text(size = 20), axis.text = element_text(size = 18), legend.position = "None") +
    scale_x_continuous(expand = c(0.005, 0.005), labels = scales::comma) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
    scale_fill_gradient2(low = low_colour, high = high_colour, mid = "white", midpoint = 0) +
    labs(x = xlabel, y = ylabel, title = title, color = NULL)
  
  return(meth_site_plot)
}

#' Create a plot with genomic annotation for a plot returned by plot_meth_site_values()
#'
#' Can combine the meth site values plot and genomic annotation together into a single plot or return the annotation plot separately. 
#'
#' @param meth_site_plot A plot of methylation site values (generally methylation level or correlation of methylation with transcription) around a TSS
#' @param reference_region An optional GRanges object with a single range. 
#' If provided, the x-axis will show the distance of methylation sites to the start of this region with methylation sites upstream 
#' relative to the reference_region shown first. If not, the x-axis will show the start site coordinate of the methylation site.
#' @param annotation_gr A GRanges object giving the locations of different classes of regions. 
#' There must be a metadata column named region_type giving the class of each region. 
#' If this is a factor, the levels will be used to order the region classes in the plot. 
#' @param class_colours An optional named vector of colours to use with different region classes. Names of vector match colours to region classes. 
#' @param annotation_line_size Linewidth for annotation plot. Default is 5. 
#' @param annotation_plot_height A value giving the proportion of the height of the plot devoted to the annotation. Default is 0.5. 
#' @param keep_meth_site_plot_legend A logical value indicating whether to retain the legend of meth_site_plot, if it has one. Default value is FALSE. 
#' @param annotation_plot_only A logical value indicating whether to return only the annotation plot. Default is to combine meth_site_plot with the annotation. 
#' @return A ggplot object
#' @export
annotate_meth_site_plot = function(meth_site_plot, annotation_gr, reference_region = NULL, region_class_colours = NULL, 
  annotation_line_size = 5, annotation_plot_height = 0.5, keep_meth_site_plot_legend = FALSE, annotation_plot_only = FALSE){
  
  # Check that if reference_region is provided, if has a length of 1
  if(!is.null(reference_region) & length(reference_region) > 1){stop("reference_region should have length of 1 if provided")}
  
  # Check that annotation_gr contains a metadata column called region_type
  if(is.null(annotation_gr$region_type)){
    stop("annotation_gr must contain a metadata column called region_type")
  }
  
  # Check that annotation_plot_height is between 0 and 1
  if(annotation_plot_height < 0 | annotation_plot_height > 1){
    stop("annotation_plot_height should be between 0 and 1")
  }
  
  # Get most extreme methylation sites in plot
  meth_site_min = row.names(meth_site_plot$data)[which.min(start(GRanges(row.names(meth_site_plot$data))))]
  meth_site_max = row.names(meth_site_plot$data)[which.max(start(GRanges(row.names(meth_site_plot$data))))]
  
  # Create a GRanges object which covers the plot
  plot_region = reduce(GRanges(c(meth_site_min, meth_site_max)), min.gapwidth = .Machine$integer.max)
  
  # Filter annotation regions for those which overlap plot_region
  annotation_gr = subsetByOverlaps(annotation_gr, plot_region)
  
  # Update start and end of annotation_gr so that they lie within plot_region
  start(annotation_gr) = pmax(start(annotation_gr), start(plot_region))
  end(annotation_gr) = pmin(end(annotation_gr), end(plot_region))
  
  # Decide x-axis values for methylation sites depending on whether reference_region provided
  if(!is.null(reference_region)){
    meth_site_plot$data$meth_site_plot_position = methodical::stranded_distance(query_gr = GRanges(row.names(meth_site_plot$data)), subject_gr = reference_region)
    annotation_gr = methodical::relative_ranges(genomic_regions = annotation_gr, reference_positions = reference_region)
  } else {
    meth_site_plot$data$meth_site_plot_position = meth_site_plot$data$meth_site_start 
  }
  
  # Convert annotation_gr to a data.frame
  annotation_df = data.frame(annotation_gr)
  
  # Ensure region_type is a factor and reverse their order so that they are displayed from top to bottom in the plot
  annotation_df$region_type = factor(annotation_df$region_type)
  annotation_df$region_type = factor(annotation_df$region_type, levels = rev(levels(annotation_df$region_type)))
  
  # Extract axis text size, axis title size and x-axis title from meth_site_plot
  axis_text_size = theme(meth_site_plot)[[1]]$theme$axis.text$size
  axis_title_size = theme(meth_site_plot)[[1]]$theme$axis.title$size
  x_axis_title = meth_site_plot$labels$x
  
  # Create colours for region classes
  if(is.null(region_class_colours)){
    palette = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", 
      "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
    region_class_colours = setNames(colorRampPalette(palette)(length(levels(annotation_df$region_type))), 
      levels(annotation_df$region_type))
  }
  
  # Create a linerange plot showing the positions of different genomic elements
  annotation_plot = ggplot(annotation_df, aes(xmin = start, xmax = end, x = NULL, y = region_type,  group = region_type, color = region_type)) + 
    geom_linerange(linewidth = annotation_line_size, linetype = "dashed", position = position_dodge(0.06)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 24),
      axis.title = element_text(size = axis_title_size), axis.text = element_text(size = axis_text_size ), legend.position = "None")  +
    labs(x = x_axis_title, y = "Region Classification") +
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
  meth_site_plot_legend = cowplot::get_legend(meth_site_plot)
  legends = cowplot::plot_grid(meth_site_plot_legend, NULL, nrow = 2, rel_heights = c(1 - annotation_plot_height, annotation_plot_height))
  
  # Combine meth_site_plot and annotation_plot
  annotated_meth_site_plot = cowplot::plot_grid(meth_site_plot + theme(legend.position = "none", 
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()), annotation_plot, 
    nrow = 2, align = "v", rel_heights = c(1 - annotation_plot_height, annotation_plot_height))
  
  # Add legend if specified
  if(keep_meth_site_plot_legend){
    annotated_meth_site_plot = cowplot::plot_grid(annotated_meth_site_plot, legends, rel_widths = c(1, 0.2))
  }
  
  # Return annotated_meth_site_plot
  return(annotated_meth_site_plot)
}

#' Add TMRs to a methylation site value plot
#'
#' @param meth_site_plot A plot of methylation site values around a TSS.
#' @param tmrs_gr A GRanges object giving the position of TMRs.
#' @param reference_region An optional GRanges object with a single range. If provided, the x-axis will show the distance of methylation sites to the start of this region with methylation sites upstream. 
#' relative to the reference_region shown first. If not, the x-axis will show the start site coordinate of the methylation site.
#' @param transcript_id An optional transcript ID. If provided, will attempt to filter tmrs_gr and reference_region using a metadata column called transcript_id with 
#' a value identical to the provided one. 
#' @param tmr_colours A vector with colours to use for negative and positive TMRs. Defaults to "#7B5C90" for negative and "#BFAB25" for positive TMRs. 
#' @param linewidth A numeric value to be provided as linewidth for geom_segment(). 
#' @return A ggplot object
#' @export
plot_tmrs = function(meth_site_plot, tmrs_gr, reference_region = NULL, transcript_id = NULL, tmr_colours = c("#A28CB1", "#D2C465"), linewidth = 5){
  
  # Filter tmrs_gr and reference_region for transcript_id if provided
  if(!is.null(transcript_id)){
    tmrs_gr = tmrs_gr[tmrs_gr$transcript_id == transcript_id]
    reference_region = reference_region[reference_region$transcript_id == transcript_id]
  }
  
  # Check that if reference_region is provided, it has a length of 1
  if(!is.null(reference_region) & length(reference_region) > 1){
    stop("reference_region should have length of 1 if provided")
  }
  
  # Decide positions for tmrs depending on whether reference_region provided
  if(!is.null(reference_region)){
      tmrs_df = data.frame(methodical::relative_ranges(
        genomic_regions = tmrs_gr, reference_positions = reference_region))
  } else {
      tmrs_df = data.frame(tmrs_gr)
  }
  
  
  # Add TMRs to meth_site_plot
  meth_site_plot_with_tmrs = meth_site_plot +
  geom_segment(data = tmrs_df, aes(x = start, xend = end, y = 0, yend = 0, color = direction), 
    linewidth = linewidth) + scale_color_manual(values = setNames(tmr_colours, levels(tmrs_df$direction))) + labs(color = "TMR Direction")
  
  # Return meth_site_plot_with_tmrs
  return(meth_site_plot_with_tmrs)
  
}

#' Create plot of Methodical score values for methylation sites around a TSS
#'
#' @param meth_site_values A data.frame with correlation values for methylation sites. There should be one column called "cor".
#' and another called "p_val" which are used to calculate the Methodical score. row.names should be the names of methylation sites and all methylation sites must be located on the same sequence. 
#' @param reference_region An optional GRanges object with a single range. If provided, the x-axis will show the distance of methylation sites to the start of this region with methylation sites upstream.
#' relative to the reference_region shown first. If not, the x-axis will show the start site coordinate of the methylation site. 
#' @param title Title of the plot. Default is no title. 
#' @param xlabel Label for the X axis in the plot. Default is "Genomic Position".
#' @param low_colour Colour to use for low values. Default value is "#7B5C90".
#' @param high_colour Colour to use for high values. Default value is "#BFAB25".
#' @param p_value_threshold The p-value threshold used to identify TMRs. Default value is 0.005. Set to NULL to not display thresholds. 
#' @return A ggplot object 
#' @export
plot_methodical_scores = function(meth_site_values, reference_region = NULL, 
  title = NULL, xlabel = "Genomic Position", low_colour = "#7B5C90", high_colour = "#BFAB25", p_value_threshold = 0.005) {
  
  # Check that if reference_region is provided, it has a length of 1
  if(!is.null(reference_region) & length(reference_region) > 1){stop("reference_region should have length of 1 if provided")}
  
  # Check that all methylation sites are on the same sequence
  if(length(seqlevels(GenomicRanges::GRanges(row.names(meth_site_values)))) > 1){
    stop("All methylation sites must be located on the same sequence")
  }
  
  # Add meth_site_start position to meth_site_values
  meth_site_values$meth_site_start = GenomicRanges::start(GenomicRanges::GRanges(row.names(meth_site_values)))
  
  # Decide x-axis values for methylation sites depending on whether reference_region provided
  if(!is.null(reference_region)){
    meth_site_values$meth_site_plot_position = methodical::stranded_distance(query_gr = GRanges(row.names(meth_site_values)), subject_gr = reference_region)
  } else {
    meth_site_values$meth_site_plot_position = meth_site_values$meth_site_start 
  }
  
  # Convert p-values into methodical score
  meth_site_values$methodical_score = log10(meth_site_values$p_val) * -sign(meth_site_values$cor)
  
  # Subset meth_site_values for necessary columns
  meth_site_values = dplyr::select(meth_site_values, meth_site_start, meth_site_plot_position, methodical_score, cor)
  
  # Create a scatter plot of methylation site values and return
  meth_site_plot = ggplot(data = meth_site_values, mapping = aes(x = meth_site_plot_position, y = methodical_score)) +
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
    meth_site_plot = meth_site_plot + 
      geom_hline(yintercept = log10(p_value_threshold), linetype = "dashed", colour = low_colour) +
      geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", colour = high_colour) 
  }
  
  return(meth_site_plot)
}

###

#' Add smoothed methodical scores curve to a methylation site plot
#'
#' @param meth_site_plot A plot of methylation site values around a TSS.
#' @param smoothed_methodical_scores A vector with smoothed methodical scores for each methylation site in meth_site_plot.
#' @param colour Colour of the smoothed curve. Default is "black"
#' @param linewidth Line width of the smoothed curve. Default value is 1.
#' @param alpha Alpha value for the curve. Default value is 0.75. 
#' @return A ggplot object 
#' @export
plot_smoothed_methodical_scores = function(meth_site_plot, smoothed_methodical_scores, colour = "black", linewidth = 1, alpha = 0.75) {
  
  # Add smoothed Methodical scores to meth_site_plot and return
  meth_site_plot_with_curve = meth_site_plot + 
    geom_line(mapping = aes(y = smoothed_methodical_scores), 
      color = colour, alpha = alpha, linewidth = linewidth)
  
  return(meth_site_plot_with_curve)
}
