#' Create plot of methylation values
#'
#' @param cpg_values A data.frame with the methylation values of CpGs for different samples such as returned by extract_cpg_values_from_methrix. 
#' All CpG sites must be located on the same sequence.
#' @param reference_region An optional GRanges object with a single range. If provided, the x-axis will show the distance of CpG sites to the start of this region with CpGs upstream 
#' relative to the reference_region shown first. If not, the x-axis will show the start site coordinate of the CpG. 
#' @param samples_subset An optional vector of sample names to use to subset cpg_values.
#' @param geom geom to display. Can be one of geom_point, geom_line, geom_boxplot or geom_linepoint which indicates to use both geom_point and geom_line together.
#' @param sample_groups An optional factor vector giving the names of the groups that samples belong to. Alternatively, can set to "all_samples" to plot and colour each sample separately.
#' @param group_summary_function The name of a function to use to summarize values from each group before plotting. Can be any function which returns a single value from a numeric vector e.g. mean.
#' If not set, individual samples are plotted instead of group summaries.
#' @param group_difference An optional vector of length two with the names of two groups from sample_groups. Providing group_difference specifies to plot the difference between 
#' the summary values calculated for each group using group_summary_function by subtracting the second group from the first.
#' @param group_colors An optional vector of colours to be used for each group. Default is to use default ggplot2 colours. If group_difference is provided, should just be a single colour. 
#' @param title Title of the plot. Default is no title. 
#' @param xlabel Label for the X axis in the plot. Default is "Genomic Position".
#' @param ylabel Label for the Y axis in the plot. Default is "CpG Methylation".
#' @return A ggplot object of the mean methylation values of the CpGs in cpg_values. 
#' @export
plot_cpg_methylation = function(cpg_values, geom, samples_subset = NULL, sample_groups = NULL, group_summary_function = NULL, group_difference = NULL, group_colours = NULL, reference_region = NULL,
  title = NULL, xlabel =  "Genomic Position", ylabel = "CpG Methylation") {
  
  # Check that allowed values are provided for geom
  match.arg(arg = geom, choices = c("geom_point", "geom_line", "geom_linepoint", "geom_boxplot"), several.ok = F)
  
  # Make sure group_summary_function is a function
  if(!is.null(group_summary_function) & !is.function(group_summary_function)){
    stop("If provided, group_summary_function should be the unquoted name of a function which returns a single value for each row in a data.frame e.g. rowMeans")
  }
  
  # If samples_subset not provided, set to all samples in cpg_values
  if(is.null(samples_subset)){samples_subset = names(cpg_values)}
  
  # Subset cpg_values for specified samples if provided
  cpg_values = dplyr::select(cpg_values, dplyr::all_of(samples_subset))
  
  # Set sample_groups if they are not set
  if(!is.null(sample_groups)){
    if(all(sample_groups == "all_samples")){sample_groups = factor(samples_subset)} 
  } else {sample_groups = factor(rep("base", length(samples_subset)))}
  
  # Check that sample_groups is the same length as samples if they are provided
  if(!is.null(sample_groups) & length(samples_subset) != length(sample_groups)){stop("samples and sample_groups are not the same length")}
  if(!is.factor(sample_groups)){sample_groups = factor(sample_groups)}
  
  # Check if group_colours is equal to the number of groups unless using group_difference.
  if(!is.null(group_colours) & is.null(group_difference) & length(levels(sample_groups)) != length(group_colours)){stop("group_colours must provide one colour for each level in sample_groups")}
  
  # Check that group_difference is a vector of length two with the names of two groups in sample_groups if it is provided
  if(!is.null(group_difference)){
    if(length(group_difference) != 2 | any(!group_difference %in% sample_groups)){stop("group_difference should be a vector of length two with the names of two groups in sample_groups")}
  }
  
  # Create a vector matching samples to groups
  sample_to_group_vec = setNames(sample_groups, samples_subset)
  sample_to_group_list = split(samples_subset, sample_groups)
  
  # Check that all CpGs are on the same sequence
  if(length(seqlevels(GenomicRanges::GRanges(row.names(cpg_values)))) > 1){stop("All CpG sites should be located on the same sequence")}
  
  # Summarize values in groups using group_summary_function if specified
  if(!is.null(group_summary_function)){
    cpg_values = data.frame(lapply(sample_to_group_list, function(x) apply(dplyr::select(cpg_values, dplyr::any_of(x)), 1, function(x) group_summary_function(x, na.rm = T))))
    
    # If group_difference supplied, calculate the difference between the two provided groups
    if(!is.null(group_difference)){
      cpg_values = dplyr::transmute(cpg_values, difference = get(group_difference[1]) - get(group_difference[2]))
      sample_groups = rep("difference", length(samples_subset))
    }
  }
  
  # Add CpG chromosome_coordinate depending on whether reference_region provided
  if(!is.null(reference_region)){
    reference_region_start = resize(reference_region, 1, fix = "start")
    cpg_values$chromosome_coordinate = cpg_distances_to_region(reference_region_start, cpg_names = row.names(cpg_values))
  } else {cpg_values$chromosome_coordinate = GenomicRanges::start(GenomicRanges::GRanges(row.names(cpg_values)))}
  
  # Convert cpg_values to long format
  cpg_values = reshape2::melt(tibble::rownames_to_column(cpg_values, "position_name"), id.var = c("position_name", "chromosome_coordinate"))
  if(is.null(group_summary_function)){
    cpg_values$sample_group = sample_to_group_vec[cpg_values$variable]
  } else {cpg_values$sample_group =  cpg_values$variable}
  
  # Create a scatter plot of CpG Methylation values
  methyl_plot = ggplot(data = cpg_values, mapping = aes(x = chromosome_coordinate, y = value, color = sample_group, label = position_name)) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 12),
    axis.title = element_text(size = 20), axis.text = element_text(size = 18)) +
    scale_x_continuous(expand = c(0.005, 0.005), labels = scales::comma) +
    labs(x = xlabel, y = ylabel, title = title, color = NULL)
  
  if(!is.null(group_colours)){methyl_plot = methyl_plot + scale_color_manual(values = group_colours)}
  
  if(geom == "geom_boxplot"){methyl_plot = methyl_plot + geom_boxplot(aes(group = interaction(chromosome_coordinate, sample_group)))}
  if(geom == "geom_point"){methyl_plot = methyl_plot + geom_point(alpha = 0.75)}
  if(geom == "geom_line"){methyl_plot = methyl_plot + geom_line(alpha = 0.75)}
  if(geom == "geom_linepoint"){methyl_plot = methyl_plot + geom_point(alpha = 0.75) + geom_line(alpha = 0.75, show.legend = F)}
  
  # If there is only one sample_groups remove the legend
  if(length(unique(sample_groups)) == 1){methyl_plot = methyl_plot + theme(legend.position = "None")}
  
  return(methyl_plot)
}

#' Create a scatter plot of methylation and another feature
#'
#' Calculates the specified correlation coefficient of two named numeric vectors, matching elements by common names and ignoring elements whose name is not present in both vectors.
#'
#' @param cpg_values A data.frame with the methylation values of CpGs for different samples such as returned by extract_cpg_values_from_methrix. 
#' @param cpg_name The name of a CpG site in cpg_values.
#' @param feature A named numeric vector for a feature of interest.
#' @param sample_groups An optional factor giving the names of the groups that samples belong to. 
#' @param group_colors An optional vector of colours to be used for each group. Default is to use default ggplot2 colours. 
#' @param use_names A logical value indicating whether elements in cpg_values and feature should be matched by their names. Default is TRUE.
#' @param regression_line A logical value indicating whether to add a regression line to the plot. Default is FALSE. 
#' @param add_correlation A logical value indicating whether to add correlation to the plot. Default if FALSE. 
#' @param method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor().
#' @param position A numeric vector of length 2 giving the x and y coordinates on a 0 to 1 scale of the correlation label. 
#' Default is 90% of max methylation value for the x-coordinate and 0.99 for the y coordinate.
#' @param label_size Size of the add_correlation. Default is 3.
#' @param xlab x-axis label. Default is "Methylation Value".
#' @param ylab y-axis label. Default is  
#' @param title Title of the plot
#' @param alpha Value specifying the alpha
#' @return A ggplot object
#' @export
methylation_feature_scatter_plot = function(cpg_values, cpg_name, feature, sample_groups = NULL, group_colours = NULL, 
  use_names = T, regression_line = F, add_correlation = F, method = "p", position = NULL,
  label_size = 3, xlab = "Methylation Value", ylab = "Feature Value", title = "Methylation Correlation", alpha = 1){
  
  # Check that an allowed value is used for position
  if(!is.null(position)){
    if(!is.numeric(position) | length(position) != 2 | any(position < 0) | any(position > 1)){
      stop("position should be a numeric vector of length two with each value between 0 and 1")}}
  
  # Check that cpg_name is in cpg_values 
  if(length(cpg_name) != 1 || !cpg_name %in% row.names(cpg_values)){stop("cpg_name must be the name of a single cpg_site in cpg_values")}
  cpg_methylation = setNames(unlist(cpg_values[cpg_name, ]), names(cpg_values))
  
  # If position not provided set to 90% of max methylation value and 0.8
  if(is.null(position)){position = c(max(cpg_methylation, na.rm = T) * 0.95, 0.99)}
  
  # Create a scatterplot of two features, matching them by names (or optionally using a preficpg_methylation of names)
  if(use_names){
    common_names = intersect(names(cpg_methylation), names(feature))
    cpg_methylation = cpg_methylation[common_names]
    feature = feature[common_names]
  }
  
  # Get sample names
  samples = names(cpg_values)
  
  # Check that sample_colours and samples are the same length
  if(!is.null(sample_groups) & length(samples) != length(sample_groups)){stop("samples and sample_groups are not the same length")}
  
  # If sample_groups are NULL, set them as sample names. Then check if group_colours is equal to the number of groups.
  if(!is.null(sample_groups)){if(!is.factor(sample_groups)){sample_groups = factor(sample_groups)}}
  if(!is.null(group_colours) & length(levels(sample_groups)) != length(group_colours)){stop("group_colours must provide one colour for each level in sample_groups")}
  
  # Create a vector matching samples to groups
  if(!is.null(sample_groups)){sample_to_group_vec = setNames(sample_groups, samples)}
  
  feature_df = tibble::rownames_to_column(data.frame(cpg_methylation = cpg_methylation, feature = feature), "sample")
  if(!is.null(sample_groups)){feature_df$sample_groups = sample_to_group_vec[feature_df$sample]}
  scatter_plot = ggplot(feature_df, aes(x = cpg_methylation, y = feature, color = sample_groups, label = sample)) +
      geom_point(alpha = alpha) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 18)) +
        labs(x = xlab, y = ylab, title = title, colour = NULL)
  if(regression_line){scatter_plot = scatter_plot + geom_smooth(method = "lm", se = F, show.legend = F)}
  if(add_correlation){
    correlation_label = dplyr::case_when(method %in% c("pearson", "p") ~ "R", method %in% c("spearman", "s") ~ "rho",
      method %in% c("kendall", "k") ~ "tau")
    scatter_plot = scatter_plot + ggpubr::stat_cor(
      method = method, label.sep = "", p.digits = NA, label.x.npc = position[1], label.y.npc = position[2], size = label_size, show.legend = F, 
      cor.coef.name = correlation_label)
  }
  if(!is.null(group_colours)){scatter_plot = scatter_plot + scale_color_manual(values = group_colours)}
  
  return(scatter_plot)
}

